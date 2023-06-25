function analytic_force_per_unit_length(coil::CoilCircularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    I = coil.current
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * I * I / (4π * R) * (log(8 * R / a) - 0.75)
end

function analytic_force_per_unit_length(coil::CoilRectangularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    I = coil.current
    R = coil.curve.R0
    a = coil.a
    b = coil.b
    thirteen_over_twelve = 1.0833333333333333
    return μ0 * I * I / (4π * R) * (
        log(8 * R / sqrt(a * b)) 
        + thirteen_over_twelve - 0.5 * rectangular_xsection_k(a, b)
    )
end

# Data are from
# circular_coil_high_fidelity_over_analytic_force_rtol_1.0e-7_atol_1.0e-7_2023-02-17T02:57:13.639.dat
# Columns are:
# a/R, (high fidelity force)/(analytic force), time
const HIFI_FORCE_DATA = [
    0.0031622776601683794  0.999998730739337  1096.289042417
    0.003651741272548377  0.9999983152945607  1127.141660709
    0.004216965034285822  0.9999977482671788  1057.835328666
    0.004869675251658631  0.9999969959001087  1056.93196675
    0.005623413251903491  0.9999959937707692  999.009765791
    0.006493816315762113  0.9999946822554593  1040.893010542
    0.007498942093324558  0.9999928764404703  796.6792955
    0.008659643233600654  0.9999904993074297  708.913499792
    0.010000000000000002  0.9999873239392976  651.006336333
    0.011547819846894581  0.9999830941557939  628.86885525
    0.01333521432163324  0.9999774442518211  577.512847792
    0.01539926526059492  0.9999699087927443  552.565091208
    0.01778279410038923  0.9999598567839436  555.340665
    0.02053525026457146  0.99994644292227  623.065368208
    0.023713737056616554  0.9999285507284997  634.579416334
    0.027384196342643614  0.999904675957565  682.98895125
    0.03162277660168379  0.9998728183361064  651.714489375
    0.03651741272548377  0.9998303094927741  600.891050083
    0.042169650342858224  0.9997735811769564  569.889291
    0.04869675251658631  0.9996978793628768  522.663210375
    0.05623413251903491  0.9995968484241344  509.849693458
    0.06493816315762113  0.9994620027393627  491.207724375
    0.07498942093324558  0.9992820127733213  520.81247425
    0.08659643233600653  0.9990417426616669  502.397032791
    0.1  0.9987209678389167  468.94140175
    0.11547819846894582  0.9982926654765886  476.875976625
    0.1333521432163324  0.9977207091250009  405.783202667
    0.1539926526059492  0.9969567867489528  449.556182208
    0.1778279410038923  0.9959362763144891  461.15434325
    0.2053525026457146  0.9945726841820015  441.427283792
    0.23713737056616552  0.9927501881983675  415.201804042
    0.27384196342643613  0.9903136131816459  441.330066292
    0.31622776601683794  0.9870549175394838  462.2300085
    0.3651741272548377  0.9826950255515474  418.171884958
    0.4216965034285822  0.9768593903757312  400.739608167
    0.4869675251658631  0.9690452485822869  414.374111333
    0.5623413251903491  0.958578122966602  434.020034667
    0.6493816315762113  0.9445551128266463  429.259728875
    0.7498942093324559  0.9257743163026707  446.194293334
    0.8659643233600653  0.9006577658956955  484.066217833
    1.0  0.8672096485513339  534.621137417
]

"""
For a circular coil, compute the force per unit length by interpolating
pretabulated data from high-fidelity calculations.
"""
function interpolated_force_per_unit_length(coil::CoilCircularXSection)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    log10_a_over_R = (-2.5):0.0625:0
    interpolant = cubic_spline_interpolation(
        log10_a_over_R, 
        log.(1 ./ HIFI_FORCE_DATA[:, 2] .- 1), 
        extrapolation_bc=Line()
    )
    analytic_over_hifi_minus_1 = exp(interpolant(log10(coil.aminor / coil.curve.R0)))
    return analytic_force_per_unit_length(coil) / (analytic_over_hifi_minus_1 + 1)
end

"""
Compute the self-force using the regularized 1D filament model, evaluated using
adaptive quadrature.
"""
function force_filament_adaptive(coil::Coil, ϕ; reltol=1e-8, abstol=1e-14)
    regularization = compute_regularization(coil)
    r_eval, tangent = position_and_tangent(coil.curve, ϕ)
    B = B_filament_adaptive(coil::Coil, r_eval; regularization=regularization, reltol=reltol, abstol=abstol)
    return coil.current * cross(tangent, B)
end

"""
Compute the force-per-unit-length for a finite-thickness circular-cross-section coil.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness(coil::CoilCircularXSection, ϕ; reltol=1e-3, abstol=1e-5)
    dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)

    function force_cubature_func(xp)
        ρ = xp[1]
        θ = xp[2]
        r = ρ * coil.aminor
        sinθ, cosθ = sincos(θ)
        sqrtg = (1 - κ * r * cosθ) * ρ
        r_eval = r0 + r * cosθ * normal + r * sinθ * binormal
        B = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            θ_shift=θ,
        )
        return sqrtg * cross(tangent, B)
    end

    force_xmin = [0, 0]
    force_xmax = [1, 2π]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol,
        rtol=reltol
    )
    Biot_savart_prefactor = coil.current * μ0 / (4 * π^2)
    force_prefactor = coil.current / π
    return Biot_savart_prefactor * force_prefactor * val
end

"""
Compute the force-per-unit-length for a finite-thickness rectangular-cross-section coil.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness(coil::CoilRectangularXSection, ϕ; reltol=1e-3, abstol=1e-5)
    dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    p, q = get_frame(coil.frame, ϕ, r0, tangent, normal)
    κ1, κ2 = CoilForces.get_κ1_κ2(p, q, normal, κ)

    function force_cubature_func(xp)
        u = xp[1]
        v = xp[2]
        u_a_over_2 = 0.5 * u * coil.a
        v_b_over_2 = 0.5 * v * coil.b
        sqrtg = 1 - u_a_over_2 * κ1 - v_b_over_2 * κ2
        r_eval = r0 + u_a_over_2 * p + v_b_over_2 * q
        B = B_finite_thickness_normalized(
            coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
        )
        return sqrtg * cross(tangent, B)
    end

    force_xmin = [-1, -1]
    force_xmax = [1, 1]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol,
        rtol=reltol
    )
    Biot_savart_prefactor = μ0 * coil.current / (16 * π)
    force_prefactor = 0.25 * coil.current
    return Biot_savart_prefactor * force_prefactor * val
end

"""
Compute the force-per-unit-length for a finite-thickness coil. This version of
the function uses a single call to HCubature to handle all 5 dimensions.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness_5D(coil::CoilCircularXSection, ϕ; reltol=1e-3, abstol=1e-5)
    dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)

    function force_cubature_func(xp)
        ρ = xp[1]
        θ = xp[2]
        ρp = xp[3]
        θp = xp[4]
        ϕp = xp[5]
        dℓdϕp, κp, τp, rp_minus_r, tangentp, normalp, binormalp = Frenet_frame(coil.curve, ϕp)
        s = ρ * coil.aminor
        sp = ρp * coil.aminor
        sinθ, cosθ = sincos(θ)
        sinθp, cosθp = sincos(θp)
        @. rp_minus_r += (
            (sp * cosθp) * normalp
            + (sp * sinθp) * binormalp
            - (s * cosθ) * normal
            - (s * sinθ) * binormal
            - r0
        )
        temp = 1 / (normsq(rp_minus_r) + 1e-30)
        return (
            (ρ * ρp * (1 - κ * s * cosθ) * (1 - κp * sp * cosθp)
            * dℓdϕp * temp * sqrt(temp))
            * cross(tangent, cross(rp_minus_r, tangentp))
        )
    end

    #force_xmin = [0, 0, 0, 0, 0]
    #force_xmax = [1, 2π, 1, 2π, 2π]
    force_xmin = [0, 0, 0, 0, ϕ]
    force_xmax = [1, 2π, 1, 2π, ϕ + 2π]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol,
        rtol=reltol
    )
    Biot_savart_prefactor = coil.current * μ0 / (4 * π^2)
    force_prefactor = coil.current / π
    return Biot_savart_prefactor * force_prefactor * val
end

"""
Compute the force-per-unit-length for a finite-thickness coil. This version uses
Siena's trick of subtracting and adding the integrand for the best-fit circular coil.

ϕ: Curve parameter at which the force-per-unit-length will be computed.
"""
function force_finite_thickness_singularity_subtraction(coil::CoilCircularXSection, ϕ; reltol=1e-3, abstol=1e-5)
    dℓdϕ, κ, τ, r0, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)

    best_fit_circle = fit_circle(coil.curve, ϕ)
    best_fit_circular_coil = CoilCircularXSection(best_fit_circle, coil.current, coil.aminor)

    function force_cubature_func(xp)
        ρ = xp[1]
        θ = xp[2]
        r = ρ * coil.aminor
        sinθ, cosθ = sincos(θ)
        sqrtg = (1 - κ * r * cosθ) * ρ
        r_eval = r0 + r * cosθ * normal + r * sinθ * binormal
        B = B_finite_thickness_singularity_subtraction(
            coil,
            best_fit_circular_coil,
            r_eval,
            reltol=reltol,
            abstol=abstol,
            ϕ_shift=ϕ,
            θ_shift=θ,
        )
        return sqrtg * cross(tangent, B)
    end

    force_xmin = [0, 0]
    force_xmax = [1, 2π]
    
    val, err = hcubature(
        force_cubature_func, 
        force_xmin,
        force_xmax;
        atol=abstol,
        rtol=reltol
    )
    Biot_savart_prefactor = coil.current * μ0 / (4 * π^2)
    force_prefactor = coil.current / π
    integral = Biot_savart_prefactor * force_prefactor * val
    force_from_best_fit_circle = (
        -normal * 
        interpolated_force_per_unit_length(CoilCircularXSection(CurveCircle(1 / κ), coil.current, coil.aminor))
    )
    total = integral + force_from_best_fit_circle
    @show integral
    @show force_from_best_fit_circle
    @show total
    return integral, force_from_best_fit_circle, total
end

"""
Compute the self-force per unit length, using the locally circular
approximation, Garren & Chen eq (34).
"""
function force_locally_circular_approximation(coil::CoilCircularXSection, ϕ)
    differential_arclength, curvature, torsion, position, tangent, normal, binormal = Frenet_frame(coil.curve, ϕ)
    circular_coil = CoilCircularXSection(CurveCircle(1 / curvature), coil.current, coil.aminor)
    return -normal * analytic_force_per_unit_length(circular_coil)
end