struct Coil
    curve::Curve
    current::Float64
    aminor::Float64
end

function analytic_force_per_unit_length(coil::Coil)
    # Assert that curve type is a CurveCircle:
    coil.curve::CurveCircle
    I = coil.current
    R = coil.curve.R0
    a = coil.aminor
    return μ0 * I * I / (4π * R) * (log(8 * R / a) - 0.75)
end