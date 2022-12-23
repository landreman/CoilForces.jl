# Data taken from simsopt/src/simsopt/configs/

# The file is expected to have :mod:`6*num_coils` many columns, and :mod:`order+1` many rows.                                    
# The columns are in the following order,                                                                                        
#                                                                                                                                      
#    sin_x_coil1, cos_x_coil1, sin_y_coil1, cos_y_coil1, sin_z_coil1, cos_z_coil1, sin_x_coil2, cos_x_coil2, sin_y_coil2, cos_y_coil2, sin_z_coil2, cos_z_coil2,  ...                                                                                                 

hsx_data = [0.0 1.450098613454403 0.0 0.08136481407962212 0.0 0.05214770022078918 0.0 1.364247137086304 0.0 0.2291881892946259 0.0 0.141555575302018 0.0 1.227084940141409 0.0 0.3586931962017947 0.0 0.1812133132220058 0.0 1.082927842630906 0.0 0.4758899202540535 0.0 0.1555085659916593 0.0 0.9355779318109653 0.0 0.5673221855979947 0.0 0.1058974345131323 0.0 0.7823601346256718 0.0 0.6604489927781515 0.0 0.03900805329348292;
    -0.2474148718390936 -0.04461965537300054 -0.01674239120693434 -0.178507840201461 -0.07886606457559459 0.3085349115227458 -0.1974233329122063 -0.06669735349369871 -0.08913598450558881 -0.1786154034259098 -0.1292589342455278 0.3199953239803489 -0.1256511919009446 -0.1161431876512416 -0.1301606105486485 -0.191713627052413 -0.1595654123593888 0.3145254445241458 -0.102050611782696 -0.1535694402675996 -0.1904267555957938 -0.1733303226854608 -0.1061484295706132 0.2988994978977337 -0.03573876176683196 -0.2083415271495969 -0.1721760219506088 -0.2016856729870607 -0.1629565237806347 0.1893335208216264 0.0887865022410676 -0.2096829037069601 -0.1060663936328356 -0.2299734515957785 -0.2027704656850164 0.03500584354797522;
    0.008639155493489753 -0.01692342738512217 -0.01290808090814132 -0.002791575075490101 -0.002978206172492419 -0.0006344328911988881 0.01527912552065864 -0.02175008019785673 -0.01729199381752107 0.004372210684036217 -0.005517960825466581 -0.007823005001863326 0.02756159111810172 -0.02378121602804106 -0.0245611316827341 0.00842360383616373 -0.01117032245970599 -0.01381169269372599 0.02924447647306906 -0.04204330412442347 -0.03247475930326364 0.02536260299442666 -0.02120123811426113 -0.01443918590581418 0.03428636235121949 -0.02114059341851812 -0.03450342033342272 0.02572691824805218 -0.006530133227560186 -0.03417683938813579 0.02254915988293732 0.02567851154978407 -0.03607472912837986 0.02583036706219939 0.01429624531492022 -0.01644701317182896;
    -0.01782651830957457 -0.01274394989177814 -0.00510964909294264 0.0117200057205984 -0.01140971106185065 0.01730796933182038 -0.01481605023934182 -0.0252097381760331 -0.001783820262763865 0.01551555293370957 -0.017002973602373 0.01556382835128874 -0.01965487471167211 -0.03533803100230879 0.007427001967799341 0.01523868961941819 -0.01878714811941281 0.00727964954297312 -0.04564091633344888 -0.0115786670537373 0.0187789529743671 -0.0094878368221247 -0.01798909416200898 0.006134558250282037 -0.05666414887466208 -0.009450306710632612 0.02699992766353758 -0.01161512202717966 -0.03620884906608093 -0.002493389824819134 -0.04788589502019225 -0.02406231942387799 0.03765166957373955 0.00890159824433479 -0.0557677668575913 -0.02267381187523191;
    -0.001713612259021574 -0.001329234214928187 0.001086180429829237 0.003851665101090454 0.004649779007888845 0.005580223812342257 -0.005575157149165955 -0.005533788813564074 0.003523331921677385 0.004948819537642362 0.004142935338351451 -0.0009649967357910503 -0.003656819696196109 -0.00658550834859699 -0.004122277268822814 0.001333080251928467 -0.005976532519555423 -0.009680801185971106 -0.003566297707787713 -0.01302678745117686 0.003111427768135027 0.01715019528369476 -0.0170929527924024 5.313324006256378e-05 -0.008350054464149013 -0.0118974257508249 0.0101040607268655 0.01926365824421457 -0.01624724700814385 -0.005181103925880991 -0.01044335236185535 -0.007827068847797295 0.0082044147764334 0.0158436391348558 -0.009645499906782559 -0.01148708053573716;
    -0.0001005823284171784 -0.003361073109820117 -0.002697214237357881 0.007708455798920935 -0.0007688213046379142 0.001193485150782259 -0.002590784189928517 -0.01258402429136424 0.005922827615547483 0.01385793870517058 -0.001353280337312436 0.002862858878026338 -0.01225110774344057 -0.01332615675468807 0.01683055799172453 0.009760359932945375 -0.0008202411158904526 -0.0007372613777052953 -0.005897417361324276 0.008358250409763394 0.004209492439939518 -0.01235861909151007 -0.00439109348622957 0.002610746577210686 -0.003678458872270454 -0.00144986983153489 0.005807018093089629 -0.004242657493404638 -0.005638896048838326 0.001692761588856276 -4.947012069535397e-05 -0.009503185073399924 0.01042617621619768 0.0002783454760244866 -0.002090919875529879 0.0002137808068046487;
    -0.002779551385348589 0.0007964801155168525 0.00205279160071562 -7.698854587557718e-06 0.001189881528429651 0.00106888912362874 -0.002500698264242391 0.002107028484299841 0.000592354443011868 -0.001745870888836288 -0.0001215045070470313 -0.001288060341062149 -0.000309894354778557 -0.001241549750255795 -0.002801887335754751 0.0008843996817371436 -0.002824972281414373 -0.002964174820369573 0.001032310461887095 0.003296426024793039 0.0006720825392066729 0.0001307766182670596 -0.004050328910271539 0.00520008141733631 0.007496416125782325 -0.001458264557268063 0.0004423800679911082 0.001657636436674581 0.002609528561208005 0.005241213815861743 0.006056394980210278 -0.001940488356487237 0.002216926181476822 -0.002801322986631159 -0.0004643300651748926 0.005750182143472103;
    0.0006263752364409071 0.000129940774886324 -0.00163125192988831 0.002785032061091969 -0.000284674973216022 0.0008694056392647126 -0.002158943507583175 -0.0007872101938002648 0.001884688529598002 0.002247429728758916 -0.0009892634543917495 0.001254584850019921 -0.001071803012984214 0.001291637058264399 0.001814189615140381 0.001440319899452729 -0.001536831007524525 0.003162484128655432 -0.001286638602494532 0.001362422327938524 0.001239202385108663 -0.004111166929599547 0.001137926200100212 0.001445640997087717 0.003383187538259416 0.003024398476110939 -0.0009600866175324597 -0.00270922285203278 0.0008881105046378972 0.000998329154419136 0.001308322590765889 0.0007682214474447779 0.00179846602187101 -0.0008734672783231939 0.002221790826586932 0.003095639792316356;
    -0.0007091831341319518 -0.0002887381259876946 0.002122890161836885 0.001173377991986841 0.0002964944145617011 0.000715670204997118 -0.001190378710850545 0.0003404123275896269 0.001470160064840406 -0.0007890256076894956 -0.0002148109469443335 0.0008623555735137005 -0.00108823524135118 0.001655893867141413 0.0003921822728934221 -0.002371968247859196 0.0005357276703431997 -5.212924307717625e-05 0.001612086576571139 -0.001968910846479319 0.001686281995967648 0.002089679345878055 0.003518606759861904 -5.078851137416556e-05 -0.002509633568768933 0.0006596228018896409 0.001224699131977734 0.001319576381541508 0.002072882840815181 0.001622869914653905 -0.001252326872163296 -0.0003205740497308463 -0.0007625689316982401 0.0001410294334229526 0.0005695829043695982 0.0003739678058129052;
    -0.0003985332323633957 7.780520749469628e-06 0.0002338203305055565 -0.0002429424084434014 -5.05642925570683e-05 0.0004561710800827271 -0.0008149086684595709 0.001262910578200743 -0.0009018927872027746 -0.0008235242066047514 -0.001166724981167916 0.0007221123418686274 0.000538492179785068 -0.0003877680252081154 -0.001206340728889653 0.001909579817956973 -0.001059350037691672 0.003030327583891618 0.0009309734675347164 -0.001402825213888862 -0.001879042720235417 -0.0005745632366927774 -0.0003229201151405956 -0.001214373269623841 -0.001339080795433758 0.0008869557187344778 0.0003471888782743237 -0.0005652715481477609 0.0006653815524971561 -0.002070844989207776 -0.001642162320901585 -0.000111618397365698 -0.001477859599048111 -0.0006400283139432591 0.0002189419776456112 -0.000816675485623738;
    -7.816153465351116e-05 -0.0002009847452342995 0.0008066240866545759 0.0008255436371530414 -7.807687546976041e-05 0.0004164853745811902 -0.0005523393805816584 -2.481935603354514e-05 0.001068182575791408 -0.0001796510186862261 0.0005830177580454353 0.001037449652217471 -3.481181167448667e-05 0.001224700761266774 -0.0005964133742786682 -0.001201968652425583 0.0004584825434668106 0.0008277320380496165 -0.0001635139911672491 0.0009570524773483278 0.0007313806204762912 7.451489265429876e-05 0.001088344781917036 2.195367034077586e-05 -0.0007958553060898448 -0.0008206404810253051 4.110418379880252e-05 -0.0008067982905227205 0.000290915451684831 -5.816471482555755e-05 -0.001179934184194162 -0.0003118344793459686 -0.0006285268670711766 -0.0003463904065092548 -0.0005427880030082414 -0.000274092136038884;
    -0.0005143536944340813 -0.0002653804742816889 0.0008432390106404718 -0.0006395291636203084 -1.885868142699404e-05 0.0001600178306523682 -0.0001079444298802007 0.0001277182335311201 -0.0006678046868030895 -0.0002031801878276445 -0.000465764066917379 0.0004752364701527337 -0.000897943901887332 0.0001030895987594508 0.0007741191253401746 0.0001050518553732099 0.0003018119095117473 0.000469768267408582 -0.0003849172956155816 0.0003235392804654553 -4.460931324305901e-05 -0.0008787075038936054 -0.0009325144595989179 -0.001254629538565482 2.880036928857217e-05 -0.0009426420778578429 -0.0004160723461747537 -0.000165188182301544 -0.001057146557965285 0.0008056165301685989 -0.0005533166696607232 0.0004519268108321537 -0.0001860173689221575 -0.000343212260066954 -0.0002493025774549491 0.0003056097763019185;
    -0.0002249696866484157 7.805760915540403e-05 8.18676545547803e-06 3.346093675889216e-06 -0.000172992119171605 0.0002046044960843123 -0.000138147484032206 0.000353788909390999 0.0001192390201640588 5.764802453573258e-05 0.0002485577511343327 0.0003634965826352155 -0.0001807397492647728 -0.0007787058616840801 9.461973304643845e-05 3.055966420335239e-05 0.0001814463246416648 0.0003136669593223852 6.582763262167127e-05 -0.0006788981862036736 1.811810993827302e-05 0.0003620652311326771 -1.430939719224026e-05 0.0001732670665024719 8.238891006000932e-05 0.0003855537510820345 0.0001839120592283135 -0.0004312197265379137 8.182854639908471e-05 0.0006411090093275426 0.0003017080047081547 1.588132153691171e-05 0.0003972735397212805 0.0002006485115939592 0.0002264387681962563 7.06503581087673e-05;
    -0.0002063318974367847 -0.0001249385251030611 0.0004453938734729784 -0.0001937250997854732 -4.808984834003181e-07 9.771903826810657e-05 -4.316157847671333e-05 -0.0002325012629261361 5.157623328009618e-05 0.0001553063328202026 0.0001746716554652558 9.859165114730268e-05 -0.0001453112882061044 0.0005010456663944343 0.0001488895313200456 -5.721432125520176e-05 -0.0002087524037308005 1.497745248440838e-05 -4.817812479347767e-06 -0.0002034223835519717 -0.0001802720281604667 0.0005830545568873165 -0.0002255529900934968 0.0004470674243802263 0.0004198928363364515 1.28087250487727e-05 9.348540646234649e-05 0.0001464152415004175 -0.0004142652946724504 7.247086362075446e-05 0.0002718741779593693 -4.934741920553801e-05 0.0002299427964981626 0.0001331083801620539 -8.332189386295343e-06 0.00018123585109773;
    -0.0003022331540395053 6.680722019955035e-05 -7.759400441490608e-05 -0.0002606117777834632 -5.865285824465475e-05 1.902339492824271e-05 -2.849696512246128e-05 0.000189391557998611 -9.977634143645547e-06 0.000170430001791046 -2.857942668359551e-05 3.619917479099263e-05 -0.0001617847678940659 -0.0001410318145757181 0.0002568343323106022 -0.0003087333061012181 0.0002110001009470078 -0.0003790910290479775 0.0003217745078527494 9.189953174293771e-05 -4.106886137757557e-05 -0.0001158802046133196 2.075653703520431e-05 0.0001359738959075968 -0.0002009437823238153 9.407771053133783e-05 -0.0001148928752906213 0.0001230815527775522 7.25051251997194e-05 -5.612240698397028e-05 -0.0002205645583349634 0.0002298889303352109 0.0003224120921006623 -0.0001006733374779524 -0.0004502249100044474 0.0005136810995596343;
    -7.666183160489453e-05 0.0001136223997908341 3.069085515968562e-07 2.96890093442947e-05 -1.028679305543576e-05 8.668517139169007e-06 -5.437154749210739e-05 1.283500615253788e-05 0.0001018590133126054 1.656525546473645e-05 0.000133687295960366 -0.0001858982193899168 -8.154044695306804e-05 -0.0001698563244539748 2.624082998513663e-05 0.0002977102274138075 -0.0001384903749562772 4.92073595853411e-05 -0.000243900370436856 0.0002706638259667449 7.361366235023523e-05 -8.51103040405823e-05 -0.0002413901343175848 9.874931637333237e-05 -0.0002003356389102572 0.0003620062305069666 0.0001633096104663104 8.432143777860114e-05 0.0002213289836753511 0.000155975466256031 -0.000234467228430851 0.0003524549545886347 0.0003309458625387029 -0.0001547923718336705 -0.0003302058453151174 0.0003368083194438543;
    -0.0001298238455857419 3.907158056675016e-05 -2.348146360607329e-05 -8.239288767455298e-05 3.720161589039513e-05 -6.153182697136335e-05 -4.754938710183249e-05 2.047011952541541e-05 9.406034782419945e-05 0.0001028785298997861 1.078918544660184e-05 -5.934307167302391e-05 0.0001843476758627078 0.0001721921213608038 -0.0001726992877228089 6.606758335251639e-06 8.922368136072931e-06 -7.938773575792159e-05 -0.0001072589422894671 2.939566518210593e-05 0.0001541989694809523 -5.394429170317401e-05 0.0002824901565846145 -0.0001060091619014223 9.61640428293076e-05 0.0002953095261748043 5.818678486792215e-05 -0.0001511864613527624 0.0001097228573362343 -0.0001554801536926139 -4.135660214455158e-05 -4.427481130821969e-05 0.0001959685127252026 -0.0001195053746517015 0.0001027846684114793 -2.262305715484092e-05]

function get_curve(config, coil_num)
    if config == "hsx"
        data = hsx_data
    else
        throw(ArgumentError("Unrecognized configuration."))
    end

    if coil_num < 1 || 6 * coil_num > size(data, 2)
        throw(ArgumentError("Invalid coil_num"))
    end
    data = data[:, 1 + (coil_num - 1) * 6: coil_num * 6]
    xs = data[:, 1]
    xc = data[:, 2]
    ys = data[:, 3]
    yc = data[:, 4]
    zs = data[:, 5]
    zc = data[:, 6]
    return CurveXYZFourier(xc, xs, yc, ys, zc, zs)
end
