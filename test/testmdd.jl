## --- Import an MDD dataset

    datapath = joinpath("..", "examples", "ol13-mdd.csv")
    mdds = importdataset(datapath, importas=:Tuple)

    agesteps = 995:-10.:5
    tsteps = reverse(agesteps)
    Tsteps = [fill(320., 25); fill(130., 75)]    

## --- Test Multiple Domain Diffusion with with PlanarAr
    mdd = MultipleDomain(Float64, PlanarAr;
        age = mdds.age_Ma,
        age_sigma = mdds.age_sigma_Ma,
        fraction_released = mdds.fraction_degassed,
        tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
        Tsteps_experimental = mdds.temperature_C,
        fit = mdds.fit,
        volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
        agesteps,
    )
    @test mdd isa MultipleDomain{Float64, PlanarAr{Float64}}
    show(mdd)
    println()
    display(mdd)

    tdomains = .!isnan.(mdds.lnD0_a_2)
    dm = MDDiffusivity(
        D0 = (Float64.(exp.(mdds.lnD0_a_2[tdomains]).*(100/10000)^2)...,),
        D0_logsigma = (Float64.(haskey(mdds, :lnD0_a_2_sigma) ? mdds.lnD0_a_2_sigma[tdomains] : fill(log(2)/2, count(tdomains)))...,),
        Ea = Float64(nanmean(mdds.Ea_kJ_mol)),
        Ea_logsigma = Float64(haskey(mdds, :Ea_logsigma) ? nanmean(mdds.Ea_logsigma) : log(2)/2),
    )

    age, fraction = modelage(mdd, Tsteps, dm)
    # println(round.(age, sigdigits=5))
    # println(round.(fraction, sigdigits=5))
    @test round.(age, sigdigits=5) == [154.94, 217.94, 289.32, 358.84, 414.61, 460.5, 509.29, 547.24, 583.79, 616.32, 640.09, 662.0, 679.38, 692.71, 705.09, 714.04, 722.54, 729.61, 733.62, 736.62, 738.7, 740.32, 742.11, 743.55, 744.96, 746.07, 746.61, 747.06, 747.81, 749.02, 750.22, 750.56, 750.11, 750.41, 750.53, 749.47, 748.91, 749.65, 750.43, 750.88, 750.97, 751.02, 751.43, 752.13, 752.94, 753.67, 754.36, 754.96, 755.55, 756.07, 756.58, 757.05, 757.53, 757.96, 758.42, 758.84, 759.3, 759.71, 760.17, 760.59, 761.05, 761.47, 761.94, 762.35, 762.82, 763.23, 763.7, 764.11, 764.58, 764.98, 765.45, 765.84, 766.31, 766.69, 767.15, 767.53, 767.99, 768.36, 768.81, 769.17, 769.62, 769.97, 770.43, 770.85, 771.46, 772.1, 772.94, 773.87, 775.02, 776.31, 777.87, 779.73, 781.96, 784.52, 787.62, 791.45, 796.25, 802.46, 810.99, 892.04]
    @test round.(fraction, sigdigits=5) == [0.00017134, 0.00049258, 0.0009698, 0.0015041, 0.0019328, 0.0025828, 0.0032327, 0.0039057, 0.0049035, 0.0057523, 0.0068191, 0.0081277, 0.0092747, 0.010963, 0.012691, 0.014649, 0.01806, 0.020947, 0.024381, 0.028374, 0.031873, 0.0376, 0.043318, 0.049084, 0.057303, 0.063745, 0.071628, 0.082742, 0.097153, 0.11363, 0.1299, 0.14469, 0.15855, 0.17196, 0.18456, 0.1972, 0.21123, 0.22671, 0.24298, 0.25906, 0.27435, 0.28964, 0.30532, 0.32034, 0.33309, 0.34413, 0.35384, 0.36249, 0.3703, 0.3774, 0.38392, 0.38995, 0.39555, 0.4008, 0.40572, 0.41037, 0.41477, 0.41895, 0.42292, 0.42672, 0.43034, 0.43381, 0.43714, 0.44034, 0.44342, 0.44638, 0.44924, 0.45199, 0.45465, 0.45723, 0.45972, 0.46213, 0.46448, 0.46675, 0.46895, 0.47109, 0.47318, 0.4752, 0.47718, 0.4791, 0.48098, 0.48281, 0.48471, 0.48699, 0.48967, 0.49281, 0.49644, 0.5006, 0.5053, 0.51056, 0.51656, 0.52367, 0.53131, 0.53987, 0.54981, 0.56162, 0.57526, 0.59266, 0.61401, 1.0]

    @test Thermochron.model_ll(mdd) ≈ -7337.040785211569

## --- Test Multiple Domain Diffusion with SphericalAr
    mdd = MultipleDomain(Float64, SphericalAr;
        age = mdds.age_Ma,
        age_sigma = mdds.age_sigma_Ma,
        fraction_released = mdds.fraction_degassed,
        tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
        Tsteps_experimental = mdds.temperature_C,
        fit = mdds.fit,
        volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
        agesteps,
    )
    @test mdd isa MultipleDomain{Float64, SphericalAr{Float64}}

    age, fraction = modelage(mdd, Tsteps, dm)
    # println(round.(age, sigdigits=5))
    # println(round.(fraction, sigdigits=5))
    @test round.(age, sigdigits=5) == [146.5, 206.93, 275.91, 343.66, 398.52, 443.89, 492.63, 530.82, 567.9, 601.29, 625.87, 648.75, 667.16, 681.36, 694.79, 704.61, 714.1, 722.43, 727.57, 732.12, 735.92, 738.96, 741.87, 743.52, 744.95, 746.25, 747.28, 748.31, 749.04, 749.08, 749.32, 749.53, 748.93, 748.65, 749.07, 749.46, 749.61, 749.52, 749.69, 750.23, 750.89, 751.46, 751.95, 752.41, 752.95, 753.51, 754.11, 754.63, 755.17, 755.63, 756.11, 756.51, 756.96, 757.33, 757.76, 758.1, 758.53, 758.84, 759.26, 759.56, 759.97, 760.25, 760.66, 760.92, 761.32, 761.56, 761.96, 762.18, 762.58, 762.78, 763.18, 763.36, 763.75, 763.93, 764.32, 764.48, 764.87, 765.01, 765.41, 765.54, 765.93, 766.06, 766.46, 766.65, 767.13, 767.49, 768.1, 768.67, 769.47, 770.3, 771.38, 772.6, 774.13, 775.82, 777.91, 780.4, 783.46, 787.15, 791.71, 809.16]
    @test round.(fraction, sigdigits=5) == [0.00051295, 0.0014691, 0.0028787, 0.0044427, 0.0056867, 0.0075558, 0.0094035, 0.011295, 0.014059, 0.016372, 0.019231, 0.022664, 0.025607, 0.029827, 0.03401, 0.038581, 0.046138, 0.052108, 0.058747, 0.065895, 0.071743, 0.080805, 0.089484, 0.098115, 0.11023, 0.11934, 0.12955, 0.14194, 0.15536, 0.16944, 0.18429, 0.1994, 0.21523, 0.23292, 0.25239, 0.27304, 0.2946, 0.31734, 0.34162, 0.36628, 0.38935, 0.41106, 0.43168, 0.45023, 0.46537, 0.47819, 0.48927, 0.49902, 0.50772, 0.51556, 0.52271, 0.52927, 0.53533, 0.54098, 0.54626, 0.55123, 0.55591, 0.56035, 0.56456, 0.56858, 0.57241, 0.57608, 0.5796, 0.58299, 0.58625, 0.58939, 0.59242, 0.59536, 0.5982, 0.60095, 0.60362, 0.60622, 0.60875, 0.61121, 0.6136, 0.61594, 0.61821, 0.62044, 0.62261, 0.62474, 0.62681, 0.62885, 0.63097, 0.63352, 0.63654, 0.64009, 0.64422, 0.64897, 0.65439, 0.6605, 0.6675, 0.67587, 0.6849, 0.69507, 0.70688, 0.72093, 0.73712, 0.7577, 0.78295, 1.0]

    @test Thermochron.model_ll(mdd) ≈ -8171.693349187755

## --- End of File