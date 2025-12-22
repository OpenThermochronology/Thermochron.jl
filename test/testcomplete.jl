## --- Prepare problem

    # Read in data from file using StatGeochem
    datapath = joinpath("..", "examples",)
    ds = importdataset(joinpath(datapath, "exampledata", "minnesota.csv"), ',', importas=:Tuple);

    using LinearAlgebra
    BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

    params = (
        burnin = 350,               # [n] How long should we wait for MC to converge (become stationary)
        nsteps = 350,               # [n] How many steps of the Markov chain should we run after burn-in?
        dr = 1.0,                   # [μ] Radius step size
        dt = 10.0,                  # [Ma] time step size
        dTmax = 25.0,               # [C/step] Maximum reheating/burial per params timestep
        Tinit = 400.0,              # [C] initial model temperature (i.e., crystallization temperature)
        ΔTinit = -50.0,             # [C] Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0,                 # [C] Current surface temperature
        ΔTnow = 10.0,               # [C] Tnow may vary from Tnow to Tnow+ΔTnow
        tnow = 0.0 ,                # [Ma] Today
        minpoints = 1,              # [n] Minimum allowed number of t-T points
        maxpoints = 40,             # [n] Maximum allowed number of t-T points
        npoints = 5,                # [n] Initial number of t-T points
        simplified = false,         # Prefer simpler tT paths?
        dynamicsigma = true,        # Update model uncertainties?
        # Diffusivity and annealing models
        zdm = ZRDAAM(), 
        adm = RDAAM(), 
        zftm = Yamada2007PC(), 
        mftm = Jones2021FA(), 
        aftm = Ketcham2007FC(),
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
        σmodel = 25.0,              # [Ma] assumed model uncertainty (resampled if dynamicsigma)
        T0annealing = 1,            # [unitless] initial annealing "temperature"
    )

    # Populate data NamedTuple from imported dataset
    data = deepcopy(ds)
    data.raw_He_age_sigma_Ma .= 0.1*data.raw_He_age_Ma

    # Crystallization ages and start time
    tinit = ceil(maximum(data.crystallization_age_Ma)/params.dt) * params.dt
    params = (params...,
        tinit = tinit,
        agesteps = (tinit-params.dt/2 : -params.dt : 0+params.dt/2),
        tsteps = (0+params.dt/2 : params.dt : tinit-params.dt/2),
    )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = [params.tnow, params.tinit],   # [Ma] Final and initial time
        T₀ = [params.Tnow, params.Tinit],          # [C] Final and initial temperature
        ΔT = [params.ΔTnow, params.ΔTinit],        # [C] Final and initial temperature range (positive or negative)
        tboundary = :reflecting, # Reflecting time boundary conditions
        Tboundary = :reflecting, # Reflecting temperature boundary conditions
    )

    # Default: No unconformity is imposed
    unconf = Constraint()

## --- Test generation and modelling of Chronometer objects, with StepRangeLen for timesteps/agesteps

    # Modern input format, generic.csv
    tsteps = (params.dt/2 : params.dt : 3000-params.dt/2)
    agesteps = (3000-params.dt/2 : -params.dt : params.dt/2)
    Tsteps = range(650, 0, length=length(tsteps))
    params = (params..., agesteps = agesteps, tsteps = tsteps)

    dsg = importdataset(joinpath(datapath, "generic.csv"), ',', importas=:Tuple)
    chrons, damodels = chronometers(dsg, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == 27
    const FloatRange = typeof(1.0:1.0:10.0)
    @test count(x->isa(x,SphericalHe{Float64, FloatRange}), chrons) == 4
    @test count(x->isa(x,PlanarHe{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,SphericalAr{Float64, FloatRange}), chrons) == 2
    @test count(x->isa(x,PlanarAr{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,ZirconHe{Float64, FloatRange}), chrons) == 3
    @test count(x->isa(x,ApatiteHe{Float64, FloatRange}), chrons) == 4
    @test count(x->isa(x,ZirconFT{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,MonaziteFT{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,ApatiteFT{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,ZirconTrackLength{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,MonaziteTrackLength{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLength{Float64, FloatRange}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLengthOriented{Float64, FloatRange}), chrons) == 3
    @test count(x->isa(x,SingleDomain), chrons) == 1
    @test count(x->isa(x,MultipleDomain), chrons) == 2
    @test get_age(chrons) ≈ [150.37, 263.92, 150.37, 263.92, 263.92, 917.84, 1023.73, 1023.73, 380., 380., 120., 120., 1080., 300., 100., 150., 180., 4.1, 0.9194109843673132, 808.3268143245239, 808.3268143245239,] 
    @test get_age_sigma(chrons) ≈ [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0.53,0.20877153500779683,28.52408719185519,28.52408719185519,]

    # Test model calculations
    calc, calcuncert, ll = model(chrons, damodels, Tsteps)
    @test ll ≈ -121502.38423678318
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -127582.33244603685
    @test round.(calc[1:18], sigdigits=7) ≈ [100.512, 196.5576, 110.1727, 199.4224, 195.2399, 868.0376, 969.4693, 962.8585, 286.9455, 289.894, 242.1764, 276.2001, 1085.555, 304.6573, 95.84216, 149.8249, 297.8784, 262.766]
    @test calc[19] ≈ 0.8 atol=0.5
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.7578982633970572, 1.1785910438098226, 1.1389520917140208, 1.2018361658877996, 1.1302107318562702, 0.6070538659171328]

    # Test again after swapping out annealing models
    damodels = Thermochron.Model[damodels...,]
    damodels[isa.(damodels, Thermochron.ZirconAnnealingModel)] .= Guenthner2013FC()
    damodels[isa.(damodels, Thermochron.ApatiteAnnealingModel)] .= Ketcham1999FC()
    damodels = unionize(damodels)
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -121520.27640274604 
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -127600.2246119997
    @test round.(calc[1:18], sigdigits=7) ≈ [100.512, 196.5576, 110.1727, 199.4224, 195.2399, 868.0376, 969.4693, 962.8585, 286.9455, 289.894, 242.1764, 276.2001, 1110.379, 304.2772, 95.84216, 149.8249, 297.8784, 262.766]
    @test calc[19] ≈ 0.8 atol=0.5
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.8368172844202661, 1.1896389981502726, 1.1448424397109467, 1.2154485905638788, 1.1896389981502726, 0.6070538659171328]

    # Test kintetic_ll! and updatekinetics! on all chronometer types
    updatekinetics = falses(length(damodels))
    @test Thermochron.kinetic_ll!(updatekinetics, damodels, damodels) ≈ 19.64330130235549
    damodelsₚ = copy(damodels)
    Thermochron.movekinetics!(damodelsₚ, updatekinetics)
    Thermochron.movekinetics!(damodelsₚ, updatekinetics)
    @test 10 < Thermochron.kinetic_ll!(updatekinetics, damodelsₚ, damodels) < 19.64330130235549

    # Test again with partitiondaughter=true
    chrons, damodels = chronometers(dsg, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; partitiondaughter=true) ≈ -121512.04499535152
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; partitiondaughter=true, redegastracer=true) ≈ -127591.99320460518
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; partitiondaughter=true, redegastracer=true, stepwisetracerfraction=true) ≈ -121635.85117007693
    @test round.(calc[1:18], sigdigits=7) ≈ [100.512, 196.5576, 110.7153, 200.3158, 196.0983, 868.0376, 969.4693, 962.8585, 286.9455, 290.1031, 242.1764, 277.9584, 1085.555, 304.6573, 95.84216, 149.9689, 298.9059, 262.766]
    @test calc[19] ≈ 0.8 atol=0.5
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.7578982633970572, 1.1785910438098226, 1.1389520917140208, 1.2018361658877996, 1.1302107318562702, 0.6070538659171328]

    # Test empirical uncertainty estimation
    σcalc = zeros(length(chrons))
    empiricaluncertainty!(σcalc, chrons, ZirconHe, sigma_offset=15)
    empiricaluncertainty!(σcalc, chrons, ApatiteHe, sigma_offset=15)
    t = isa.(chrons, ZirconHe) .| isa.(chrons, ApatiteHe)
    @test σcalc[t] ≈ [57.014755105003104, 57.014755105003104, 14.022789552152355, 14.022789552152355, 89.77854252945355, 22.081049122650942, 18.76032805406918]

    # Test again with Minnesota dataset
    tsteps = (params.dt/2 : params.dt : params.tinit-params.dt/2)
    agesteps = (params.tinit-params.dt/2 : -params.dt : params.dt/2)
    Tsteps = range(params.Tinit, params.Tnow, length=length(tsteps))
    params = (params..., agesteps = agesteps, tsteps = tsteps)

    chrons, damodels = chronometers(ds, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == length(ds.mineral)
    @test get_age(chrons) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11, 234, 98, 233, 339, 378, 258, 158, 269, 313, 309, 392]
    @test get_age_sigma(chrons) ≈ [15.0, 51.0, 13.0, 12.6, 18.4, 10.0, 12.0, 10.0, 43.0, 7.0, 7.0, 5.2, 4.6, 4.6, 6.0, 5.0, 4.0, 4.0, 2.0, 1.1, 1.1, 0.8, 0.2, 9.4, 2.0, 5.0, 6.0, 9.3, 5.0, 3.3, 8.1, 5.8, 6.0, 8.0]

    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -52946.43831930029
    # @info round.(calc, sigdigits=5)
    @test round.(calc, sigdigits=5) ≈ [1125.8, 1123.3, 954.89, 1046.3, 1010.6, 1146.5, 956.84, 871.49, 984.5, 706.49, 574.3, 139.5, 319.64, 62.847, 183.99, 329.73, 1.712, 3.6475, 0.00037435, 16.555, 0.0041183, 33.839, 1.2216, 405.81, 551.61, 550.55, 571.48, 496.9, 534.06, 519.35, 536.78, 608.72, 598.43, 650.42]
    @test calcuncert ≈ zeros(length(chrons))
    @test first(calc) ≈ modelage(first(chrons), Tsteps, ZRDAAM())
    @test last(calc) ≈ modelage(last(chrons), Tsteps, RDAAM())

    # Test empirical uncertainty estimation
    σcalc = zeros(length(chrons))
    empiricaluncertainty!(σcalc, chrons, ZirconHe)
    @test get_age(chrons, ZirconHe) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11]
    σtotal = sqrt.(σcalc[isa.(chrons, ZirconHe)].^2 + get_age_sigma(chrons, ZirconHe).^2)
    @test σtotal ≈ [71.76511200022767, 88.69939042178478, 76.20513740085536, 49.74341557807492, 52.05638860615669, 58.80850040722076, 57.32853093748527, 97.59191128458991, 68.09272937102548, 82.41401338162268, 104.5147581608337, 72.10558616290515, 79.47676757701001, 68.01965448205827, 72.8583946436481, 73.94974380762895, 63.02840819654363, 63.478989858090436, 45.255278015066644, 66.3078727295436, 52.90379008918621, 69.14186488208391, 62.394435852786344]

    # Test again with Manitoba dataset
    tsteps = (params.dt/2 : params.dt : 2790-params.dt/2)
    agesteps = (2790-params.dt/2 : -params.dt : params.dt/2)
    Tsteps = range(650, 0, length=length(tsteps))
    params = (params..., agesteps = agesteps, tsteps = tsteps)

    dsm = importdataset(joinpath(datapath, "exampledata", "manitoba.csv"), ',', importas=:Tuple)
    chrons, damodels = chronometers(dsm, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == 333
    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -8.166180623895764e8

## --- Test generation and modelling of Chronometer objects, with nonuniform (log-spaced) Vector for agesteps/tsteps

    agesteps = cntr(logrange(3000+10, 10, length=Int(3000/10)) .- 10)
    tsteps = (first(agesteps) - Thermochron.step_at(agesteps,1)/2) .- agesteps
    Tsteps = agesteps * 650/3000
    params = (params..., agesteps=agesteps, tsteps=tsteps)

    chrons, damodels = chronometers(dsg, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == 27
    @test count(x->isa(x,SphericalHe{Float64, Vector{Float64}}), chrons) == 4
    @test count(x->isa(x,PlanarHe{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,SphericalAr{Float64, Vector{Float64}}), chrons) == 2
    @test count(x->isa(x,PlanarAr{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ZirconHe{Float64, Vector{Float64}}), chrons) == 3
    @test count(x->isa(x,ApatiteHe{Float64, Vector{Float64}}), chrons) == 4
    @test count(x->isa(x,ZirconFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,MonaziteFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ZirconTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,MonaziteTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLengthOriented{Float64, Vector{Float64}}), chrons) == 3
    @test count(x->isa(x,SingleDomain), chrons) == 1
    @test count(x->isa(x,MultipleDomain), chrons) == 2
    @test get_age(chrons) ≈ [150.37, 263.92, 150.37, 263.92, 263.92, 917.84, 1023.73, 1023.73, 380., 380., 120., 120., 1080., 300., 100., 150., 180., 4.1, 0.9194109843673132, 808.3268143245239, 808.3268143245239,] 
    @test get_age_sigma(chrons) ≈ [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0.53,0.20877153500779683,28.52408719185519,28.52408719185519,]

    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -118298.29723826758
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -124355.64361412608 
    @test round.(calc[1:18], sigdigits=7) ≈ [97.38443, 192.5678, 106.6764, 195.3374, 191.1125, 865.807, 967.5609, 960.928, 287.2698, 290.2198, 240.5851, 274.3691, 1078.846, 302.8478, 95.39966, 149.9576, 296.5752, 259.3508]
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.7573177921647993, 1.1845573491521355, 1.1457607696410228, 1.2049589019595262, 1.1414549313278457, 0.6188094901594741]

    damodels = Thermochron.Model[damodels...,]
    damodels[isa.(damodels, Thermochron.ZirconAnnealingModel)] .= Guenthner2013FC()
    damodels[isa.(damodels, Thermochron.ApatiteAnnealingModel)] .= Ketcham1999FC()
    damodels = unionize(damodels)
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -118309.37754858492
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -124366.72392444342
    @test round.(calc[1:18], sigdigits=7) ≈ [97.38443, 192.5678, 106.6764, 195.3374, 191.1125, 865.807, 967.5609, 960.928, 287.2698, 290.2198, 240.5851, 274.3691, 1103.519, 302.1657, 95.39966, 149.9576, 296.5752, 259.3508]
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.8360236234430574, 1.1942271553159216, 1.149391800085125, 1.2195745732151688, 1.1942271553159216, 0.6188094901594741]

    # Test again with Minnesota dataset
    agesteps = cntr(logrange(tinit+10, 10, length=Int(tinit/10)) .- 10)
    tsteps = (first(agesteps) - Thermochron.step_at(agesteps,1)/2) .- agesteps
    Tsteps = agesteps * (params.Tinit-params.Tnow)/length(tsteps) .+ params.Tnow
    params = (params..., agesteps=agesteps, tsteps=tsteps)

    chrons, damodels = chronometers(ds, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == length(ds.mineral)
    @test get_age(chrons) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11, 234, 98, 233, 339, 378, 258, 158, 269, 313, 309, 392]
    @test get_age_sigma(chrons) ≈ [15.0, 51.0, 13.0, 12.6, 18.4, 10.0, 12.0, 10.0, 43.0, 7.0, 7.0, 5.2, 4.6, 4.6, 6.0, 5.0, 4.0, 4.0, 2.0, 1.1, 1.1, 0.8, 0.2, 9.4, 2.0, 5.0, 6.0, 9.3, 5.0, 3.3, 8.1, 5.8, 6.0, 8.0]

    tsteps = maximum(params.agesteps)+minimum(params.agesteps) .- params.agesteps
    Tsteps = params.agesteps * (params.Tinit - params.Tnow)/tinit .+ params.Tnow

    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -52195.37929684592
    # @info round.(calc, sigdigits=5)
    @test round.(calc, sigdigits=5) ≈ [1124.5, 1120.2, 951.61, 1046.8, 1010.9, 1146.0, 958.95, 872.19, 986.33, 708.15, 579.2, 144.31, 324.63, 64.507, 189.0, 336.73, 2.0527, 4.3125, 0.00035165, 18.467, 0.0052616, 36.526, 1.332, 403.37, 548.11, 547.24, 566.45, 493.98, 530.91, 515.74, 533.55, 604.88, 594.81, 645.67]
    @test calcuncert ≈ zeros(length(chrons))
    @test first(calc) ≈ modelage(first(chrons), Tsteps, ZRDAAM())
    @test last(calc) ≈ modelage(last(chrons), Tsteps, RDAAM())

## --- Invert for t-T path via MCMC

    # Run Markov Chain
    @time "\nCompiling MCMC" MCMC(data, params, boundary, unconf)
    @time "\nRunning MCMC" tT = MCMC(data, params, boundary, unconf; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 300
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -800 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT : $(mean(tT.jTdist))"

## --- Test generation and use of Chronometer objects, switching back to a uniform Vector for agesteps/tsteps

    tsteps = collect(params.dt/2 : params.dt : 3000-params.dt/2)
    agesteps = collect(3000-params.dt/2 : -params.dt : params.dt/2)
    Tsteps = range(650, 0, length=length(tsteps))
    params = (params..., agesteps = agesteps, tsteps = tsteps)

    chrons, damodels = chronometers(dsg, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == 27
    @test count(x->isa(x,SphericalHe{Float64, Vector{Float64}}), chrons) == 4
    @test count(x->isa(x,PlanarHe{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,SphericalAr{Float64, Vector{Float64}}), chrons) == 2
    @test count(x->isa(x,PlanarAr{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ZirconHe{Float64, Vector{Float64}}), chrons) == 3
    @test count(x->isa(x,ApatiteHe{Float64, Vector{Float64}}), chrons) == 4
    @test count(x->isa(x,ZirconFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,MonaziteFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteFT{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ZirconTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,MonaziteTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLength{Float64, Vector{Float64}}), chrons) == 1
    @test count(x->isa(x,ApatiteTrackLengthOriented{Float64, Vector{Float64}}), chrons) == 3
    @test count(x->isa(x,SingleDomain), chrons) == 1
    @test count(x->isa(x,MultipleDomain), chrons) == 2
    @test get_age(chrons) ≈ [150.37, 263.92, 150.37, 263.92, 263.92, 917.84, 1023.73, 1023.73, 380., 380., 120., 120., 1080., 300., 100., 150., 180., 4.1, 0.9194109843673132, 808.3268143245239, 808.3268143245239,] 
    @test get_age_sigma(chrons) ≈ [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0.53,0.20877153500779683,28.52408719185519,28.52408719185519,]

    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -121502.38423678318
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -127582.33244603685
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true, stepwisetracerfraction=true) ≈ -121626.19041150859
    @test round.(calc[1:18], sigdigits=7) ≈ [100.512, 196.5576, 110.1727, 199.4224, 195.2399, 868.0376, 969.4693, 962.8585, 286.9455, 289.894, 242.1764, 276.2001, 1085.555, 304.6573, 95.84216, 149.8249, 297.8784, 262.766]
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.7578982633970572, 1.1785910438098226, 1.1389520917140208, 1.2018361658877996, 1.1302107318562702, 0.6070538659171328]

    damodels = Thermochron.Model[damodels...,]
    damodels[isa.(damodels, Thermochron.ZirconAnnealingModel)] .= Guenthner2013FC()
    damodels[isa.(damodels, Thermochron.ApatiteAnnealingModel)] .= Ketcham1999FC()
    damodels = unionize(damodels)
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -121520.27640274604
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps; redegastracer=true) ≈ -127600.2246119997
    @test round.(calc[1:18], sigdigits=7) ≈ [100.512, 196.5576, 110.1727, 199.4224, 195.2399, 868.0376, 969.4693, 962.8585, 286.9455, 289.894, 242.1764, 276.2001, 1110.379, 304.2772, 95.84216, 149.8249, 297.8784, 262.766]
    @test calc[20] ≈ 735 atol=85
    @test calc[21] ≈ 755 atol=65
    @test round.(calc[22:end], sigdigits=3) ≈ [9, 14.3, 14.3, 14.3, 14.3, 7] atol = 10
    @test calcuncert[1:21] ≈ zeros(21)
    @test calcuncert[22:end] ≈ [1.8368172844202661, 1.1896389981502726, 1.1448424397109467, 1.2154485905638788, 1.1896389981502726, 0.6070538659171328]

    # Modern input format, Minnesota dataset
    tsteps = collect(params.dt/2 : params.dt : params.tinit-params.dt/2)
    agesteps = collect(params.tinit-params.dt/2 : -params.dt : params.dt/2)
    Tsteps = range(params.Tinit, params.Tnow, length=length(tsteps))
    params = (params..., agesteps = agesteps, tsteps = tsteps)

    chrons, damodels = chronometers(ds, params, zirconvolumeweighting=:spherical, apatitevolumeweighting=:spherical)
    @test chrons isa Vector{<:Chronometer}
    @test length(chrons) == length(ds.mineral)
    @test get_age(chrons) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11, 234, 98, 233, 339, 378, 258, 158, 269, 313, 309, 392]
    @test get_age_sigma(chrons) ≈ [15.0, 51.0, 13.0, 12.6, 18.4, 10.0, 12.0, 10.0, 43.0, 7.0, 7.0, 5.2, 4.6, 4.6, 6.0, 5.0, 4.0, 4.0, 2.0, 1.1, 1.1, 0.8, 0.2, 9.4, 2.0, 5.0, 6.0, 9.3, 5.0, 3.3, 8.1, 5.8, 6.0, 8.0]

    calc = zeros(length(chrons))
    calcuncert = zeros(length(chrons))
    @test Thermochron.model!(calc, calcuncert, chrons, damodels, Tsteps) ≈ -52946.43831930029
    # @info round.(calc, sigdigits=5)
    @test round.(calc, sigdigits=5) ≈ [1125.8, 1123.3, 954.89, 1046.3, 1010.6, 1146.5, 956.84, 871.49, 984.5, 706.49, 574.3, 139.5, 319.64, 62.847, 183.99, 329.73, 1.712, 3.6475, 0.00037435, 16.555, 0.0041183, 33.839, 1.2216, 405.81, 551.61, 550.55, 571.48, 496.9, 534.06, 519.35, 536.78, 608.72, 598.43, 650.42]
    @test calcuncert ≈ zeros(length(chrons))
    @test first(calc) ≈ modelage(first(chrons), Tsteps, ZRDAAM())
    @test last(calc) ≈ modelage(last(chrons), Tsteps, RDAAM())

## --- Invert for t-T path, as above, but with variable kinetic parameters

    # Run Markov Chain
    @time "\nCompiling MCMC_varkinetics" MCMC_varkinetics(data, params, boundary, unconf)
    @time "\nRunning MCMC_varkinetics" tT, kinetics = MCMC_varkinetics(data, params, boundary, unconf; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 150
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -450 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"


    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT: $(mean(tT.jTdist))"

    # Kinetics
    admdist = kinetics.dmdist[end, :]
    D0Lmean = mean(admdist .|> x-> x.D0L)
    @test 0 < D0Lmean
    @info "Mean apatite D0L: $D0Lmean"
    EaLmean = mean(admdist .|> x-> x.EaL)
    @test 0 < EaLmean 
    @info "Mean apatite EaL: $EaLmean"
    EaTrapmean = mean(admdist .|> x-> x.EaTrap)
    @test 0 < EaTrapmean 
    @info "Mean apatite EaTrap: $EaTrapmean"
    rmr0mean = mean(admdist .|> x-> x.rmr0)
    @test 0 < rmr0mean < 1
    @info "Mean apatite rmr0: $rmr0mean"

    zdmdist = kinetics.dmdist[1, :]
    DzD0mean = mean(zdmdist .|> x-> x.DzD0)
    @test 0 < DzD0mean
    @info "Mean zircon DzD0: $DzD0mean"
    DzEamean = mean(zdmdist .|> x-> x.DzEa)
    @test 0 < DzEamean
    @info "Mean zircon DzEa: $DzEamean"
    DN17D0mean = mean(zdmdist .|> x-> x.DN17D0)
    @test 0 < DN17D0mean
    @info "Mean zircon DN17D0: $DN17D0mean"
    DN17Eamean = mean(zdmdist .|> x-> x.DN17Ea)
    @test 0 < DN17Eamean
    @info "Mean zircon DN17Ea: $DN17Eamean"
    rminmean = mean(zdmdist .|> x-> x.rmin)
    @test 0 < rminmean < 1
    @info "Mean zircon rmin: $rminmean"

## --- MCMC with Detail interval
    detail = DetailInterval(
        agemin = 0, # Youngest end of detail interval
        agemax = 1000, # Oldest end of detail interval
        minpoints = 3, # Minimum number of points in detail interval
    )
    @time "\nMCMC with Detail interval" tT = MCMC(data, params, boundary, unconf, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 150
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -450 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT: $(mean(tT.jTdist))"

## --- MCMC_varkinetics with Detail interval

    @time "\nMCMC_varkinetics with Detail interval" tT, kinetics = MCMC_varkinetics(data, params, boundary, unconf, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 200
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -450 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"


    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT: $(mean(tT.jTdist))"

    # Kinetics
    admdist = kinetics.dmdist[end, :]
    D0Lmean = mean(admdist .|> x-> x.D0L)
    @test 0 < D0Lmean
    @info "Mean apatite D0L: $D0Lmean"
    EaLmean = mean(admdist .|> x-> x.EaL)
    @test 0 < EaLmean 
    @info "Mean apatite EaL: $EaLmean"
    EaTrapmean = mean(admdist .|> x-> x.EaTrap)
    @test 0 < EaTrapmean 
    @info "Mean apatite EaTrap: $EaTrapmean"
    rmr0mean = mean(admdist .|> x-> x.rmr0)
    @test 0 < rmr0mean < 1
    @info "Mean apatite rmr0: $rmr0mean"
    zdmdist = kinetics.dmdist[1, :]
    DzD0mean = mean(zdmdist .|> x-> x.DzD0)
    @test 0 < DzD0mean
    @info "Mean zircon DzD0: $DzD0mean"
    DzEamean = mean(zdmdist .|> x-> x.DzEa)
    @test 0 < DzEamean
    @info "Mean zircon DzEa: $DzEamean"
    DN17D0mean = mean(zdmdist .|> x-> x.DN17D0)
    @test 0 < DN17D0mean
    @info "Mean zircon DN17D0: $DN17D0mean"
    DN17Eamean = mean(zdmdist .|> x-> x.DN17Ea)
    @test 0 < DN17Eamean
    @info "Mean zircon DN17Ea: $DN17Eamean"
    rminmean = mean(zdmdist .|> x-> x.rmin)
    @test 0 < rminmean < 1
    @info "Mean zircon rmin: $rminmean"

    ## --- Add dynamic jumping, IGB partitioning, and a constraint box
    params = (params...,
        dynamicjumping=true,
        partitiondaughter=true,
    )
    unconf = Constraint(
        agedist = [Uniform(500,580),],  # [Ma] Age distribution
        Tdist =   [   Uniform(0,50),],  # [C] Temperature distribution
    )

    @time "\nMCMC with Detail interval & dynamicjumping" tT = MCMC(data, params, boundary, unconf, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 175
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -450 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test params.dt < mean(tT.jtdist) < params.tinit
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test 0 < mean(tT.jTdist) < params.Tinit
    @info "Mean σjT: $(mean(tT.jTdist))"

## --- MCMC_varkinetics with Detail interval & dynamicjumping

    @time "\nMCMC_varkinetics with Detail interval & dynamicjumping" tT, kinetics = MCMC_varkinetics(data, params, boundary, unconf, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 200
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -450 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test params.dt < mean(tT.jtdist) < params.tinit
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test 0 < mean(tT.jTdist) < params.Tinit
    @info "Mean σjT: $(mean(tT.jTdist))"

    # Kinetics
    admdist = kinetics.dmdist[end, :]
    D0Lmean = mean(admdist .|> x-> x.D0L)
    @test 0 < D0Lmean
    @info "Mean apatite D0L: $D0Lmean"
    EaLmean = mean(admdist .|> x-> x.EaL)
    @test 0 < EaLmean 
    @info "Mean apatite EaL: $EaLmean"
    EaTrapmean = mean(admdist .|> x-> x.EaTrap)
    @test 0 < EaTrapmean 
    @info "Mean apatite EaTrap: $EaTrapmean"
    rmr0mean = mean(admdist .|> x-> x.rmr0)
    @test 0 < rmr0mean < 1
    @info "Mean apatite rmr0: $rmr0mean"
    zdmdist = kinetics.dmdist[1, :]
    DzD0mean = mean(zdmdist .|> x-> x.DzD0)
    @test 0 < DzD0mean
    @info "Mean zircon DzD0: $DzD0mean"
    DzEamean = mean(zdmdist .|> x-> x.DzEa)
    @test 0 < DzEamean
    @info "Mean zircon DzEa: $DzEamean"
    DN17D0mean = mean(zdmdist .|> x-> x.DN17D0)
    @test 0 < DN17D0mean
    @info "Mean zircon DN17D0: $DN17D0mean"
    DN17Eamean = mean(zdmdist .|> x-> x.DN17Ea)
    @test 0 < DN17Eamean
    @info "Mean zircon DN17Ea: $DN17Eamean"
    rminmean = mean(zdmdist .|> x-> x.rmin)
    @test 0 < rminmean < 1
    @info "Mean zircon rmin: $rminmean"

## --- Test conversion of t-T images 

    ttimage, xq, yq = image_from_paths(tT, xresolution=200, yresolution=100, yrange=boundary.T₀)
    @test xq == 8.75:17.5:3491.25
    @test yq == 2.0:4.0:398.0
    @test ttimage isa Matrix
    @test axes(ttimage, 1) == 1:100
    @test axes(ttimage, 2) == 1:200

    ttimage2, xq2, yq2 = image_from_paths!(tT, xresolution=200, yresolution=100, yrange=boundary.T₀)
    @test xq2 == 8.75:17.5:3491.25
    @test yq2 == 2.0:4.0:398.0
    @test ttimage2 isa Matrix
    @test axes(ttimage2, 1) == 1:100
    @test axes(ttimage2, 2) == 1:200

    @test ttimage == ttimage2

## --- End of File