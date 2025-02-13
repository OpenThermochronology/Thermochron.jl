
# Read in data from file using StatGeochem
datapath = joinpath("..", "examples")
ds = importdataset(joinpath(datapath, "minnesota.csv"), ',', importas=:Tuple);

using LinearAlgebra
BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

## --- Prepare problem

model = (
    burnin = 350,               # [n] How long should we wait for MC to converge (become stationary)
    nsteps = 350,               # [n] How many steps of the Markov chain should we run after burn-in?
    dr = 1.0,                   # [μ] Radius step size
    dt = 10.0,                  # [Ma] time step size
    dTmax = 25.0,               # [C/step] Maximum reheating/burial per model timestep
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
    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
    σmodel = 25.0,              # [Ma] assumed model uncertainty (resampled if dynamicsigma)
    T0annealing = 1,            # [unitless] initial annealing "temperature"
)

# Populate data NamedTuple from imported dataset
data = (
    halfwidth = ds.halfwidth_um,            # Crystal half-width, in microns
    U = ds.U238_ppm,                        # U concentration, in PPM
    Th = ds.Th232_ppm,                      # Th-232 concentration, in PPM
    Sm = ds.Sm147_ppm,                      # Sm-147 concentration, in PPM (optional)
    HeAge = ds.raw_He_age_Ma,               # He age, in Ma
    HeAge_sigma = ds.raw_He_age_Ma.*0.1,    # He age uncertainty (1-sigma), in Ma
    crystAge = ds.crystallization_age_Ma,   # Crystallization age, in Ma
    mineral = ds.mineral                    # zircon or apatite
)

# Crystallization ages and start time
tinit = ceil(maximum(data.crystAge)/model.dt) * model.dt
model = (model...,
    tinit = tinit,
    agesteps = (tinit-model.dt/2 : -model.dt : 0+model.dt/2),
    tsteps = (0+model.dt/2 : model.dt : tinit-model.dt/2),
)

# Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
boundary = Boundary(
    agepoints = [model.tnow, model.tinit],   # [Ma] Final and initial time
    T₀ = [model.Tnow, model.Tinit],          # [C] Final and initial temperature
    ΔT = [model.ΔTnow, model.ΔTinit],        # [C] Final and initial temperature range (positive or negative)
    tboundary = :reflecting, # Reflecting time boundary conditions
    Tboundary = :reflecting, # Reflecting temperature boundary conditions
)

# Default: No unconformity is imposed
unconf = Constraint()

## --- Test generation of Chronometer objects

# Modern input format, generic
dsg = importdataset(joinpath(datapath, "generic.csv"), ',', importas=:Tuple);
chrons = chronometers(dsg, model)
@test chrons isa Vector{<:Chronometer}
@test length(chrons) == 18
@test count(x->isa(x,GenericHe), chrons) == 4
@test count(x->isa(x,GenericAr), chrons) == 2
@test count(x->isa(x,ZirconHe), chrons) == 2
@test count(x->isa(x,ApatiteHe), chrons) == 2
@test count(x->isa(x,ZirconFT), chrons) == 1
@test count(x->isa(x,MonaziteFT), chrons) == 1
@test count(x->isa(x,ApatiteFT), chrons) == 1
@test count(x->isa(x,ZirconTrackLength), chrons) == 1
@test count(x->isa(x,MonaziteTrackLength), chrons) == 1
@test count(x->isa(x,ApatiteTrackLength), chrons) == 3
@test get_age(chrons) ≈ [150.37, 263.92, 150.37, 263.92, 917.84, 1023.73, 380., 380., 120., 120., 680., 300.,100.] 
@test get_age_sigma(chrons) ≈ [5,5,5,5,5,5,5,5,5,5,5,5,5]

dt = 10.0
tsteps = (0+dt/2 : dt : 3000-dt/2)
Tsteps = range(650, 0, length=length(tsteps))

calc = zeros(length(chrons))
calcuncert = zeros(length(chrons))
@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Yamada2007PC(), Jones2021FA(), Ketcham2007FC(); trackhist=false) ≈ -95.97338155773869
@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Yamada2007PC(), Jones2021FA(), Ketcham2007FC(); trackhist=true) ≈ -97.47230160573113 
@test round.(calc[1:end-5], sigdigits=7) ≈ [138.4124, 232.8114, 144.2487, 233.9706, 902.567, 1010.98, 386.8558, 388.6112, 122.8127, 130.7147, 688.0081, 304.6573, 95.84216] 
@test round.(calc[end-4:end], sigdigits=3) ≈ [8, 14.3, 14.3, 14.3, 6] atol = 9
@test calcuncert[1:end-5] ≈ zeros(length(chrons)-5)
@test calcuncert[end-4:end] ≈ [1.8364436281549124, 1.1785910438098226, 1.1389520917140208, 1.2018361658877996, 0.6070538659171328]

@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Guenthner2013FC(), Jones2021FA(), Ketcham2007FC(); trackhist=false) ≈ -3799.27575648043
@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Guenthner2013FC(), Jones2021FA(), Ketcham2007FC(); trackhist=true) ≈ -3800.7746765284223
@test round.(calc[1:end-5], sigdigits=7) ≈ [138.4124, 232.8114, 144.2487, 233.9706, 902.567, 1010.98, 386.8558, 388.6112, 122.8127, 130.7147, 1110.379, 304.6573, 95.84216]
@test round.(calc[end-4:end], sigdigits=3) ≈ [8, 14.3, 14.3, 14.3, 6] atol = 9
@test calcuncert[1:end-5] ≈ zeros(length(chrons)-5)
@test calcuncert[end-4:end] ≈ [1.8368172844202661, 1.1785910438098226, 1.1389520917140208, 1.2018361658877996, 0.6070538659171328]

# Modern input format, Minnesota dataset
chrons = chronometers(ds, model)
@test chrons isa Vector{<:Chronometer}
@test length(chrons) == length(ds.mineral)
@test get_age(chrons) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11, 234, 98, 233, 339, 378, 258, 158, 269, 313, 309, 392]
@test get_age_sigma(chrons) ≈ [15.0, 51.0, 13.0, 12.6, 18.4, 10.0, 12.0, 10.0, 43.0, 7.0, 7.0, 5.2, 4.6, 4.6, 6.0, 5.0, 4.0, 4.0, 2.0, 1.1, 1.1, 0.8, 0.2, 9.4, 2.0, 5.0, 6.0, 9.3, 5.0, 3.3, 8.1, 5.8, 6.0, 8.0]

tsteps = model.tsteps
Tsteps = range(model.Tinit, model.Tnow, length=length(tsteps))

calc = zeros(length(chrons))
calcuncert = zeros(length(chrons))
@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Yamada2007PC(), Jones2021FA(), Ketcham2007FC(); trackhist=false) ≈ -26867.66124989401
@test Thermochron.model!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), Yamada2007PC(), Jones2021FA(), Ketcham2007FC(); trackhist=true) ≈ -26867.66124989401

@test round.(calc, sigdigits=5) ≈ [1313.9, 1320.7, 1185.0, 1243.0, 1216.1, 1335.4, 1141.7, 1094.4, 1170.2, 923.8, 723.59, 201.76, 429.67, 95.576, 259.05, 419.15, 2.9065, 6.1464, 0.00063415, 27.545, 0.007082, 55.056, 2.0682, 174.81, 283.3, 287.74, 266.09, 240.19, 267.84, 244.36, 274.74, 328.26, 322.88, 352.43]
@test calcuncert ≈ zeros(length(chrons))

# Test an individual zircon
@test first(calc) ≈ modelage(first(chrons), Tsteps, ZRDAAM())
# Test an individual apatite
@test last(calc) ≈ modelage(last(chrons), Tsteps, RDAAM())

# Test empirical uncertainty estimation
σcalc = zeros(length(chrons))
empiricaluncertainty!(σcalc, chrons, ZirconHe)
@test get_age(chrons, ZirconHe) ≈ [770, 659, 649, 638, 620, 557, 545, 500, 493, 357, 329, 253, 241, 225, 225, 217, 193, 190, 72, 57, 42, 29, 11]
σtotal = sqrt.(σcalc[isa.(chrons, ZirconHe)].^2 + get_age_sigma(chrons, ZirconHe).^2)
@test σtotal ≈ [71.76511200022767, 88.69939042178478, 76.20513740085536, 49.74341557807492, 52.05638860615669, 58.80850040722076, 57.32853093748527, 97.59191128458991, 68.09272937102548, 82.41401338162268, 104.5147581608337, 72.10558616290515, 79.47676757701001, 68.01965448205827, 72.8583946436481, 73.94974380762895, 63.02840819654363, 63.478989858090436, 45.255278015066644, 66.3078727295436, 52.90379008918621, 69.14186488208391, 62.394435852786344]

set_age!(chrons, zeros(23), ZirconHe)
set_age_sigma!(chrons, zeros(23), ZirconHe)
@test get_age(chrons, ZirconHe) ≈ zeros(23)
@test get_age_sigma(chrons, ZirconHe) ≈ zeros(23)

# Legacy input format
chrons2 = chronometers(data, model)
@test typeof.(chrons) == typeof.(chrons2)
Thermochron.model!(calc, calcuncert, chrons2, Tsteps, ZRDAAM(), RDAAM(), Yamada2007PC(), Jones2021FA(), Ketcham2007FC())
@test round.(calc, sigdigits=5) ≈ [1313.9, 1320.7, 1185.0, 1243.0, 1216.1, 1335.4, 1141.7, 1094.4, 1170.2, 923.8, 723.59, 201.76, 429.67, 95.576, 259.05, 419.15, 2.9065, 6.1464, 0.00063415, 27.545, 0.007082, 55.056, 2.0682, 174.81, 283.3, 287.74, 266.09, 240.19, 267.84, 244.36, 274.74, 328.26, 322.88, 352.43]
@test calcuncert ≈ zeros(length(chrons))
# println(round.(calc, sigdigits=5))

## --- Invert for maximum likelihood t-T path

# Run Markov Chain
@time "\nCompiling MCMC" MCMC(data, model, boundary, unconf)
@time "\nRunning MCMC" tT = MCMC(data, model, boundary, unconf)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"


@test mean(tT.jtdist) ≈ model.tinit/60
@info "Mean σjt: $(mean(tT.jtdist))"

@test mean(tT.jTdist) ≈ model.Tinit/60
@info "Mean σjT : $(mean(tT.jTdist))"

## --- As above, but with variable kinetic parameters

# Run Markov Chain
@time "\nCompiling MCMC_varkinetics" MCMC_varkinetics(data, model, boundary, unconf)
@time "\nRunning MCMC_varkinetics" tT, kinetics = MCMC_varkinetics(data, model, boundary, unconf)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"


@test mean(tT.jtdist) ≈ model.tinit/60
@info "Mean σjt: $(mean(tT.jtdist))"

@test mean(tT.jTdist) ≈ model.Tinit/60
@info "Mean σjT: $(mean(tT.jTdist))"

# Kinetics
D0Lmean = mean(kinetics.admdist .|> x-> x.D0L)
@test 0 < D0Lmean
@info "Mean apatite D0L: $D0Lmean"
EaLmean = mean(kinetics.admdist .|> x-> x.EaL)
@test 0 < EaLmean 
@info "Mean apatite EaL: $EaLmean"
EaTrapmean = mean(kinetics.admdist .|> x-> x.EaTrap)
@test 0 < EaTrapmean 
@info "Mean apatite EaTrap: $EaTrapmean"
rmr0mean = mean(kinetics.admdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean apatite rmr0: $rmr0mean"
DzD0mean = mean(kinetics.zdmdist .|> x-> x.DzD0)
@test 0 < DzD0mean
@info "Mean zircon DzD0: $DzD0mean"
DzEamean = mean(kinetics.zdmdist .|> x-> x.DzEa)
@test 0 < DzEamean
@info "Mean zircon DzEa: $DzEamean"
DN17D0mean = mean(kinetics.zdmdist .|> x-> x.DN17D0)
@test 0 < DN17D0mean
@info "Mean zircon DN17D0: $DN17D0mean"
DN17Eamean = mean(kinetics.zdmdist .|> x-> x.DN17Ea)
@test 0 < DN17Eamean
@info "Mean zircon DN17Ea: $DN17Eamean"
rmr0mean = mean(kinetics.zdmdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean zircon rmr0: $rmr0mean"

## ---
detail = DetailInterval(
    agemin = 0, # Youngest end of detail interval
    agemax = 1000, # Oldest end of detail interval
    minpoints = 3, # Minimum number of points in detail interval
)
@time "\nMCMC with Detail interval" tT = MCMC(data, model, boundary, unconf, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"

@test mean(tT.jtdist) ≈ model.tinit/60
@info "Mean σjt: $(mean(tT.jtdist))"

@test mean(tT.jTdist) ≈ model.Tinit/60
@info "Mean σjT: $(mean(tT.jTdist))"

## ---

@time "\nMCMC_varkinetics with Detail interval" tT, kinetics = MCMC_varkinetics(data, model, boundary, unconf, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"


@test mean(tT.jtdist) ≈ model.tinit/60
@info "Mean σjt: $(mean(tT.jtdist))"

@test mean(tT.jTdist) ≈ model.Tinit/60
@info "Mean σjT: $(mean(tT.jTdist))"

# Kinetics
D0Lmean = mean(kinetics.admdist .|> x-> x.D0L)
@test 0 < D0Lmean
@info "Mean apatite D0L: $D0Lmean"
EaLmean = mean(kinetics.admdist .|> x-> x.EaL)
@test 0 < EaLmean 
@info "Mean apatite EaL: $EaLmean"
EaTrapmean = mean(kinetics.admdist .|> x-> x.EaTrap)
@test 0 < EaTrapmean 
@info "Mean apatite EaTrap: $EaTrapmean"
rmr0mean = mean(kinetics.admdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean apatite rmr0: $rmr0mean"
DzD0mean = mean(kinetics.zdmdist .|> x-> x.DzD0)
@test 0 < DzD0mean
@info "Mean zircon DzD0: $DzD0mean"
DzEamean = mean(kinetics.zdmdist .|> x-> x.DzEa)
@test 0 < DzEamean
@info "Mean zircon DzEa: $DzEamean"
DN17D0mean = mean(kinetics.zdmdist .|> x-> x.DN17D0)
@test 0 < DN17D0mean
@info "Mean zircon DN17D0: $DN17D0mean"
DN17Eamean = mean(kinetics.zdmdist .|> x-> x.DN17Ea)
@test 0 < DN17Eamean
@info "Mean zircon DN17Ea: $DN17Eamean"
rmr0mean = mean(kinetics.zdmdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean zircon rmr0: $rmr0mean"

## --- Add dynamic jumping and a constraint box
model = (model...,
    dynamicjumping=true
)
unconf = Constraint(
    agedist = [Uniform(500,580),],  # [Ma] Age distribution
    Tdist =   [   Uniform(0,50),],  # [C] Temperature distribution
)

@time "\nMCMC with Detail interval & dynamicjumping" tT = MCMC(data, model, boundary, unconf, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"

@test model.dt < mean(tT.jtdist) < model.tinit
@info "Mean σjt: $(mean(tT.jtdist))"

@test 0 < mean(tT.jTdist) < model.Tinit
@info "Mean σjT: $(mean(tT.jTdist))"

## ---

@time "\nMCMC_varkinetics with Detail interval & dynamicjumping" tT, kinetics = MCMC_varkinetics(data, model, boundary, unconf, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
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
@test minimum(tT.ndist) >= 0
@test maximum(tT.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(tT.ndist))"

@test model.dt < mean(tT.jtdist) < model.tinit
@info "Mean σjt: $(mean(tT.jtdist))"

@test 0 < mean(tT.jTdist) < model.Tinit
@info "Mean σjT: $(mean(tT.jTdist))"

# Kinetics
D0Lmean = mean(kinetics.admdist .|> x-> x.D0L)
@test 0 < D0Lmean
@info "Mean apatite D0L: $D0Lmean"
EaLmean = mean(kinetics.admdist .|> x-> x.EaL)
@test 0 < EaLmean 
@info "Mean apatite EaL: $EaLmean"
EaTrapmean = mean(kinetics.admdist .|> x-> x.EaTrap)
@test 0 < EaTrapmean 
@info "Mean apatite EaTrap: $EaTrapmean"
rmr0mean = mean(kinetics.admdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean apatite rmr0: $rmr0mean"
DzD0mean = mean(kinetics.zdmdist .|> x-> x.DzD0)
@test 0 < DzD0mean
@info "Mean zircon DzD0: $DzD0mean"
DzEamean = mean(kinetics.zdmdist .|> x-> x.DzEa)
@test 0 < DzEamean
@info "Mean zircon DzEa: $DzEamean"
DN17D0mean = mean(kinetics.zdmdist .|> x-> x.DN17D0)
@test 0 < DN17D0mean
@info "Mean zircon DN17D0: $DN17D0mean"
DN17Eamean = mean(kinetics.zdmdist .|> x-> x.DN17Ea)
@test 0 < DN17Eamean
@info "Mean zircon DN17Ea: $DN17Eamean"
rmr0mean = mean(kinetics.zdmdist .|> x-> x.rmr0)
@test 0 < rmr0mean < 1
@info "Mean zircon rmr0: $rmr0mean"

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