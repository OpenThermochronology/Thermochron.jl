
# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesota.csv")
ds = importdataset(datapath, ',', importas=:Tuple);

using LinearAlgebra
BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

## --- Prepare problem

model = (
    burnin = 350, # How long should we wait for MC to converge (become stationary)
    nsteps = 250, # How many steps of the Markov chain should we run after burn-in?
    dr = 1.0,    # Radius step, in microns
    dt = 10.0,   # time step size in Myr
    dTmax = 25.0, # Maximum reheating/burial per model timestep
    Tinit = 400.0, # initial model temperature (in C) (i.e., crystallization temperature)
    ΔTinit = -50.0, # Tinit can vary from Tinit to Tinit+ΔTinit
    Tnow = 0.0, # Current surface temperature (in C)
    ΔTnow = 10.0, # Tnow may vary from Tnow to Tnow+ΔTnow
    tnow = 0.0 ,  # Today
    tinitMax = 4000.0, # Ma -- forbid anything older than this
    minpoints = 1,  # Minimum allowed number of t-T points
    maxpoints = 40, # Maximum allowed number of t-T points
    npoints = 5, # Initial number of t-T points
    Tr = 250., # Residence temperature of initial proposal
    simplified = false, # Prefer simpler tT paths?
    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
    σmodel = 25.0, # Ma
    σannealing = 35.0, # initial annealing uncertainty [Ma]
    λannealing = 10 ./ 200 # annealing decay [1/n]
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

# Sort out crystallization ages and start time
map!(x->min(x, model.tinitMax), data.crystAge, data.crystAge)
tinit = ceil(maximum(data.crystAge)/model.dt) * model.dt
model = (model...,
    tinit = tinit,
    agesteps = Array{Float64}(tinit-model.dt/2 : -model.dt : 0+model.dt/2),
    tsteps = Array{Float64}(0+model.dt/2 : model.dt : tinit-model.dt/2),
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

# Modern input format
chrons = chronometers(ds, model)
@test chrons isa Vector{<:Chronometer}
@test length(chrons) == length(ds.mineral)

tsteps = model.tsteps
Tsteps = range(model.Tinit, model.Tnow, length=length(tsteps))

calc = zeros(length(chrons))
calcuncert = zeros(length(chrons))
Thermochron.modelages!(calc, calcuncert, chrons, Tsteps, ZRDAAM(), RDAAM(), FCKetcham2007)
@test round.(calc, sigdigits=5) ≈ [1313.9, 1320.7, 1185.0, 1243.0, 1216.1, 1335.4, 1141.7, 1094.4, 1170.2, 923.8, 723.59, 201.76, 429.67, 95.576, 259.05, 419.15, 2.9065, 6.1464, 0.00063415, 27.545, 0.007082, 55.056, 2.0682, 174.81, 283.3, 287.74, 266.14, 240.19, 267.84, 244.37, 274.74, 328.26, 322.88, 352.43]
@test calcuncert ≈ zeros(length(chrons))

# Test an individual zircon
@test first(calc) ≈ modelage(first(chrons), Tsteps, ZRDAAM())
# Test an individual apatite
@test last(calc) ≈ modelage(last(chrons), Tsteps, RDAAM())

# Legacy input format
chrons2 = chronometers(data, model)
@test typeof.(chrons) == typeof.(chrons2)
Thermochron.modelages!(calc, calcuncert, chrons2, Tsteps, ZRDAAM(), RDAAM(), FCKetcham2007)
@test round.(calc, sigdigits=5) ≈ [1313.9, 1320.7, 1185.0, 1243.0, 1216.1, 1335.4, 1141.7, 1094.4, 1170.2, 923.8, 723.59, 201.76, 429.67, 95.576, 259.05, 419.15, 2.9065, 6.1464, 0.00063415, 27.545, 0.007082, 55.056, 2.0682, 174.81, 283.3, 287.74, 266.14, 240.19, 267.84, 244.37, 274.74, 328.26, 322.88, 352.43]
@test calcuncert ≈ zeros(length(chrons))

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
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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
    agemax = 541, # Oldest end of detail interval
    minpoints = 5, # Minimum number of points in detail interval
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
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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
@test 0 < abserr < 150
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -420 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(tT.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(tT.acceptancedist), 0.5, atol=0.4)
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