
# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesotazrn.csv")
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
    npoints = 5, # Number of initial t-T points
    Tr = 250., # Residence temperature in initial proposal
    simplified = false, # Prefer simpler tT paths?
    tboundary = :reflecting, # Reflecting time boundary conditions
    Tboundary = :reflecting, # Reflecting temperature boundary conditions
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
    Th = ds.Th232_ppm,                      # Th concentration, in PPM
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

# Default: No constraint is imposed
constraint = Constraint()


## --- Invert for maximum likelihood t-T path

# Run Markov Chain
@time "\nCompiling MCMC" MCMC(data, model, boundary, constraint)
@time "\nRunning MCMC" tT = MCMC(data, model, boundary, constraint)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.HeAgedist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 50
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -200 < llmean < 0
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
detail = DetailInterval(
    agemin = 0, # Youngest end of detail interval
    agemax = 541, # Oldest end of detail interval
    minpoints = 5, # Minimum number of points in detail interval
)
@time "\nMCMC with Detail interval" tT = MCMC(data, model, boundary, constraint, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.HeAgedist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 50
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -200 < llmean < 0
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
model = (model...,
    dynamicjumping=true
)
@time "\nMCMC with Detail interval & dynamicjumping" tT = MCMC(data, model, boundary, constraint, detail)

@test isa(tT.Tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.Tinit
@test nanminimum(tT.Tpointdist) >= model.Tnow

@test isa(tT.tpointdist, AbstractMatrix)
@test nanmaximum(tT.Tpointdist) <= model.tinit
@test nanminimum(tT.Tpointdist) >= 0

@test isa(tT.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.HeAgedist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 50
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -200 < llmean < 0
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