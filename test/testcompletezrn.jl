
# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesotazrn.csv")
ds = importdataset(datapath, ',', importas=:Tuple);

using LinearAlgebra
BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

## --- Prepare problem

model = (
    burnin = 350,               # [n] How long should we wait for MC to converge (become stationary)
    nsteps = 250,               # [n] How many steps of the Markov chain should we run after burn-in?
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
    Tr = 250.,                  # [C] Residence temperature of initial proposal
    simplified = false,         # Prefer simpler tT paths?
    dynamicsigma = true,        # Update model uncertainties?
    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
    σmodel = 25.0,              # [Ma] assumed model uncertainty (resampled if dynamicsigma)
    σannealing = 35.0,          # [Ma] initial annealing uncertainty
    λannealing = 10 ./ 200,     # [1/n] annealing decay
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

# Crystallization ages and start time
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

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 150
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
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

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 150
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
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

@test isa(tT.resultdist, AbstractMatrix)
abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 150
@info "Mean absolute error: $abserr"

@test isa(tT.lldist, AbstractVector)
llmean = mean(tT.lldist)
@test -300 < llmean < 0
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