
# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesota.csv")
ds = importdataset(datapath, ',', importas=:Tuple);

using LinearAlgebra
BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

## --- Prepare problem

model = (
    nsteps = 500, # How many steps of the Markov chain should we run?
    burnin = 300, # How long should we wait for MC to converge (become stationary)
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
    halfwidth = ds.Halfwidth_um,            # Crystal half-width, in microns
    U = ds.U238_ppm,                        # U concentration, in PPM
    Th = ds.Th232_ppm,                      # Th-232 concentration, in PPM
    Sm = ds.Sm147_ppm,                      # Sm-147 concentration, in PPM (optional)
    HeAge = ds.HeAge_Ma_raw,                # He age, in Ma
    HeAge_sigma = ds.HeAge_Ma_sigma_10pct,  # He age uncertainty (1-sigma), in Ma
    crystAge = ds.CrystAge_Ma,              # Crystallization age, in Ma
    mineral = ds.Mineral                    # zircon or apatite
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


## --- Invert for maximum likelihood t-T path

# Run Markov Chain
@time "\nCompiling MCMC" MCMC(data, model, boundary, unconf)
@time "\nRunning MCMC" result = MCMC(data, model, boundary, unconf)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"


@test mean(result.jtdist[model.burnin:end]) ≈ model.tinit/60
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test mean(result.jTdist[model.burnin:end]) ≈ model.Tinit/60
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"

## --- As above, but with variable kinetic parameters

# Run Markov Chain
@time "\nCompiling MCMC_varkinetics" MCMC_varkinetics(data, model, boundary, unconf)
@time "\nRunning MCMC_varkinetics" result = MCMC(data, model, boundary, unconf)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"


@test mean(result.jtdist[model.burnin:end]) ≈ model.tinit/60
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test mean(result.jTdist[model.burnin:end]) ≈ model.Tinit/60
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"

## ---
detail = DetailInterval(
    agemin = 0, # Youngest end of detail interval
    agemax = 541, # Oldest end of detail interval
    minpoints = 5, # Minimum number of points in detail interval
)
@time "\nMCMC with Detail interval" result = MCMC(data, model, boundary, unconf, detail)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"

@test mean(result.jtdist[model.burnin:end]) ≈ model.tinit/60
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test mean(result.jTdist[model.burnin:end]) ≈ model.Tinit/60
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"

## ---

@time "\nMCMC_varkinetics with Detail interval" result = MCMC_varkinetics(data, model, boundary, unconf, detail)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"


@test mean(result.jtdist[model.burnin:end]) ≈ model.tinit/60
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test mean(result.jTdist[model.burnin:end]) ≈ model.Tinit/60
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"

## ---
model = (model...,
    dynamicjumping=true
)
@time "\nMCMC with Detail interval & dynamicjumping" result = MCMC(data, model, boundary, unconf, detail)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"

@test model.dt < mean(result.jtdist[model.burnin:end]) < model.tinit
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test 0 < mean(result.jTdist[model.burnin:end]) < model.Tinit
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"

## ---

@time "\nMCMC_varkinetics with Detail interval & dynamicjumping" result = MCMC_varkinetics(data, model, boundary, unconf, detail)

@test isa(result.Tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.Tinit
@test nanminimum(result.Tpointdist) >= model.Tnow

@test isa(result.tpointdist, AbstractMatrix)
@test nanmaximum(result.Tpointdist) <= model.tinit
@test nanminimum(result.Tpointdist) >= 0

@test isa(result.HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(result.HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(result.lldist, AbstractVector)
llmean = mean(@view(result.lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(result.acceptancedist, AbstractVector{Bool})
@test isapprox(mean(result.acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(result.acceptancedist[model.burnin:end]))"

@test isa(result.ndist, AbstractVector{Int})
@test minimum(result.ndist) >= 0
@test maximum(result.ndist) <= model.maxpoints
@info "Mean npoints: $(mean(result.ndist[model.burnin:end]))"

@test model.dt < mean(result.jtdist[model.burnin:end]) < model.tinit
@info "Mean jₜ: $(mean(result.jtdist[model.burnin:end]))"

@test 0 < mean(result.jTdist[model.burnin:end]) < model.Tinit
@info "Mean jT: $(mean(result.jTdist[model.burnin:end]))"
