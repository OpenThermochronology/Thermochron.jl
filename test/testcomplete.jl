
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
    tinitMax = 4000.0, # Ma -- forbid anything older than this
    minpoints = 1,  # Minimum allowed number of t-T points
    maxpoints = 40, # Maximum allowed number of t-T points
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
map!(x->max(x, model.tinitMax), data.crystAge, data.crystAge)
tinit = ceil(maximum(data.crystAge)/model.dt) * model.dt
model = (model...,
    tinit = tinit,
    agesteps = Array{Float64}(tinit-model.dt/2 : -model.dt : 0+model.dt/2),
    tsteps = Array{Float64}(0+model.dt/2 : model.dt : tinit-model.dt/2),
)

# Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
boundary = Boundary(
    agepoints = Float64[model.Tnow, model.tinit],  # Ma
    Tpoints = Float64[model.Tnow, model.Tinit],    # Degrees C
    T₀ = Float64[model.Tnow, model.Tinit],
    ΔT = Float64[model.ΔTnow, model.ΔTinit],
)

# Default: No unconformity is imposed
unconf = Unconformity()


## --- Invert for maximum likelihood t-T path

# This is where the "transdimensional" part comes in
agepoints = Array{Float64}(undef, model.maxpoints+1) # Array of fixed size to hold all optional age points
Tpoints = Array{Float64}(undef, model.maxpoints+1) # Array of fixed size to hold all optional age points
# Fill some intermediate points to give the MCMC something to work with
Tr = 250 # Residence temperature
npoints = 5
agepoints[1:npoints] .= (model.tinit/30,model.tinit/4,model.tinit/2,model.tinit-model.tinit/4,model.tinit-model.tinit/30) # Ma
Tpoints[1:npoints] .= Tr  # Degrees C

# Run Markov Chain
@time "Compiling MCMC_vartcryst" MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, boundary, unconf)
@time "Running MCMC_vartcryst" (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist) = MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, boundary, unconf)

@test isa(Tpointdist, AbstractMatrix)
@test maximum(Tpointdist) <= model.Tinit
@test minimum(Tpointdist) >= model.Tnow

@test isa(tpointdist, AbstractMatrix)
@test maximum(tpointdist) <= model.tinit
@test minimum(tpointdist) >= 0

@test isa(HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(lldist, AbstractVector)
llmean = mean(@view(lldist[model.burnin:end]))
@test -300 < llmean < 0
@info "Mean ll: $llmean"

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(acceptancedist[model.burnin:end]))"

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxpoints
@info "Mean npoints: $(mean(ndist[model.burnin:end]))"

@info "Mean σⱼₜ: $(mean(σⱼtdist[model.burnin:end]))"
@info "Mean σⱼT: $(mean(σⱼTdist[model.burnin:end]))"

## ---
detail = DetailInterval(
    agemin = 0, # Youngest end of detail interval
    agemax = 541, # Oldest end of detail interval
    minpoints = 5, # Minimum number of points in detail interval
)
@time "MCMC_vartcryst with Detail interval" (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist) = MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, boundary, unconf, detail)

@test isa(Tpointdist, AbstractMatrix)
@test maximum(Tpointdist) <= model.Tinit
@test minimum(Tpointdist) >= model.Tnow

@test isa(tpointdist, AbstractMatrix)
@test maximum(tpointdist) <= model.tinit
@test minimum(tpointdist) >= 0

@test isa(HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(lldist, AbstractVector)
llmean = mean(@view(lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(acceptancedist[model.burnin:end]))"

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxpoints
@info "Mean npoints: $(mean(ndist[model.burnin:end]))"

@info "Mean σⱼₜ: $(mean(σⱼtdist[model.burnin:end]))"
@info "Mean σⱼT: $(mean(σⱼTdist[model.burnin:end]))"

## ---
model = (model...,
    dynamicjumping=true
)
@time "MCMC_vartcryst with Detail interval & dynamicjumping" (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist) = MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, boundary, unconf, detail)

@test isa(Tpointdist, AbstractMatrix)
@test maximum(Tpointdist) <= model.Tinit
@test minimum(Tpointdist) >= model.Tnow

@test isa(tpointdist, AbstractMatrix)
@test maximum(tpointdist) <= model.tinit
@test minimum(tpointdist) >= 0

@test isa(HeAgedist, AbstractMatrix)
abserr = abs(sum(nanmean(HeAgedist[:,model.burnin:end], dims=2) - data.HeAge)/length(data.HeAge))
@test 0 < abserr < 100
@info "Mean absolute error: $abserr"

@test isa(lldist, AbstractVector)
llmean = mean(@view(lldist[model.burnin:end]))
@test -200 < llmean < 0
@info "Mean ll: $llmean"

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.5, atol=0.4)
@info "Mean acceptance rate: $(mean(acceptancedist[model.burnin:end]))"

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxpoints
@info "Mean npoints: $(mean(ndist[model.burnin:end]))"

@info "Mean σⱼₜ: $(mean(σⱼtdist[model.burnin:end]))"
@info "Mean σⱼT: $(mean(σⱼTdist[model.burnin:end]))"
