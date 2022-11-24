using Statistics
using StatGeochem

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
    dTmax = 10.0, # Maximum reheating/burial per model timestep
    TInit = 400.0, # Initial model temperature (in C) (i.e., crystallization temperature)
    ΔTInit = -50.0, # TInit can vary from TInit to TInit+ΔTinit
    TNow = 0.0, # Current surface temperature (in C)
    ΔTNow = 10.0, # TNow may vary from TNow to TNow+ΔTNow
    tInitMax = 4000.0, # Ma -- forbid anything older than this
    minPoints = 1,  # Minimum allowed number of t-T points
    maxPoints = 40, # Maximum allowed number of t-T points
    simplified = false, # Prefer simpler tT paths?
    # Diffusion parameters
    DzEa = 165.0, # kJ/mol
    DzD0 = 193188.0, # cm^2/sec
    DN17Ea = 71.0, # kJ/mol
    DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
    σModel = 25.0, # Ma
    σAnnealing = 35.0, # Initial annealing uncertainty [Ma]
    λAnnealing = 10 ./ 200 # Annealing decay [1/n]
)

# Populate data NamedTuple from imported dataset
data = (
    halfwidth = ds.Halfwidth_um,            # Crystal half-width, in microns
    U = ds.U238_ppm,                        # U concentration, in PPM
    Th = ds.Th232_ppm,                      # Th concentration, in PPM
    HeAge = ds.HeAge_Ma_raw,                # He age, in Ma
    HeAge_sigma = ds.HeAge_Ma_sigma_raw,    # He age uncertainty (1-sigma), in Ma
    CrystAge = ds.CrystAge_Ma,              # Crystallization age, in Ma
)

# Sort out crystallization ages and start time
map!(x->max(x, model.tInitMax), data.CrystAge, data.CrystAge)
model = (model...,
    tInit = ceil(maximum(data.CrystAge)/model.dt) * model.dt,
)

# Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
boundary = (
    agePoints = Float64[model.TNow, model.tInit],  # Ma
    TPoints = Float64[model.TNow, model.TInit],    # Degrees C
    T₀ = Float64[model.TNow, model.TInit],
    ΔT = Float64[model.ΔTNow, model.ΔTInit],
)

# Default: No unconformity is imposed
unconf = (
    agePoints = Float64[],  # Ma
    TPoints = Float64[],    # Degrees C
)

# Add additional vectors for proposed unconformity and boundary points
unconf = (unconf...,
    agePointsₚ = similar(unconf.agePoints),
    TPointsₚ = similar(unconf.TPoints),
)
boundary = (boundary...,
    agePointsₚ = similar(boundary.agePoints),
    TPointsₚ = similar(boundary.TPoints),
)
model = (model...,
    ageSteps = Array{Float64}(model.tInit-model.dt/2 : -model.dt : 0+model.dt/2),
    tSteps = Array{Float64}(0+model.dt/2 : model.dt : model.tInit-model.dt/2),
)

## --- Invert for maximum likelihood t-T path

# This is where the "transdimensional" part comes in
agePoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
TPoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
# Fill some intermediate points to give the MCMC something to work with
Tr = 250 # Residence temperature
nPoints = 5
agePoints[1:nPoints] .= (model.tInit/30,model.tInit/4,model.tInit/2,model.tInit-model.tInit/4,model.tInit-model.tInit/30) # Ma
TPoints[1:nPoints] .= Tr  # Degrees C

# Run Markov Chain
@time "Compiling MCMC_vartcryst" MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
@time "Running MCMC_vartcryst" (TStepdist, HeAgedist, ndist, lldist, acceptancedist) = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)

@test isa(TStepdist, AbstractMatrix)
@test maximum(TStepdist) <= model.TInit
@test minimum(TStepdist) >= model.TNow

@test isa(HeAgedist, AbstractMatrix)
@test abs(sum(nanmean(HeAgedist, dims=2) - data.HeAge)/length(data.HeAge)) < 300

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxPoints

@test isa(lldist, AbstractVector)
@test -100 < mean(@view(lldist[model.burnin:end])) < 0

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.5, atol=0.35)

detail = (
    agemin = 0, # Youngest end of detail interval
    agemax = 541, # Oldest end of detail interval
    minpoints = 6, # Minimum number of points in detail interval
)
@time "Running MCMC_vartcryst with Detail interval" (TStepdist, HeAgedist, ndist, lldist, acceptancedist) = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary, detail)

@test isa(TStepdist, AbstractMatrix)
@test maximum(TStepdist) <= model.TInit
@test minimum(TStepdist) >= model.TNow

@test isa(HeAgedist, AbstractMatrix)
@test abs(sum(nanmean(HeAgedist, dims=2) - data.HeAge)/length(data.HeAge)) < 300

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxPoints

@test isa(lldist, AbstractVector)
@test -100 < mean(@view(lldist[model.burnin:end])) < 0

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.5, atol=0.35)
