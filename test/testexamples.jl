using Statistics, LinearAlgebra
using ProgressMeter: @showprogress
using StatGeochem

# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesota.csv")
ds = importdataset(datapath, ',', importas=:Tuple);

## --- Prepare problem

annealingburnin = 100
model = (
    nsteps = 200, # How many steps of the Markov chain should we run?
    burnin = 100, # How long should we wait for MC to converge (become stationary)
    dr = 1.0,    # Radius step, in microns
    dt = 10.0,   # time step size in Myr
    dTmax = 10.0, # Maximum reheating/burial per model timestep (to prevent numerical overflow)
    TInit = 400.0, # Initial model temperature (in C) (i.e., crystallization temperature)
    ΔTInit = -50.0, # TInit can vary from TInit to TInit+ΔTinit
    TNow = 0.0, # Current surface temperature (in C)
    ΔTNow = 10.0, # TNow may vary from TNow to TNow+ΔTNow
    tInitMax = 4000.0, # Ma -- forbid anything older than this
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
    λAnnealing = 10 ./ annealingburnin # Annealing decay [1/n]
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

tSteps = Array{Float64,1}(0+model.dt/2 : model.dt : model.tInit-model.dt/2)
ageSteps = reverse(tSteps)

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

# # Uncomment this section if you wish to impose an unconformity at any point in the record
# # Uniform distributions from Age₀ to Age₀+ΔAge, T₀ to T₀+ΔT,
# unconf = (;
#     agePoints = Float64[560.0,],  # Ma
#     TPoints = Float64[20.0,],     # Degrees C
#     Age₀ = Float64[500,],
#     ΔAge = Float64[80,],
#     T₀ = Float64[0,],
#     ΔT = Float64[50,],
# )
## --- Invert for maximum likelihood t-T path

# Add additional vectors for proposed unconformity and boundary points
unconf = (unconf...,
    agePointsₚ = similar(unconf.agePoints),
    TPointsₚ = similar(unconf.TPoints),
)
boundary = (boundary...,
    agePointsₚ = similar(boundary.agePoints),
    TPointsₚ = similar(boundary.TPoints),
)

# This is where the "transdimensional" part comes in
agePoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
TPoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
# Fill some intermediate points to give the MCMC something to work with
Tr = 250 # Residence temperature
nPoints = 5
agePoints[1:nPoints] .= [model.tInit/30,model.tInit/4,model.tInit/2,model.tInit-model.tInit/4,model.tInit-model.tInit/30] # Ma
TPoints[1:nPoints] .+ Tr  # Degrees C

function MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
    # Calculate model ages for initial proposal
    TSteps = linterp1s([view(agePoints, 1:nPoints) ; boundary.agePoints ; unconf.agePoints],
                       [view(TPoints, 1:nPoints) ; boundary.TPoints ; unconf.TPoints], ageSteps)
    calcHeAges = Array{Float64}(undef, size(data.HeAge))
    pr = anneal(model.dt, tSteps, TSteps, :zrdaam) # Damage annealing history

    # Prepare a Mineral object for each analysis
    zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
    for i=1:length(zircons)
        # Iterate through each grain, calculate the modeled age for each
        first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
        zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, ageSteps[first_index:end])
        calcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), model)
    end

    # Simulated annealing of uncertainty
    σₐ = simannealsigma.(1, data.HeAge_sigma; params=model)
    σₙ = simannealsigma.(model.nsteps, data.HeAge_sigma; params=model)

    # Log-likelihood for initial proposal
    ll = normpdf_ll(data.HeAge, σₐ, calcHeAges)
    if model.simplified
        ll -= log(nPoints)
    end

    # Variables to hold proposals
    llₚ = ll
    nPointsₚ = nPoints
    agePointsₚ = similar(agePoints)
    TPointsₚ = similar(TPoints)
    TStepsₚ = similar(TSteps)
    calcHeAgesₚ = similar(calcHeAges)

    # distributions to populate
    HeAgedist = Array{Float64}(undef, length(data.HeAge), model.nsteps)
    TStepdist = Array{Float64}(undef, length(tSteps), model.nsteps)
    lldist = Array{Float64}(undef, model.nsteps)
    ndist = zeros(Int, model.nsteps)
    acceptancedist = zeros(Bool, model.nsteps)

    # Standard deviations of Gaussian proposal distributions for temperature and time
    t_sigma = model.tInit/60
    T_sigma = model.TInit/60

    # Proposal probabilities (must sum to 1)
    move = 0.65
    birth = 0.15
    death = 0.15 # Should equal birth
    boundarymove = 0.05
    maxattempts = 1000

    @showprogress 10 "Running MCMC..." for n=1:model.nsteps

        # Copy proposal from last accepted solution
        nPointsₚ = nPoints
        copyto!(agePointsₚ, agePoints)
        copyto!(TPointsₚ, TPoints)
        copyto!(TStepsₚ, TSteps)
        copyto!(unconf.agePointsₚ, unconf.agePoints)
        copyto!(unconf.TPointsₚ, unconf.TPoints)
        copyto!(boundary.TPointsₚ, boundary.TPoints)

        # Adjust the proposal
        r = rand()
        if r < move
            # Move the age of one model point
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                k = ceil(Int, rand() * nPoints)

                agePointsₚ[k] += randn() * t_sigma
                if agePointsₚ[k] < model.dt
                    # Don't let any point get too close to 0
                    agePointsₚ[k] += (model.dt - agePointsₚ[k])
                elseif agePointsₚ[k] > (model.tInit - model.dt)
                    # Don't let any point get too close to tInit
                    agePointsₚ[k] -= (agePointsₚ[k] - (model.tInit - model.dt))
                end
                # Move the Temperature of one model point
                if TPointsₚ[k] < 0
                    # Don't allow T<0
                    TPointsₚ[k] = 0
                elseif TPointsₚ[k] > model.TInit
                    # Don't allow T>TInit
                    TPointsₚ[k] = model.TInit
                end

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < model.dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                copyto!(agePointsₚ, agePoints)
                copyto!(TPointsₚ, TPoints)
            end
        elseif (r < move+birth) && (nPointsₚ < model.maxPoints)
            # Birth: add a new model point
            nPointsₚ += 1
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                agePointsₚ[nPointsₚ] = rand()*model.tInit
                TPointsₚ[nPointsₚ] = rand()*model.TInit

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < model.dTmax && break
            end
        elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > 1)
            # Death: remove a model point
            nPointsₚ -= 1 # Delete last point in array from proposal
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                k = ceil(Int, rand()*nPoints) # Choose point to delete
                agePointsₚ[k] = agePointsₚ[nPoints]
                TPointsₚ[k] = TPointsₚ[nPoints]

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < model.dTmax && break
            end
        else
            # Move boundary conditions
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                # Move the temperatures of the starting and ending boundaries
                @. boundary.TPointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if length(unconf.agePoints) > 0
                    @. unconf.agePointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.TPointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end


                # Recalculate interpolated proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < model.dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                copyto!(unconf.agePointsₚ, unconf.agePoints)
                copyto!(unconf.TPointsₚ, unconf.TPoints)
                copyto!(boundary.TPointsₚ, boundary.TPoints)
            end
        end

        # Recalculate interpolated proposed t-T path
        TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                            [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

         # Calculate model ages for each grain
        anneal!(pr, model.dt, tSteps, TStepsₚ, :zrdaam)
        for i=1:length(zircons)
            first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
            calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(TStepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), model)
        end

        # Calculate log likelihood of proposal
        σₐ .= simannealsigma.(n, data.HeAge_sigma; params=model)
        llₚ = normpdf_ll(data.HeAge, σₐ, calcHeAgesₚ)
        llₗ = normpdf_ll(data.HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
        if model.simplified # slightly penalize more complex t-T paths
            llₚ -= log(nPointsₚ)
            llₗ -= log(nPoints)
        end

        # Accept or reject proposal based on likelihood
        # To avoid numerical problems with diffusion code, also reject proposal
        # if maximum proposed heating rate is greater than 25C per timestep.
        # (Fast cooling should not be a problem, however)
        if log(rand()) < (llₚ - llₗ)
            ll = llₚ
            nPoints = nPointsₚ
            copyto!(agePoints, agePointsₚ)
            copyto!(TPoints, TPointsₚ)
            copyto!(unconf.agePoints, unconf.agePointsₚ)
            copyto!(unconf.TPoints, unconf.TPointsₚ)
            copyto!(boundary.TPoints, boundary.TPointsₚ)
            copyto!(calcHeAges, calcHeAgesₚ)

            # These are saved for ouput, but not critical to the function of the MCMC loop
            copyto!(TSteps, TStepsₚ)
            acceptancedist[n] = true
        end

        # Record results for analysis and troubleshooting
        lldist[n] = normpdf_ll(data.HeAge, σₙ, calcHeAges) # Recalculated to constant baseline
        ndist[n] = nPoints # distribution of # of points
        HeAgedist[:,n] = calcHeAges # distribution of He ages

        # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
        TStepdist[:,n] = TSteps # distribution of T paths

    end
    return (TStepdist, HeAgedist, ndist, lldist, acceptancedist)
end


# Run Markov Chain
@time MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
@time (TStepdist, HeAgedist, ndist, lldist, acceptancedist) = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)

@test isa(TStepdist, AbstractMatrix)
@test maximum(TStepdist) <= model.TInit
@test minimum(TStepdist) >= model.TNow

@test isa(HeAgedist, AbstractMatrix)
@test abs(sum(nanmean(HeAgedist, dims=2) - data.HeAge)/length(data.HeAge)) < 100

@test isa(ndist, AbstractVector{Int})
@test minimum(ndist) >= 0
@test maximum(ndist) <= model.maxPoints

@test isa(lldist, AbstractVector)
@test isapprox(mean(@view(lldist[model.burnin:end])), -50, atol=10)

@test isa(acceptancedist, AbstractVector{Bool})
@test isapprox(mean(acceptancedist), 0.58, atol=0.3)
