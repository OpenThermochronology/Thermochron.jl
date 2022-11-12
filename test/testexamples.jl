using Statistics, LinearAlgebra
using StatsBase: fit, Histogram
using ProgressMeter: @showprogress
using StatGeochem

# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesota.csv")
data = importdataset(datapath, ',', importas=:Tuple)

# How many steps of the Markov chain should we run?
nsteps = 200

# How long should we wait for MC to converge (become stationary)
burnin = 100

## --- Prepare problem

simplified = false # Prefer simpler tT paths?

dt = 10 # time step size in Myr
dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

# Other model parameters
CrystAgeMax_Ma = 4000.0 # Ma -- forbid anything older than this
TCryst = 400.0 # Temperature (in C)
dr = 1 # Radius step, in microns

diffusionparams = (;
    DzEa = 165.0, # kJ/mol
    DzD0 = 193188.0, # cm^2/sec
    DN17Ea = 71.0, # kJ/mol
    DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
)

# Model uncertainty is not well known (depends on annealing parameters,
# decay constants, diffusion parameters, etc.), but is certainly non-zero.
# Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
simannealparams = (
    25.0, # Model uncertainty [Ma]
    35.0, # Initial uncertainty [Ma]
    10 ./ burnin, # lambda [1/n]
)

# Default: No unconformity is imposed
unconf_agePoints = Float64[]
unconf_TPoints = Float64[]

# # Uncomment this section if you wish to impose an unconformity at any point in the record
# unconf_agePoints = Float64[600.0] # Ma
# unconf_TPoints = Float64[20.0] # Degrees C
# unconf_params = Float64[560,80,0,50] # Age0, DAge, T0, DT
# name = "$(name)Imposed-lowff"

# Populate local variables from data frame with specified options
Halfwidth = data.Halfwidth_um
U_ppm = data.U238_ppm
Th_ppm = data.Th232_ppm
HeAge_Ma = data.HeAge_Ma_raw
HeAge_Ma_sigma = data.HeAge_Ma_sigma_raw
CrystAge_Ma = data.CrystAge_Ma

CrystAge_Ma[CrystAge_Ma .> CrystAgeMax_Ma] .= CrystAgeMax_Ma
tCryst = ceil(maximum(CrystAge_Ma)/dt) * dt

tSteps = Array{Float64,1}(0+dt/2 : dt : tCryst-dt/2)
ageSteps = reverse(tSteps)
ntSteps = length(tSteps) # Number of time steps
eU = U_ppm+.238*Th_ppm # Used only for plotting


## --- Invert for maximum likelihood t-T path

# Boundary conditions (10C at present and 650 C at the time of zircon
# formation). Optional: Apppend required unconformities at base of Cambrian,
# or elsewhere
boundary_agePoints = Float64[0, tCryst] # Ma
boundary_TPoints = Float64[0, TCryst] # Degrees C

# This is where the "transdimensional" part comes in
nPoints = 0
maxPoints = 40 # Maximum allowed number of age points
agePoints = Array{Float64}(undef, maxPoints+1) # Array of fixed size to hold all optional age points
TPoints = Array{Float64}(undef, maxPoints+1) # Array of fixed size to hold all optional age points

# Fill some intermediate points to give the MCMC something to work with
Tr = 250 # Residence temperature of initial proposal (value should not matter too much)
for t in [tCryst/30,tCryst/4,tCryst/2,tCryst-tCryst/4,tCryst-tCryst/30]
    global nPoints += 1
    agePoints[nPoints] = t # Ma
    TPoints[nPoints] = Tr # Degrees C
end

# # (Optional) Start with something close to the expected path
#Tr = 150
#T0 = 30
#agePoints[1:4] = Float64[tCryst*29/30,   720, 510, 250]) # Age (Ma)
#TPoints[1:4]  =  Float64[       Tr+T0, Tr+T0,  T0,  70]) # Temp. (C)
#nPoints+=4

function MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)
    # Calculate model ages for initial proposal
    TSteps = linterp1s([view(agePoints, 1:nPoints) ; boundary_agePoints ; unconf_agePoints],
                       [view(TPoints, 1:nPoints) ; boundary_TPoints ; unconf_TPoints], ageSteps)
    CalcHeAges = Array{Float64}(undef, size(HeAge_Ma))
    pr = DamageAnnealing(dt,tSteps,TSteps) # Damage annealing history
    for i=1:length(Halfwidth)
        # Iterate through each grain, calculate the modeled age for each
        first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
        CalcHeAges[i] = ZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i], diffusionparams)
    end
    AnnealedSigma = simannealsigma.(1, HeAge_Ma_sigma; params=simannealparams)
    UnAnnealedSigma = simannealsigma.(nsteps, HeAge_Ma_sigma; params=simannealparams)

    # Log-likelihood for initial proposal
    ll = normpdf_ll(HeAge_Ma, AnnealedSigma, CalcHeAges)
    if simplified
        ll -= log(nPoints)
    end

    # Variables to hold proposals
    llₚ = ll
    nPointsₚ = nPoints
    agePointsₚ = similar(agePoints)
    TPointsₚ = similar(TPoints)
    TStepsₚ = similar(TSteps)
    unconf_agePointsₚ = similar(unconf_agePoints)
    unconf_TPointsₚ = similar(unconf_TPoints)
    boundary_TPointsₚ = similar(boundary_TPoints)
    CalcHeAgesₚ = similar(CalcHeAges)

    # Distributions to populate
    HeAgeDist = Array{Float64}(undef, length(HeAge_Ma), nsteps)
    TStepDist = Array{Float64}(undef, ntSteps, nsteps)
    llDist = Array{Float64}(undef, nsteps)
    nDist = zeros(Int, nsteps)
    acceptanceDist = zeros(Bool, nsteps)

    # Standard deviations of Gaussian proposal distributions for temperature and time
    t_sigma = tCryst/60
    T_sigma = TCryst/60

    # Proposal probabilities (must sum to 1)
    move = 0.64
    birth = 0.15
    death = 0.15 # Should equal birth
    boundary = 0.06
    maxattempts = 1000

    @showprogress 10 "Running MCMC..." for n=1:nsteps

        # Copy proposal from last accepted solution
        nPointsₚ = nPoints
        copyto!(agePointsₚ, agePoints)
        copyto!(TPointsₚ, TPoints)
        copyto!(TStepsₚ, TSteps)
        copyto!(unconf_agePointsₚ, unconf_agePoints)
        copyto!(unconf_TPointsₚ, unconf_TPoints)
        copyto!(boundary_TPointsₚ, boundary_TPoints)

        # Adjust the proposal
        r = rand()
        if r < move
            # Move the age of one model point
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                k = ceil(Int, rand() * nPoints)

                agePointsₚ[k] += randn() * t_sigma
                if agePointsₚ[k] < dt
                    # Don't let any point get too close to 0
                    agePointsₚ[k] += (dt - agePointsₚ[k])
                elseif agePointsₚ[k] > (tCryst - dt)
                    # Don't let any point get too close to tCryst
                    agePointsₚ[k] -= (agePointsₚ[k] - (tCryst - dt))
                end
                # Move the Temperature of one model point
                if TPointsₚ[k] < 0
                    # Don't allow T<0
                    TPointsₚ[k] = 0
                elseif TPointsₚ[k] > TCryst
                    # Don't allow T>TCryst
                    TPointsₚ[k] = TCryst
                end

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary_agePoints ; unconf_agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary_TPointsₚ ; unconf_TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                copyto!(agePointsₚ, agePoints)
                copyto!(TPointsₚ, TPoints)
            end
        elseif (r < move+birth) && (nPointsₚ < maxPoints)
            # Birth: add a new model point
            nPointsₚ += 1
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                agePointsₚ[nPointsₚ] = rand()*tCryst
                TPointsₚ[nPointsₚ] = rand()*TCryst

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary_agePoints ; unconf_agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary_TPointsₚ ; unconf_TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < dTmax && break
            end
        elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > 1)
            # Death: remove a model point
            nPointsₚ -= 1 # Delete last point in array from proposal
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                k = ceil(Int, rand()*nPoints) # Choose point to delete
                agePointsₚ[k] = agePointsₚ[nPoints]
                TPointsₚ[k] = TPointsₚ[nPoints]

                # Interpolate proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary_agePoints ; unconf_agePointsₚ],
                                    [view(TPointsₚ, 1:nPointsₚ) ; boundary_TPointsₚ ; unconf_TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < dTmax && break
            end
        else
            # Move boundary conditions
            for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                r2 = rand()
                if r2 < 0.5
                    # Allow the present temperature to vary from 0 to 10 degrees C
                    boundary_TPointsₚ[1] = 0+rand()*10
                else
                    # Allow the initial temperature to vary from TCryst to TCryst-50 C
                    boundary_TPointsₚ[2] = TCryst-rand()*50
                end
                if length(unconf_agePointsₚ) > 0
                    # If there's an imposed unconformity, adjust within parameters
                    unconf_agePointsₚ = unconf_params[1] + rand()*unconf_params[2]
                    unconf_TPointsₚ = unconf_params[3] + rand()*unconf_params[4]
                end

                # Recalculate interpolated proposed t-T path
                TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary_agePoints ; unconf_agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary_TPointsₚ ; unconf_TPointsₚ], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TStepsₚ)) < dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                copyto!(unconf_agePointsₚ, unconf_agePoints)
                copyto!(unconf_TPointsₚ, unconf_TPoints)
                copyto!(boundary_TPointsₚ, boundary_TPoints)
            end
        end

        # Recalculate interpolated proposed t-T path
        TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary_agePoints ; unconf_agePointsₚ],
                            [view(TPointsₚ, 1:nPointsₚ) ; boundary_TPointsₚ ; unconf_TPointsₚ], ageSteps)

         # Calculate model ages for each grain
        DamageAnnealing!(pr, dt, tSteps, TStepsₚ)
        for i=1:length(Halfwidth)
            first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
            @views CalcHeAgesₚ[i] = ZrnHeAgeSpherical(dt,ageSteps[first_index:end],TStepsₚ[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i],diffusionparams)
        end

        # Calculate log likelihood of proposal
        AnnealedSigma .= simannealsigma.(n, HeAge_Ma_sigma; params=simannealparams)
        llₚ = normpdf_ll(HeAge_Ma, AnnealedSigma, CalcHeAgesₚ)
        llₗ = normpdf_ll(HeAge_Ma, AnnealedSigma, CalcHeAges) # Recalulate last one too with new AnnealedSigma
        if simplified # slightly penalize more complex t-T paths
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
            copyto!(unconf_agePoints, unconf_agePointsₚ)
            copyto!(unconf_TPoints, unconf_TPointsₚ)
            copyto!(boundary_TPoints, boundary_TPointsₚ)
            copyto!(CalcHeAges, CalcHeAgesₚ)

            # These are saved for ouput, but not critical to the function of the MCMC loop
            copyto!(TSteps, TStepsₚ)
            acceptanceDist[n] = true
        end

        # Record results for analysis and troubleshooting
        llDist[n] = normpdf_ll(HeAge_Ma, UnAnnealedSigma, CalcHeAges) # Recalculated to constant baseline
        nDist[n] = nPoints # Distribution of # of points
        HeAgeDist[:,n] = CalcHeAges # Distribution of He ages

        # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
        TStepDist[:,n] = TSteps # Distribution of T paths

    end
    return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
end

# Run Markov Chain
@time (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)
@time (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)


@test isa(TStepDist, AbstractMatrix)
@test isa(HeAgeDist, AbstractMatrix)
@test isa(nDist, AbstractVector{Int})
@test isa(llDist, AbstractVector)
@test isa(acceptanceDist, AbstractVector{Bool})
