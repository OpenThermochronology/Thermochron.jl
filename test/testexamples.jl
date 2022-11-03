using Statistics, LinearAlgebra
using StatsBase: fit, Histogram
using ProgressMeter: @showprogress
using StatGeochem

# Read in data from file using StatGeochem
datapath = joinpath("..", "examples", "minnesota.csv")
data = importdataset(datapath, ',', importas=:Tuple)

# Default: No unconformity is imposed
unconf_agePoints = Array{Float64}([])
unconf_TPoints = Array{Float64}([])

# How many steps of the Markov chain should we run?
nsteps = 50

# How long should we wait for MC to converge (become stationary)
burnin = 25


# Model uncertainty is not well known (depends on annealing parameters,
# decay constants, diffusion parameters, etc.), but is certainly non-zero.
# Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
simannealparams = (
    25.0, # Model uncertainty [Ma]
    35.0, # Initial uncertainty [Ma]
    10 ./ burnin, # lambda [1/n]
)

simplified = false
CrystAgeMax_Ma = 4000.0 # Ma -- forbid anything older than this

# Other model parameters
TCryst = 400.0 # Temperature (in C)
dr = 1 # Radius step, in microns

diffusionparams = (;
    DzEa = 165.0, # kJ/mol
    DzD0 = 193188.0, # cm^2/sec
    DN17Ea = 71.0, # kJ/mol
    DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
)

# # # # # # # # # # Choice of regional thermochron data # # # # # # # # # # #

# # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
name = "MinnesotaInversion"
dt = 10 # time step size in Myr
dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Uncomment this section if you wish to impose an unconformity at any point in the record
# unconf_agePoints = [600.0] # Ma
# unconf_TPoints = [20.0] # Degrees C
# unconf_params = [560,80,0,50] # Age0, DAge, T0, DT
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

AnnealedSigma = simannealsigma.(1, HeAge_Ma_sigma; params=simannealparams)

## --- Invert for maximum likelihood t-T path

# Boundary conditions (10C at present and 650 C at the time of zircon
# formation). Optional: Apppend required unconformities at base of Cambrian,
# or elsewhere
boundary_agePoints = Array{Float64}([0; tCryst;]) # Ma
boundary_TPoints = Array{Float64}([0; TCryst;]) # Degrees C

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
#agePoints[1:4] = Array{Float64}([tCryst*29/30,   720, 510, 250]) # Age (Ma)
#TPoints[1:4]  =  Array{Float64}([       Tr+T0, Tr+T0,  T0,  70]) # Temp. (C)
#nPoints+=4

function MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)
    # Calculate model ages for initial proposal
    TSteps = linterp1s([agePoints[1:nPoints] ; boundary_agePoints ; unconf_agePoints],
                        [TPoints[1:nPoints] ; boundary_TPoints ; unconf_TPoints], ageSteps)
    CalcHeAges = Array{Float64}(undef, size(HeAge_Ma))
    pr = DamageAnnealing(dt,tSteps,TSteps) # Damage annealing history
    for i=1:length(Halfwidth)
        # Iterate through each grain, calculate the modeled age for each
        first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
        CalcHeAges[i] = ZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i], diffusionparams)
    end
    CalcHeAges_prop = copy(CalcHeAges)

    # Log-likelihood for initial proposal
    if simplified
        ll = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2)) - log(nPoints)
    else
        ll = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2))
    end
    ll_prop = copy(ll) # Initialize

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

    @showprogress 10 "Running MCMC..." for n=1:nsteps

        # Copy proposal from last accepted solution
        nPoints_prop = copy(nPoints)
        agePoints_prop = copy(agePoints)
        TPoints_prop = copy(TPoints)
        TSteps_prop = copy(TSteps)
        unconf_agePoints_prop = copy(unconf_agePoints)
        unconf_TPoints_prop = copy(unconf_TPoints)
        boundary_TPoints_prop = copy(boundary_TPoints)

        # Adjust the proposal
        r = rand()
        if r < move
            # Move the age of one model point
            for i=1:10000 # Try 10000 times to satisfy the reheating rate limit
                k = Int(ceil(rand() * nPoints))

                agePoints_prop[k] += randn() * t_sigma
                if agePoints_prop[k] < dt
                    # Don't let any point get too close to 0
                    agePoints_prop[k] += (dt - agePoints_prop[k])
                elseif agePoints_prop[k] > (tCryst - dt)
                    # Don't let any point get too close to tCryst
                    agePoints_prop[k] -= (agePoints_prop[k] - (tCryst - dt))
                end
                # Move the Temperature of one model point
                if TPoints_prop[k] < 0
                    # Don't allow T<0
                    TPoints_prop[k] = 0
                elseif TPoints_prop[k] > TCryst
                    # Don't allow T>TCryst
                    TPoints_prop[k] = TCryst
                end

                # Interpolate proposed t-T path
                TSteps_prop = linterp1s([agePoints_prop[1:nPoints_prop] ; boundary_agePoints ; unconf_agePoints_prop],
                                          [TPoints_prop[1:nPoints_prop] ; boundary_TPoints_prop ; unconf_TPoints_prop], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TSteps_prop)) < dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                agePoints_prop = copy(agePoints)
                TPoints_prop = copy(TPoints)
            end

        elseif (r < move+birth) && (nPoints_prop < maxPoints)
            # Birth: add a new model point
            nPoints_prop += 1
            for i=1:10000 # Try 10000 times to satisfy the reheating rate limit
                agePoints_prop[nPoints_prop] = rand()*tCryst
                TPoints_prop[nPoints_prop] = rand()*TCryst

                # Interpolate proposed t-T path
                TSteps_prop = linterp1s([agePoints_prop[1:nPoints_prop] ; boundary_agePoints ; unconf_agePoints_prop],
                                        [TPoints_prop[1:nPoints_prop] ; boundary_TPoints_prop ; unconf_TPoints_prop], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TSteps_prop)) < dTmax && break
            end

        elseif (r < move+birth+death) && (r >= move+birth) && (nPoints_prop > 1)
            # Death: remove a model point
            nPoints_prop -= 1 # Delete last point in array from proposal
            for i=1:10000 # Try 10000 times to satisfy the reheating rate limit
                k = Int(ceil(rand()*nPoints)) # Choose point to delete
                agePoints_prop[k] = agePoints_prop[nPoints]
                TPoints_prop[k] = TPoints_prop[nPoints]

                # Interpolate proposed t-T path
                TSteps_prop = linterp1s([agePoints_prop[1:nPoints_prop] ; boundary_agePoints ; unconf_agePoints_prop],
                                        [TPoints_prop[1:nPoints_prop] ; boundary_TPoints_prop ; unconf_TPoints_prop], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TSteps_prop)) < dTmax && break
            end

        else
            # Move boundary conditions
            for i=1:10000
                r2 = rand()
                if r2 < 0.5
                    # Allow the present temperature to vary from 0 to 10 degrees C
                    boundary_TPoints_prop[1] = 0+rand()*10
                else
                    # Allow the initial temperature to vary from TCryst to TCryst-50 C
                    boundary_TPoints_prop[2] = TCryst-rand()*50
                end
                if length(unconf_agePoints_prop) > 0
                    # If there's an imposed unconformity, adjust within parameters
                    unconf_agePoints_prop = unconf_params[1] + rand()*unconf_params[2]
                    unconf_TPoints_prop = unconf_params[3] + rand()*unconf_params[4]
                end

                # Recalculate interpolated proposed t-T path
                TSteps_prop = linterp1s([agePoints_prop[1:nPoints_prop] ; boundary_agePoints ; unconf_agePoints_prop],
                                        [TPoints_prop[1:nPoints_prop] ; boundary_TPoints_prop ; unconf_TPoints_prop], ageSteps)

                # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                maximum(diff(TSteps_prop)) < dTmax && break

                # Copy last accepted solution to re-modify if we don't break
                unconf_agePoints_prop = copy(unconf_agePoints)
                unconf_TPoints_prop = copy(unconf_TPoints)
                boundary_TPoints_prop = copy(boundary_TPoints)
            end
        end

        # Recalculate interpolated proposed t-T path
        TSteps_prop = linterp1s([agePoints_prop[1:nPoints_prop] ; boundary_agePoints ; unconf_agePoints_prop],
                                [TPoints_prop[1:nPoints_prop] ; boundary_TPoints_prop ; unconf_TPoints_prop], ageSteps)

         # Calculate model ages for each grain
        pr = DamageAnnealing(dt,tSteps,TSteps_prop)
        for i=1:length(Halfwidth)
            first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
            CalcHeAges_prop[i] = ZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps_prop[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i],diffusionparams)
        end

        # Calculate log likelihood of proposal
        AnnealedSigma .= simannealsigma.(n, HeAge_Ma_sigma; params=simannealparams)

        # Recalculate ll in case annealed sigma has changed
        if simplified
            # if simplified, slightly penalize more complex t-T paths
            ll_last = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2)) - log(nPoints)
            ll_prop = sum(-(CalcHeAges_prop - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2)) - log(nPoints_prop)
        else
            ll_last = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2))
            ll_prop = sum(-(CalcHeAges_prop - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2))
        end

        # Accept or reject proposal based on likelihood
        # To avoid numerical problems with diffusion code, also reject proposal
        # if maximum proposed heating rate is greater than 25C per timestep.
        # (Fast cooling should not be a problem, however)
        if rand()<exp(ll_prop-ll_last)
            ll = copy(ll_prop)
            nPoints = copy(nPoints_prop)
            agePoints = copy(agePoints_prop)
            TPoints = copy(TPoints_prop)
            unconf_agePoints = copy(unconf_agePoints_prop)
            unconf_TPoints = copy(unconf_TPoints_prop)
            boundary_TPoints = copy(boundary_TPoints_prop)
            CalcHeAges = copy(CalcHeAges_prop)

            # These are saved for ouput, but not critical to the function of the MCMC loop
            TSteps = copy(TSteps_prop)
            acceptanceDist[n] = 1
        end

        # Record results for analysis and troubleshooting
        llDist[n] = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* simannealsigma.(nsteps, HeAge_Ma_sigma; params=simannealparams).^2)) # Recalculated to constant baseline
        nDist[n] = nPoints # Distribution of # of points
        HeAgeDist[:,n] = CalcHeAges # Distribution of He ages

        # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
        TStepDist[:,n] = TSteps # Distribution of T paths

    end
    return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
end

# Run Markov Chain
(TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)

@test isa(TStepDist, AbstractMatrix)
@test isa(HeAgeDist, AbstractMatrix)
@test isa(nDist, AbstractVector{Int})
@test isa(llDist, AbstractVector)
@test isa(acceptanceDist, AbstractVector{Bool})
