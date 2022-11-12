# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                            ZrnHeInversion.jl                                  #
#                                                                               #
#       A transdimensional Bayesian Markov-chain Monte Carlo code for           #
#   inverting time-temperature paths from zircon U-Th/He ages, using a 1D       #
#   Crank-Nicholson finite difference diffusion model along with the zircon     #
#   damage and annealing model "ZRDAAM" from Guenthner et al, 2013 (AJS).       #
#   Requires uncorrected ZHe ages, grain size, [U] and [Th] in ppm.             #
#                                                                               #
#       This file contains the MCMC code, while "ZrnHe.jl" contains the damage  #
#   annealing function and the Crank-Nicholson code.                            #
#                                                                               #
#   © 2022 C. Brenhin Keller and Kalin McDannell                                #                                                                            #
#                                                                               #
#   If running for the first time, make sure you have instantiated the          #
#   manifest that came with this file                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## ---  Load required packages
    using Statistics, LinearAlgebra
    using StatsBase: fit, Histogram
    using ProgressMeter: @showprogress
    using StatGeochem
    using JLD: @load, @save
    using Plots; gr() # Use the GR backend for plotting

    using Thermochron

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)

    # # # # # # # # # # Choice of regional thermochron data # # # # # # # # # #

    # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    name = "MinnesotaInversion"
    data = importdataset("minnesota.csv", ',', importas=:Tuple)


## --- Prepare problem

    # How many steps of the Markov chain should we run?
    nsteps = 500000

    # How long should we wait for MC to converge (become stationary)
    burnin = 250000

    simplified = false # Prefer simpler tT paths?

    dt = 10 # time step size in Myr
    dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # Other model parameters
    CrystAgeMax = 4000.0 # Ma -- forbid anything older than this
    TStart = 400.0 # Temperature (in C)
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
    Uppm = data.U238_ppm
    Thppm = data.Th232_ppm
    HeAge = data.HeAge_Ma_raw
    HeAge_sigma = data.HeAge_Ma_sigma_raw
    CrystAge = data.CrystAge_Ma

    CrystAge[CrystAge .> CrystAgeMax] .= CrystAgeMax
    tStart = ceil(maximum(CrystAge)/dt) * dt

    tSteps = Array{Float64,1}(0+dt/2 : dt : tStart-dt/2)
    ageSteps = reverse(tSteps)
    ntSteps = length(tSteps) # Number of time steps
    eU = Uppm+.238*Thppm # Used only for plotting

## --- Test proscribed t-T paths with Neoproterozoic exhumation step

    # # Generate T path to test
    # Tr = 150
    # T0 = 30
    # agePoints = Float64[tStart, tStart*29/30,   720, 580, 250,  0] # Age (Ma)
    # TPoints  =  Float64[TStart,        Tr+T0, Tr+T0,  T0,  70, 10] # Temp. (C)
    # TSteps = linterp1s(agePoints,TPoints,ageSteps)
    #
    # # Plot t-T path
    # plot(ageSteps,TSteps,xflip=true)
    #
    # # Calculate model ages
    # CalcHeAges = Array{Float64}(undef, size(HeAge))
    # @time pr = DamageAnnealing(dt,tSteps,TSteps)
    # zircons = Array{Zircon{Float64}}(undef, length(Halfwidth))
    # @time for i=1:length(zircons)
    #     # Iterate through each grain, calculate the modeled age for each
    #     first_index = 1 + floor(Int64,(tStart - CrystAge[i])/dt)
    #     zircons[i] = Zircon(Halfwidth[i], dr, Uppm[i], Thppm[i], dt, ageSteps[first_index:end])
    #     CalcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionparams)
    # end
    #
    # # Plot Comparison of results
    # p2 = plot(eU, CalcHeAges, seriestype=:scatter,label="Model")
    # plot!(p2, eU, HeAge, yerror=HeAge_sigma*2, seriestype=:scatter, label="Data")
    # xlabel!(p2,"eU"); ylabel!(p2,"Age (Ma)")
    # display(p2)
    #
    # # Check log likelihood
    # ll = sum(-(CalcHeAges - HeAge).^2 ./ (2 .* HeAge_sigma.^2))
    # ll = sum(-(CalcHeAges - HeAge).^2 ./ (2 .* AnnealedSigma.^2))


## --- Invert for maximum likelihood t-T path

    # Boundary conditions (10C at present and 650 C at the time of zircon
    # formation). Optional: Apppend required unconformities at base of Cambrian,
    # or elsewhere
    boundary_agePoints = Float64[0, tStart] # Ma
    boundary_TPoints = Float64[0, TStart] # Degrees C

    # This is where the "transdimensional" part comes in
    nPoints = 0
    maxPoints = 40 # Maximum allowed number of age points
    agePoints = Array{Float64}(undef, maxPoints+1) # Array of fixed size to hold all optional age points
    TPoints = Array{Float64}(undef, maxPoints+1) # Array of fixed size to hold all optional age points

    # Fill some intermediate points to give the MCMC something to work with
    Tr = 250 # Residence temperature of initial proposal (value should not matter too much)
    for t in [tStart/30,tStart/4,tStart/2,tStart-tStart/4,tStart-tStart/30]
        global nPoints += 1
        agePoints[nPoints] = t # Ma
        TPoints[nPoints] = Tr # Degrees C
    end

    # # (Optional) Start with something close to the expected path
    #Tr = 150
    #T0 = 30
    #agePoints[1:4] = Float64[tStart*29/30,   720, 510, 250]) # Age (Ma)
    #TPoints[1:4]  =  Float64[       Tr+T0, Tr+T0,  T0,  70]) # Temp. (C)
    #nPoints+=4

    function MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)
        # Calculate model ages for initial proposal
        TSteps = linterp1s([view(agePoints, 1:nPoints) ; boundary_agePoints ; unconf_agePoints],
                           [view(TPoints, 1:nPoints) ; boundary_TPoints ; unconf_TPoints], ageSteps)
        CalcHeAges = Array{Float64}(undef, size(HeAge))
        pr = DamageAnnealing(dt,tSteps,TSteps) # Damage annealing history
        zircons = Array{Zircon{Float64}}(undef, length(Halfwidth))
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tStart - CrystAge[i])/dt)
            zircons[i] = Zircon(Halfwidth[i], dr, Uppm[i], Thppm[i], dt, ageSteps[first_index:end])
            CalcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionparams)
        end
        AnnealedSigma = simannealsigma.(1, HeAge_sigma; params=simannealparams)
        UnAnnealedSigma = simannealsigma.(nsteps, HeAge_sigma; params=simannealparams)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(HeAge, AnnealedSigma, CalcHeAges)
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
        HeAgeDist = Array{Float64}(undef, length(HeAge), nsteps)
        TStepDist = Array{Float64}(undef, ntSteps, nsteps)
        llDist = Array{Float64}(undef, nsteps)
        nDist = zeros(Int, nsteps)
        acceptanceDist = zeros(Bool, nsteps)

        # Standard deviations of Gaussian proposal distributions for temperature and time
        t_sigma = tStart/60
        T_sigma = TStart/60

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
                    elseif agePointsₚ[k] > (tStart - dt)
                        # Don't let any point get too close to tStart
                        agePointsₚ[k] -= (agePointsₚ[k] - (tStart - dt))
                    end
                    # Move the Temperature of one model point
                    if TPointsₚ[k] < 0
                        # Don't allow T<0
                        TPointsₚ[k] = 0
                    elseif TPointsₚ[k] > TStart
                        # Don't allow T>TStart
                        TPointsₚ[k] = TStart
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
                    agePointsₚ[nPointsₚ] = rand()*tStart
                    TPointsₚ[nPointsₚ] = rand()*TStart

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
                        # Allow the initial temperature to vary from TStart to TStart-50 C
                        boundary_TPointsₚ[2] = TStart-rand()*50
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
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(tStart - CrystAge[i])/dt)
                CalcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(TStepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionparams)
            end

            # Calculate log likelihood of proposal
            AnnealedSigma .= simannealsigma.(n, HeAge_sigma; params=simannealparams)
            llₚ = normpdf_ll(HeAge, AnnealedSigma, CalcHeAgesₚ)
            llₗ = normpdf_ll(HeAge, AnnealedSigma, CalcHeAges) # Recalulate last one too with new AnnealedSigma
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
            llDist[n] = normpdf_ll(HeAge, UnAnnealedSigma, CalcHeAges) # Recalculated to constant baseline
            nDist[n] = nPoints # Distribution of # of points
            HeAgeDist[:,n] = CalcHeAges # Distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            TStepDist[:,n] = TSteps # Distribution of T paths

        end
        return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
    end

    # Run Markov Chain
    @time (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)

    # # Save results using JLD
    @save string(name, ".jld") ageSteps tSteps TStepDist burnin nsteps TStart tStart

    # Plot sample age-eU correlations
    h = scatter(eU,HeAgeDist[:,burnin:burnin+50], label="")
    plot!(h, eU, HeAge, yerror=HeAge_sigma, seriestype=:scatter, label="Data")
    xlabel!(h,"eU (ppm)"); ylabel!(h,"Age (Ma)")
    savefig(h,string(name,"Age-eU.pdf"))

    #


## --- Create image of paths

    # Resize the post-burnin part of the stationary distribution
    TStepDistResized = Array{Float64}(undef, 2001, size(TStepDist,2)-burnin)
    xq = collect(range(0,tStart,length=2001))
    for i=1:size(TStepDist,2)-burnin
        TStepDistResized[:,i] = linterp1s(tSteps,TStepDist[:,i+burnin],xq)
    end

    # Calculate composite image
    tTimage = zeros(ceil(Int, TStart)*2, size(TStepDistResized,1))
    yq = collect(0:0.5:TStart)
    for i=1:size(TStepDistResized,1)
        hist = fit(Histogram,TStepDistResized[i,:],yq,closed=:right)
        tTimage[:,i] = hist.weights
    end

## --- Plot image with 'ylcn' custom colorscale
    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    A = imsc(reverse(tTimage,dims=2), ylcn, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)",ylabel="Temperature (°C)",yticks=0:50:400,xticks=0:500:3500,yminorticks=5,xminorticks=5,tick_dir=:out,framestyle=:box)
    plot!(k[1],xq,cntr(yq),A,yflip=true,xflip=true,legend=false,aspectratio=3500/400/1.5,xlims=(0,3500),ylims=(0,400))

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, xticks=false, framestyle=:box, yflip=false)

    #plot!([659, 717.4, 717.4, 659, 659],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Sturtian glacial
    plot!([635.5, 717.4, 717.4, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.5) #Sturtian & Marinoan glacial
    #plot!([635.5, 650.0, 650.0, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Marinoan glacial
    #plot!([480, 640, 640, 480, 480],[0, 0, 50, 50, 0], linestyle = :dot, color=:black, linewidth=1.25) # t-T box 640 to 480 Ma, 0-50°C


    savefig(k,"minnesota-mcmc-zrdaam-unc.pdf")
    display(k)

    #img = plot(reverse(xq), cntr(yq), imsc(reverse(tTimage,dims=2), viridis, 0, nanpctile(tTimage[:],97.5)), yflip=true, xflip=true, legend=false,
    #     xlabel="Age (Ma)", ylabel="Temperature (°C)")
    #savefig(img,string(name,"minnesota-mcmc.pdf"))
    #display(img)

## ---
