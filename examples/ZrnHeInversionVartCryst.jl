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
#   Last modified by C. Brenhin Keller 2018-03-13                               #
#   Last modified by K. McDannell 2022-10-19                                    #                                                                            #
#                                                                               #
#   If running for the first time, make sure you have instantiated the          #
#   manifest that came with this file                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## ---  Load required packages
    using Statistics, LinearAlgebra, SparseArrays
    using StatsBase: fit, Histogram
    using StatGeochem
    using JLD: @load, @save
    using Plots; gr() # Use the GR backend for plotting

    using Thermchron

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)


## --- Prepare problem

    # Read in data from file using StatGeochem
    data = importdataset("minnesota.csv", ',', importas=:Tuple)

    # Default: No unconformity is imposed
    unconf_agePoints = Array{Float64}([])
    unconf_TPoints = Array{Float64}([])

    # How many steps of the Markov chain should we run?
    nsteps = 500000

    # How long should we wait for MC to converge (become stationary)
    burnin = 250000

    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
    simannealparams = (;
        ModelUncertainty = 25.0, # Ma
        InitialUncertainty = 35.0, # [Ma]
        lambda = 10 ./ burnin, # [1/n]
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
     #test = (data[:HeAge_reference] .== "Guenthner et al. 2013") .| (data[:HeAge_reference] .== "Miltich 2005")
     dt = 10 # time step size in Myr
     dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)
    # # simplified = true

    # # New samples from the East Central Minnesota Batholith, collected by NLS-H
    # name = "ECMBInversionSphericalVartCryst"
    # test = (data[:Geological_Province] .== "East Central Minnesota Batholith") # .& (data[:HeAge_reference] .== "This study")
    # # Exclude two grains that are suspiciously old for their size
    # # test = test .& .~(data[:Grain_Name] .== "ECMB2-3") .& .~(data[:Grain_Name] .== "ECMB6-1")
    # dt = 10 # time step size in Myr
    # dTmax = 5.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # # New samples from the Keewenawan / Midcontinent Rift
    # name = "KeweenawanInversionSphericalVartCryst"
    # test = data[:Geological_Province] .== "Keweenawan Rift"
    # dt = 10 # time step size in Myr
    # dTmax = 5.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # ECMB plus Minnesota Archean subprovince
    #name = "AllMinnesotaVartCryst"
    #test = (data[:Key] .== "ECMB") .| (data[:Key] .== "MinnesotaArchean")
    #dt = 10 # time step size in Myr
    #ModelUncertainty = 35.0 # Ma
    #dTmax = 5.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # # ECMB, Minnesota, and Keewenawan
    # name = "AllMidwestVartCryst"
    # test = (data[:Key] .== "ECMB") .| (data[:Key] .== "MinnesotaArchean") .| (data[:Key] .== "MidcontinentRift")
    # dt = 10 # time step size in Myr
    # ModelUncertainty = 35. # Ma
    # dTmax = 5.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # # Literature samples from DeLucia et al. 2017 (Geology), Ozark Plateau
    # name = "OzarkInversionSphericalVartCryst"
    # test = data[:Geological_Province] .== "Ozark Plateau"
    # dt = 5 # time step size in Myr
    # dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

    # # Literature samples from Guenthner et al. 2017 (Tectonics), Fennoscandia
    # name = "FennoscandiaInversionSphericalVartCryst"
    # test = (data[:Geological_Province] .== "Fennoscandian Shield")
    # dt = 10 # time step size in Myr
    # dTmax = 10.0 # Maximum reheating/burial per model timestep (to prevent numerical overflow)

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

    AnnealedSigma = SimAnnealSigma.(1, HeAge_Ma_sigma, simannealparams)


## --- Test proscribed t-T paths with Neoproterozoic exhumation step

    # Generate T path to test
    # Tr = 150
    # T0 = 30
    # agePoints = Array{Float64}([tCryst, tCryst*29/30,   720, 580, 250,  0]) # Age (Ma)
    # TPoints  =  Array{Float64}([TCryst,        Tr+T0, Tr+T0,  T0,  70, 10]) # Temp. (C)
    # TSteps = linterp1s(agePoints,TPoints,ageSteps)

    # Plot t-T path
    #plot(ageSteps,TSteps,xflip=true)

    # Calculate model ages
    #CalcHeAges = Array{Float64}(undef, size(HeAge_Ma))
    #@time pr = CalcDamageAnnealing(dt,tSteps,TSteps)
    #@time for i=1:length(Halfwidth)
    #    first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
    #    CalcHeAges[i] = CalcZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i],diffusionparams)
    #end

    # Plot Comparison of results
    #p2 = plot(eU, CalcHeAges, seriestype=:scatter,label="Model")
    #plot!(p2, eU, HeAge_Ma, yerror=HeAge_Ma_sigma*2, seriestype=:scatter, label="Data")
    #xlabel!(p2,"eU"); ylabel!(p2,"Age (Ma)")
    #display(p2)

    # Check log likelihood
    #ll = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* HeAge_Ma_sigma.^2))
    #ll = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* AnnealedSigma.^2))


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

    simannealparams = lambda, InitialUncertainty, ModelUncertainty

    function MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints, simannealparams, diffusionparams)


        # Calculate model ages for initial proposal
        TSteps = linterp1s([agePoints[1:nPoints] ; boundary_agePoints ; unconf_agePoints],
                            [TPoints[1:nPoints] ; boundary_TPoints ; unconf_TPoints], ageSteps)
        CalcHeAges = Array{Float64}(undef, size(HeAge_Ma))
        pr = CalcDamageAnnealing(dt,tSteps,TSteps) # Damage annealing history
        for i=1:length(Halfwidth)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
            CalcHeAges[i] = CalcZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i], diffusionparams)
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
        nDist = zeros(nsteps, 1)
        acceptanceDist = zeros(nsteps, 1)

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
            pr = CalcDamageAnnealing(dt,tSteps,TSteps_prop)
            for i=1:length(Halfwidth)
                first_index = 1 + floor(Int64,(tCryst - CrystAge_Ma[i])/dt)
                CalcHeAges_prop[i] = CalcZrnHeAgeSpherical(dt,ageSteps[first_index:end],TSteps_prop[first_index:end],pr[first_index:end,first_index:end],Halfwidth[i],dr,U_ppm[i],Th_ppm[i],diffusionparams)
            end

            # Calculate log likelihood of proposal
            AnnealedSigma .= SimAnnealSigma.(n, HeAge_Ma_sigma, simannealparams)

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
            llDist[n] = sum(-(CalcHeAges - HeAge_Ma).^2 ./ (2 .* SimAnnealSigma.(nsteps, HeAge_Ma_sigma, simannealparams).^2)) # Recalculated to constant baseline
            nDist[n] = nPoints # Distribution of # of points
            HeAgeDist[:,n] = CalcHeAges # Distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            TStepDist[:,n] = TSteps # Distribution of T paths

        end
        return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
    end

    # Run Markov Chain
    (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(nPoints, maxPoints, agePoints, TPoints, unconf_agePoints, unconf_TPoints, boundary_agePoints, boundary_TPoints)

    # # Save results using JLD
    @save string(name, ".jld") ageSteps tSteps TStepDist burnin nsteps TCryst tCryst

    # Plot sample age-eU correlations
    h = scatter(eU,HeAgeDist[:,burnin:burnin+50], label="")
    plot!(h, eU, HeAge_Ma, yerror=HeAge_Ma_sigma, seriestype=:scatter, label="Data")
    xlabel!(h,"eU (ppm)"); ylabel!(h,"Age (Ma)")
    savefig(h,string(name,"Age-eU.pdf"))

    #


## --- Create image of paths

    # Resize the post-burnin part of the stationary distribution
    TStepDistResized = Array{Float64}(undef, 2001, size(TStepDist,2)-burnin)
    xq = collect(range(0,tCryst,length=2001))
    for i=1:size(TStepDist,2)-burnin
        TStepDistResized[:,i] = linterp1s(tSteps,TStepDist[:,i+burnin],xq)
    end

    # Calculate composite image
    tTimage = zeros(ceil(Int, TCryst)*2, size(TStepDistResized,1))
    yq = collect(0:0.5:TCryst)
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
