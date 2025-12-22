## --- Prepare problem

    # Read in data from file using StatGeochem
    datapath = joinpath("..", "examples", "exampledata", "minnesotazrn.csv")
    ds = importdataset(datapath, ',', importas=:Tuple);

    using LinearAlgebra
    BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

    params = (
        burnin = 350,               # [n] How long should we wait for MC to converge (become stationary)
        nsteps = 350,               # [n] How many steps of the Markov chain should we run after burn-in?
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
        simplified = false,         # Prefer simpler tT paths?
        dynamicsigma = true,        # Update model uncertainties?
        # Diffusivity and annealing models
        zdm = ZRDAAM(), 
        adm = RDAAM(), 
        zftm = Yamada2007PC(), 
        mftm = Jones2021FA(), 
        aftm = Ketcham2007FC(),
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 25 Ma for all age-based chronometers.
        σmodel = 25.0,              # [Ma] assumed model uncertainty (resampled if dynamicsigma)
        T0annealing = 1,            # [unitless] initial annealing "temperature"
    )

    # Populate data NamedTuple from imported dataset
    data = deepcopy(ds)
    data.raw_He_age_sigma_Ma .= 0.1*data.raw_He_age_Ma

    # Crystallization ages and start time
    tinit = ceil(maximum(data.crystallization_age_Ma)/params.dt) * params.dt
    params = (params...,
        tinit = tinit,
        agesteps = Array{Float64}(tinit-params.dt/2 : -params.dt : 0+params.dt/2),
        tsteps = Array{Float64}(0+params.dt/2 : params.dt : tinit-params.dt/2),
    )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = [params.tnow, params.tinit],   # [Ma] Final and initial time
        T₀ = [params.Tnow, params.Tinit],          # [C] Final and initial temperature
        ΔT = [params.ΔTnow, params.ΔTinit],        # [C] Final and initial temperature range (positive or negative)
        tboundary = :reflecting, # Reflecting time boundary conditions
        Tboundary = :reflecting, # Reflecting temperature boundary conditions
    )

    # Example constraint box
    constraint = Constraint(
        agedist = [  Uniform(500,541),],  # [Ma] Age distribution
        Tdist =   [     Uniform(0,50),],  # [C] Temperature distribution
    )

## --- Invert for t-T path via MCMC

    # Run Markov Chain
    @time "\nCompiling MCMC" MCMC(data, params, boundary, constraint)
    @time "\nRunning MCMC" tT = MCMC(data, params, boundary, constraint; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
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
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT: $(mean(tT.jTdist))"

## --- MCMC with Detail interval
    detail = DetailInterval(
        agemin = 0, # Youngest end of detail interval
        agemax = 1000, # Oldest end of detail interval
        minpoints = 3, # Minimum number of points in detail interval
    )
    @time "\nMCMC with Detail interval" tT = MCMC(data, params, boundary, constraint, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 200
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -300 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test mean(tT.jtdist) ≈ params.tinit/60
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test mean(tT.jTdist) ≈ params.Tinit/60
    @info "Mean σjT: $(mean(tT.jTdist))"

## --- MCMC with Detail interval & dynamicjumping
    params = (params...,
        dynamicjumping=true
    )
    @time "\nMCMC with Detail interval & dynamicjumping" tT = MCMC(data, params, boundary, constraint, detail; liveplot)

    @test isa(tT.Tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.Tinit
    @test nanminimum(tT.Tpointdist) >= params.Tnow

    @test isa(tT.tpointdist, AbstractMatrix)
    @test nanmaximum(tT.Tpointdist) <= params.tinit
    @test nanminimum(tT.Tpointdist) >= 0

    @test isa(tT.resultdist, AbstractMatrix)
    abserr = abs(sum(nanmean(tT.resultdist, dims=2) - data.raw_He_age_Ma)/length(data.raw_He_age_Ma))
    @test 0 < abserr < 200
    @info "Mean absolute error: $abserr"

    @test isa(tT.lldist, AbstractVector)
    llmean = mean(tT.lldist)
    @test -300 < llmean < 0
    @info "Mean ll: $llmean"

    @test isa(tT.acceptancedist, AbstractVector{Bool})
    @test isapprox(mean(tT.acceptancedist), 0.5, atol=0.45)
    @info "Mean acceptance rate: $(mean(tT.acceptancedist))"

    @test isa(tT.ndist, AbstractVector{Int})
    @test minimum(tT.ndist) >= params.minpoints
    @test maximum(tT.ndist) <= params.maxpoints
    @info "Mean npoints: $(mean(tT.ndist))"

    @test params.dt < mean(tT.jtdist) < params.tinit
    @info "Mean σjt: $(mean(tT.jtdist))"

    @test 0 < mean(tT.jTdist) < params.Tinit
    @info "Mean σjT: $(mean(tT.jTdist))"

## --- End of File