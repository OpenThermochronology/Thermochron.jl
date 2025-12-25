## --- Test various inversion utilities
    # Test simulated annealing
    @test Thermochron.simannealT(1, 100., 10 / 10^5) ≈ 100.99000049998334
    @test Thermochron.simannealT(10^5, 100., 10 / 10^5) ≈ 1.0045399929762484

    # Test collectto!
    buffer = rand(20)
    a = rand(3)
    b = rand(4)
    c = rand(2)
    @test Thermochron.collectto!(buffer, a, b, c) == vcat(a,b,c)

    # Test pointsininterval
    @test Thermochron.pointsininterval(collect(1:10), 10, 0.5, 5.5, 0:10) == 5
    @test Thermochron.pointsininterval(collect(1:10), 10, 0.5, 5.5, 0.5:9.5) == 5

    # Test diff_ll
    @test Thermochron.diff_ll(1:10, 1, 1) == -6.238324625039508
    @test Thermochron.diff_ll(2*(1:10), 1, 1) ≈ -16.569194805083374
    @test Thermochron.diff_ll(2*(1:10), 1, 0.5) ≈ -34.048659003138276
    x = rand(100)
    d = diff(x)
    @test Thermochron.diff_ll(x, 0, 1) ≈ sum(logccdf.(Normal(0,1),d[d.>0]))

    # Test boundary condition functions
    @test Thermochron.reflecting(101, 0, 100) == 99
    @test Thermochron.reflecting(101, 1, 100) == 99
    @test Thermochron.reflecting(-1, 0, 100) == 1
    @test Thermochron.reflecting(-1, 1, 100) == 3
    @test Thermochron.reflecting(150, 0, 100) == 50
    @test Thermochron.reflecting(200, 0, 100) == 0
    @test Thermochron.reflecting(201, 0, 100) == 1
    @test Thermochron.reflecting(-101, 0, 100) == 99

    @test Thermochron.periodic(101, 0, 100) == 1
    @test Thermochron.periodic(101, 1, 100) == 2
    @test Thermochron.periodic(-1, 0, 100) == 99
    @test Thermochron.periodic(-1, 1, 100) == 98
    @test Thermochron.periodic(150, 0, 100) == 50
    @test Thermochron.periodic(200, 0, 100) == 0
    @test Thermochron.periodic(201, 0, 100) == 1
    @test Thermochron.periodic(-101, 0, 100) == 99

    @test Thermochron.hard(101, 0, 100) == 100
    @test Thermochron.hard(101, 1, 100) == 100
    @test Thermochron.hard(-1, 0, 100) == 0
    @test Thermochron.hard(-1, 1, 100) == 1
    @test Thermochron.hard(150, 0, 100) == 100
    @test Thermochron.hard(200, 0, 100) == 100
    @test Thermochron.hard(201, 0, 100) == 100
    @test Thermochron.hard(-101, 0, 100) == 0

## --- Test Constraint, Boundary, and TTPath types
    constraint = Constraint(
        agedist = [Uniform(500,580),],  # [Ma] Age distribution
        Tdist =   [   Uniform(0,50),],  # [C] Temperature distribution
    )
    @test constraint isa Constraint{Float64}

    boundary = Boundary(
        agepoints = [0.0, 1000.0],  # [Ma] Final and initial time
        T₀ = [0.0, 400.0],          # [C] Final and initial temperature
        ΔT = [50, -100.0],          # [C] Final and initial temperature range (positive or negative)
        tboundary = :reflecting,    # Reflecting time boundary conditions
        Tboundary = :reflecting,    # Reflecting temperature boundary conditions
    )
    @test boundary isa Boundary{Float64}

    agesteps = reverse(1:10:1000.)
    path = Thermochron.TTPath(agesteps, constraint, boundary, DetailInterval(), 50)
    @test path isa Thermochron.TTPath{Float64}
    @test path.agepoints == path.agepointsₚ == zeros(50)
    @test path.Tpoints == path.Tpointsₚ == zeros(50)
    @test path.Tsteps == zeros(100)

    Thermochron.initialproposal!(path, 10)
    @test path.agepoints != path.agepointsₚ
    @test path.Tpoints != path.Tpointsₚ
    a,T = copy(path.agepoints), copy(path.Tpoints)

    Thermochron.resetproposal!(path)
    @test path.agepoints == path.agepointsₚ == a
    @test path.Tpoints == path.Tpointsₚ == T

    Thermochron.initialproposal!(path, 10)
    @test path.agepoints != path.agepointsₚ
    @test path.Tpoints != path.Tpointsₚ
    a,T = copy(path.agepointsₚ), copy(path.Tpointsₚ)

    Thermochron.acceptproposal!(path)
    @test path.agepoints == path.agepointsₚ == a
    @test path.Tpoints == path.Tpointsₚ == T

## --- End of File
