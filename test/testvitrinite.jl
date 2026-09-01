using Test
using Thermochron

@testset "VitriniteRo construction" begin
    agesteps = collect(50.0:-1.0:0.0)
    vro = VitriniteRo(Float64; Ro=0.6, Ro_sigma=0.05, agesteps=agesteps)

    @test vro.Ro ≈ 0.6
    @test vro.Ro_sigma ≈ 0.05
    @test vro.offset ≈ 0.0
    @test vro.height ≈ 0.0
    @test Thermochron.value(vro) ≈ 0.6
    @test Thermochron.stdev(vro) ≈ 0.05
    @test Thermochron.agesteps_geol(vro) == agesteps
    @test issorted(Thermochron.tsteps_geol(vro))
end

@testset "NielsenBasinRo construction" begin
    dm = NielsenBasinRo()

    @test dm.Ro0 ≈ 0.2104
    @test dm.logA ≈ 60.9856
    @test length(dm.Ea) == 20
    @test length(dm.stoich) == 20
    # Stoichiometric factors must be normalized to sum to 1
    @test sum(dm.stoich) ≈ 1.0 atol=1e-12
    @test all(>=(0), dm.stoich)

    # Custom-parameter construction should also work and stay unnormalized-input-safe
    dm2 = NielsenBasinRo(Float64; Ro0=0.25)
    @test dm2.Ro0 ≈ 0.25
    @test sum(dm2.stoich) ≈ 1.0 atol=1e-12
end

@testset "modelage(VitriniteRo) forward model" begin
    # Reference values validated against a standalone reimplementation of
    # the original calculate_final_basin_ro script (see conversation record);
    # both isothermal and burial-then-exhumation paths agreed with the
    # original script to <1% before these exact values were locked in here.
    dm = NielsenBasinRo()

    # Isothermal at 80 C for 50 Ma
    agesteps_iso = collect(50.0:-1.0:0.0)
    vro_iso = VitriniteRo(Float64; Ro=NaN, Ro_sigma=NaN, agesteps=agesteps_iso)
    Tsteps_iso = fill(80.0, length(agesteps_iso))
    Ro_iso = Thermochron.modelage(vro_iso, Tsteps_iso, dm)
    @test isapprox(Ro_iso, 0.6314790269061636, atol=1e-8)

    # Burial-then-exhumation: 100 Ma (20 C) -> 40 Ma (150 C) -> present (90 C)
    agesteps_burial = collect(100.0:-2.0:0.0)
    vro_burial = VitriniteRo(Float64; Ro=NaN, Ro_sigma=NaN, agesteps=agesteps_burial)
    Tsteps_burial = similar(agesteps_burial)
    for (i, a) in enumerate(agesteps_burial)
        Tsteps_burial[i] = a > 40 ? 20 + (100-a)/(100-40)*(150-20) : 150 + (40-a)/40*(90-150)
    end
    Ro_burial = Thermochron.modelage(vro_burial, Tsteps_burial, dm)
    @test isapprox(Ro_burial, 1.2689771793216664, atol=1e-8)

    # Sanity checks independent of the exact reference values:
    # more heating should never produce a lower predicted %Ro
    Tsteps_hotter = Tsteps_iso .+ 20.0
    Ro_hotter = Thermochron.modelage(vro_iso, Tsteps_hotter, dm)
    @test Ro_hotter > Ro_iso

    # Predicted %Ro should never fall below the immature starting value Ro0
    @test Ro_iso >= dm.Ro0
    @test Ro_burial >= dm.Ro0
end

@testset "model! dispatches VitriniteRo" begin
    dm = NielsenBasinRo()
    agesteps = collect(50.0:-1.0:0.0)
    Tsteps = fill(80.0, length(agesteps))

    # Ro/Ro_sigma must be real (non-NaN) here since norm_ll needs a finite
    # observation to compare the model prediction against
    vro = VitriniteRo(Float64; Ro=0.6, Ro_sigma=0.05, agesteps=agesteps)

    chrons = Thermochron.Chronometer{Float64}[vro]
    damodels = Thermochron.Model{Float64}[dm]
    μcalc = zeros(1)
    σcalc = zeros(1)

    ll = Thermochron.model!(μcalc, σcalc, chrons, damodels, Tsteps)

    @test isfinite(ll)
    @test isapprox(ll, 1.878607913357523, atol=1e-8)
    @test isapprox(μcalc[1], 0.6314790269061636, atol=1e-8)
end

@testset "model! with mixed chronometer types" begin
    # Confirms VitriniteRo coexists with another chronometer type in the
    # same inversion without breaking Tsteps slicing or log-likelihood summation
    dm_vro = NielsenBasinRo()
    agesteps = collect(50.0:-1.0:0.0)
    Tsteps = fill(80.0, length(agesteps))
    vro = VitriniteRo(Float64; Ro=0.6, Ro_sigma=0.05, agesteps=agesteps)

    aft = ApatiteFT(Float64; age=10.0, age_sigma=1.0, rmr0=0.83, agesteps=agesteps)
    dm_aft = Ketcham2007FC()

    chrons = Thermochron.Chronometer{Float64}[vro, aft]
    damodels = Thermochron.Model{Float64}[dm_vro, dm_aft]
    μcalc = zeros(2)
    σcalc = zeros(2)

    ll = Thermochron.model!(μcalc, σcalc, chrons, damodels, Tsteps)

    @test isfinite(ll)
    @test isapprox(μcalc[1], 0.6314790269061636, atol=1e-8)  # VitriniteRo unaffected by neighboring ApatiteFT
    @test !isnan(μcalc[2])                                    # ApatiteFT still computed
end
