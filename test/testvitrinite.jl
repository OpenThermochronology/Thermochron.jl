## --- Test vitrinite reflectance kinetic model

    @test NielsenBasinRo() isa Thermochron.VitriniteReflectanceModel{Float64}

    dm = NielsenBasinRo()
    @test dm isa Thermochron.NielsenBasinRo{Float64}
    @test dm.Ro0 ≈ 0.2104
    @test dm.logA ≈ 60.9856
    @test length(dm.Ea) == 20
    @test length(dm.stoich) == 20
    @test sum(dm.stoich) ≈ 1.0

    dm2 = NielsenBasinRo(Float64; Ro0=0.25)
    @test dm2.Ro0 ≈ 0.25
    @test sum(dm2.stoich) ≈ 1.0

## --- Test VitriniteRo construction

    vro = VitriniteRo(Ro=0.6, Ro_sigma=0.05, agesteps=reverse(cntr(0:50)))
    @test vro isa VitriniteRo{Float64}
    show(vro)
    println()
    display(vro)

    @test Thermochron.value(vro) ≈ 0.6
    @test Thermochron.stdev(vro) ≈ 0.05

## --- Test vitrinite reflectance model ages

    # Isothermal residence at 80 C for 50 Ma
    vro = VitriniteRo(agesteps=reverse(cntr(0:50)))
    @test modelage(vro, fill(80, 50), NielsenBasinRo()) ≈ 0.6306865531922565

    # Burial-then-exhumation: 100 Ma (20 C) -> 40 Ma (150 C) -> present (90 C)
    agesteps = reverse(cntr(0:100))
    Tsteps = similar(collect(agesteps))
    for (i, a) in enumerate(agesteps)
        Tsteps[i] = a > 40 ? 20 + (100-a)/(100-40)*(150-20) : 150 + (40-a)/40*(90-150)
    end
    vro = VitriniteRo(agesteps=agesteps)
    @test modelage(vro, Tsteps, NielsenBasinRo()) ≈ 1.2622215580241707

    # Sanity checks: more heating should never yield lower predicted %Ro,
    # and predicted %Ro should never drop below the immature starting value
    vro = VitriniteRo(agesteps=reverse(cntr(0:50)))
    dm = NielsenBasinRo()
    Ro_cool = modelage(vro, fill(80, 50), dm)
    Ro_hot = modelage(vro, fill(100, 50), dm)
    @test Ro_hot > Ro_cool
    @test Ro_cool >= dm.Ro0
    @test Ro_hot >= dm.Ro0

    # 20 C isothermal for 28 Ma: mild but non-negligible maturation above Ro0
    vro = VitriniteRo(agesteps=reverse(cntr(0:28)))
    @test modelage(vro, fill(20., 28), NielsenBasinRo()) ≈ 0.2610983385906111

## --- Test model! dispatch for VitriniteRo

    vro = VitriniteRo(Ro=0.6, Ro_sigma=0.05, agesteps=reverse(cntr(0:50)))
    dm = NielsenBasinRo()
    Tsteps = fill(80.0, 50)

    chrons = Thermochron.Chronometer{Float64}[vro]
    damodels = Thermochron.Model{Float64}[dm]
    μcalc = zeros(1)
    σcalc = zeros(1)
    ll = Thermochron.model!(μcalc, σcalc, chrons, damodels, Tsteps)

    @test ll ≈ 1.88846083098508
    @test μcalc[1] ≈ 0.6306865531922565

## --- Test model! with mixed chronometer types

    apatite = ApatiteFT(age=10, age_sigma=1, agesteps=reverse(cntr(0:50)), dpar=2.16)
    aftm = Ketcham2007FC()
    vro = VitriniteRo(Ro=0.6, Ro_sigma=0.05, agesteps=reverse(cntr(0:50)))
    vrom = NielsenBasinRo()

    chrons = Thermochron.Chronometer{Float64}[vro, apatite]
    damodels = Thermochron.Model{Float64}[vrom, aftm]
    μcalc = zeros(2)
    σcalc = zeros(2)
    ll = Thermochron.model!(μcalc, σcalc, chrons, damodels, fill(80.0, 50))

    @test ll ≈ -264.0174839895607
    @test μcalc[1] ≈ 0.6306865531922565
    @test isfinite(μcalc[2]) && !isnan(μcalc[2])  # ApatiteFT computed something real; exact value not pinned here

## --- End of File
