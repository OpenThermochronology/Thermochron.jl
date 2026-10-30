## --- Test vitrinite reflectance kinetic model

    @test NielsenBasinRo() isa Thermochron.VitriniteReflectanceModel{Float64}

    dm = NielsenBasinRo()
    @test dm isa Thermochron.NielsenBasinRo{Float64}
    @test dm.Ro0 ≈ 0.2104
    @test dm.logA ≈ 60.9856
    @test length(dm.Ea) == length(dm.fraction) == 20
    @test sum(dm.fraction) ≈ 1.0

    dm2 = NielsenBasinRo(Ro0=0.25)
    @test dm2.Ro0 ≈ 0.25
    @test sum(dm2.fraction) ≈ 1.0

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
    @test modelreflectance(vro, fill(80, 50), NielsenBasinRo()) ≈ 0.543989018090385 

    # Burial-then-exhumation: 100 Ma (20 C) -> 40 Ma (150 C) -> present (90 C)
    agesteps = reverse(cntr(0:100))
    Tsteps = similar(collect(agesteps))
    for (i, a) in enumerate(agesteps)
        Tsteps[i] = a > 40 ? 20 + (100-a)/(100-40)*(150-20) : 150 + (40-a)/40*(90-150)
    end
    vro = VitriniteRo(agesteps=agesteps)
    @test modelreflectance(vro, Tsteps, NielsenBasinRo()) ≈ 0.8705766330673169 

    # Sanity checks: more heating should never yield lower predicted %Ro,
    # and predicted %Ro should never drop below the immature starting value
    vro = VitriniteRo(agesteps=reverse(cntr(0:50)))
    dm = NielsenBasinRo()
    Ro_cool = modelreflectance(vro, fill(80, 50), dm)
    Ro_hot = modelreflectance(vro, fill(100, 50), dm)
    @test Ro_hot > Ro_cool
    @test Ro_cool >= dm.Ro0
    @test Ro_hot >= dm.Ro0

    # 20 C isothermal for 28 Ma: mild but non-negligible maturation above Ro0
    vro = VitriniteRo(agesteps=reverse(cntr(0:28)))
    @test modelreflectance(vro, fill(20., 28), NielsenBasinRo()) ≈ 0.2594903705509342

    # Constant heating from 10°C at 2.8 C/Ma,
    # as in Fig. 5 of Nielsen et al. 2017: 
    Tsteps = 10:2.8:290
    vro = VitriniteRo(;agesteps = length(Tsteps)-1:-1:0.0)
    @test modelreflectance(vro, Tsteps, NielsenBasinRo()) ≈ 7.85173749893207

    Tsteps = 10:2.8:201
    vro = VitriniteRo(;agesteps = length(Tsteps)-1:-1:0.0)
    @test modelreflectance(vro, Tsteps, NielsenBasinRo()) ≈ 2.748360568398062

    Tsteps = 10:2.8:150
    vro = VitriniteRo(;agesteps = length(Tsteps)-1:-1:0.0)
    @test modelreflectance(vro, Tsteps, NielsenBasinRo()) ≈ 0.8046313060705653


## --- Test model! dispatch for VitriniteRo

    vro = VitriniteRo(Ro=0.6, Ro_sigma=0.05, agesteps=reverse(cntr(0:50)))
    dm = NielsenBasinRo()
    Tsteps = fill(80.0, 50)

    chrons = Thermochron.Chronometer{Float64}[vro]
    damodels = Thermochron.Model{Float64}[dm]
    μcalc = zeros(1)
    σcalc = zeros(1)
    ll = Thermochron.model!(μcalc, σcalc, chrons, damodels, Tsteps)

    @test ll ≈ 1.4493477214534756
    @test μcalc[1] ≈ 0.543989018090385 

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

    @test ll ≈ -264.4565970990923
    @test μcalc[1] ≈ 0.543989018090385
    @test isfinite(μcalc[2]) && !isnan(μcalc[2])  # ApatiteFT computed something real; exact value not pinned here

## --- End of File
