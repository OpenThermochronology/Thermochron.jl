# Test fission track annealing

    @test Ketcham1999FC() isa Thermochron.ApatiteAnnealingModel{Float64}
    @test Ketcham2007FC() isa Thermochron.ApatiteAnnealingModel{Float64}
    @test Yamada2005PC() isa Thermochron.ZirconAnnealingModel{Float64}

## --- Ketcham et al. 1999 Fanning Curvilinear apatite
    am = Ketcham1999FC()
    @test am isa Thermochron.Ketcham1999FC{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.969772147771495
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.964689584499856
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.8757101462242018
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.6865397486321438
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.0

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0017396519912747071, 0.003298079407422488, 0.004794716652673411, 0.006252588352279524, 0.007682296289640212, 0.009089953816004663, 0.010479534909963182, 0.01185382650735338, 0.013214888894288724, 0.014564305574711722]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈  [660.4429511024003, 1395.0963084362782, 2160.6369846963353, 2946.952051746238, 3749.0850835490833, 4564.053074060813, 5389.851797105461, 6225.0371840101425, 7068.5169276849765, 7919.434602704827]
 
    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9516354364343921
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9435033351997696 
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.801136233958723 
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.3210110092650815
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Ketcham et al. 2007 Simplified Curvilinear apatite

    am = Ketcham2007FC()
    @test am isa Thermochron.Ketcham2007FC{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9749628064438557
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.9701964063931915
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.8759898250381366
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.5928200145772959
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.016540366949720275

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0015478065155001086, 0.002954611923287363, 0.00431269182539665, 0.005640066461673555, 0.006945066485586564, 0.008232508754276639, 0.009505531052511628, 0.010766337684274287, 0.012016560593707997, 0.013257455657955832]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [646.5878531291968, 1354.896862856624, 2088.5372542290797, 2839.1277381635086, 3602.579481389471, 4376.439427402521, 5159.063447504911, 5949.269302029918, 6746.164254698982, 7549.049388886254]
    
    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9599404903101693 
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9523142502291065
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.8015837200610186
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.07551154545163152
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Ketcham et al. 2007 Simplified Curvilinear apatite

    am = Yamada2005PC()
    @test am isa Thermochron.Yamada2005PC{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9964516173825159
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.995078688092618
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.940753376497744
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.587467392379713
    @test Thermochron.reltracklength(1, 500, am) ≈ 4.8685772129709805e-21

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0020577410195178228, 0.0041154820390356525, 0.006173223058553454, 0.00823096407807132, 0.010288705097588996, 0.01234644611710693, 0.014404187136625298, 0.016461928156142672, 0.018519669175660315, 0.02057741019517803]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [222.98141544261495, 445.9628308852307, 668.9442463278432, 891.9256617704631, 1114.907077213062, 1337.8884926556889, 1560.8699080982742, 1783.8513235409293, 2006.8327389835247, 2229.814154426128] 

    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9955645217281448
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9938483601157724
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.92594172062218 
    @test Thermochron.reltrackdensity(1, 200, am) ≈  0.48433424047464124
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Test "multikinetic" rmr0 model

    F = [1.75, 1.76, 1.64, 1.66, 1.64, 1.72, 1.72, 1.7, 1.66, 1.66]
    Cl = [0.01, 0.01, 0.0, 0.01, 0.0, 0.0, 0.01, 0.01, 0.0, 0.01]
    OH = [0.24, 0.23, 0.36, 0.33, 0.36, 0.28, 0.27, 0.29, 0.34, 0.33]
    @test rmr0model.(F, Cl, OH) ≈ [0.8573573076438294, 0.857484193068046, 0.8569770580132927, 0.856210256388646, 0.8569770580132927, 0.857991689484071, 0.8569759763207477, 0.8567211907298894, 0.8572313896305666, 0.856210256388646]

    rmr0 = rmr0model.(F, Cl, OH)
    rmr = Thermochron.reltracklength.(1:10, 95, Ketcham2007FC())
    @test Thermochron.rlr.(rmr, rmr0) ≈ [0.7390142328562013, 0.6849516733686434, 0.6428585390459669, 0.6061834471512542, 0.5490365748013377, 0.44555008284518977, 0.33418938122916036, 0.0, 0.0, 0.0]

## --- Test zircon fission track model ages

    zircon = ZirconFT(agesteps=reverse(cntr(0:100)))
    @test zircon isa ZirconFT{Float64}
    show(zircon)
    display(zircon)

    # Isothermal residence
    @test modelage(zircon, fill(0, 100), Yamada2005PC()) ≈ 99.0263284502223
    @test modelage(zircon, fill(50, 100), Yamada2005PC()) ≈ 95.55787511259514
    @test modelage(zircon, fill(75, 100), Yamada2005PC()) ≈ 91.39371027871266 
    @test modelage(zircon, fill(100, 100), Yamada2005PC()) ≈ 84.31219878239003

    # Linear cooling
    @test modelage(zircon, reverse(1:100), Yamada2005PC()) ≈ 95.84543737681496

    # As above but longer history
    zircon = ZirconFT(agesteps=reverse(cntr(0:200)))
    @test zircon isa ZirconFT{Float64}

    @test modelage(zircon, reverse(1:200), Yamada2005PC()) ≈ 158.732395422251
    @test modelage(zircon, reverse(1:200)./2, Yamada2005PC()) ≈ 190.5424305134446
    @test modelage(zircon, reverse(1:200).*2, Yamada2005PC()) ≈ 84.2529410694908

    # Fish Canyon Tuff zircon example
    zircon = ZirconFT(age=27, age_sigma=3, agesteps=reverse(cntr(0:28)))
    @test modelage(zircon, fill(20., 28), Yamada2005PC()) ≈ 27.600459673764092
    @test Thermochron.model_ll(zircon, fill(20., 28), Yamada2005PC()) ≈ -2.0375814785292756

## --- Test apatite fission track model ages

    apatite = ApatiteFT(agesteps=reverse(cntr(0:100)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294
    show(apatite)
    display(apatite)

    # Isothermal residence
    @test modelage(apatite, fill(0, 100), Ketcham1999FC()) ≈ 89.47899236366592
    @test modelage(apatite, fill(0, 100), Ketcham2007FC()) ≈ 91.18888272552469

    @test modelage(apatite, fill(50, 100), Ketcham1999FC()) ≈ 71.75968375010271
    @test modelage(apatite, fill(50, 100), Ketcham2007FC()) ≈ 74.20874464530966

    @test modelage(apatite, fill(75, 100), Ketcham1999FC()) ≈ 22.088815691960797
    @test modelage(apatite, fill(75, 100), Ketcham2007FC()) ≈ 21.6047849747377

    @test modelage(apatite, fill(100, 100), Ketcham1999FC()) ≈ 0.4682983384208703 
    @test modelage(apatite, fill(100, 100), Ketcham2007FC()) ≈ 0.42338708872671615

    # Linear cooling
    @test modelage(apatite, reverse(1:100), Ketcham1999FC()) ≈ 66.15470736807784 
    @test modelage(apatite, reverse(1:100), Ketcham2007FC()) ≈ 67.87471588034019

    # As above but longer history
    apatite = ApatiteFT(agesteps=reverse(cntr(0:200)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294

    @test modelage(apatite, reverse(1:200), Ketcham1999FC()) ≈ 66.15470736807784
    @test modelage(apatite, reverse(1:200), Ketcham2007FC()) ≈ 67.87471588034019

    @test modelage(apatite, reverse(1:200)./2, Ketcham1999FC()) ≈ 124.23025599587466
    @test modelage(apatite, reverse(1:200)./2, Ketcham2007FC()) ≈ 127.60575486065872

    @test modelage(apatite, reverse(1:200).*2, Ketcham1999FC()) ≈ 34.92192150860748
    @test modelage(apatite, reverse(1:200).*2, Ketcham2007FC()) ≈ 35.79582640096202 

    # Fish Canyon Tuff apatite example
    apatite = ApatiteFT(age=25, age_sigma=3, agesteps=reverse(cntr(0:28)), dpar=2.16)
    @test modelage(apatite, fill(20., 28), Ketcham2007FC()) ≈ 25.25247092840902
    @test Thermochron.model_ll(apatite, fill(20., 28), Ketcham2007FC()) ≈ -2.0210920201889886

## --- Test track lengths

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:20)), F=1.75, Cl=0.01, OH=0.24)
    show(track)
    display(track)

    l, σ = modellength(track, fill(75, 20), Ketcham1999FC()) .* 16.38
    @test l ≈ 12.432312726056672
    @test σ ≈  0.5736749779212151
    @test round.(track.r, sigdigits=4) ≈ [0.7101, 0.7134, 0.7167, 0.7202, 0.7238, 0.7275, 0.7314, 0.7354, 0.7396, 0.7441, 0.7488, 0.7539, 0.7593, 0.7653, 0.7718, 0.7792, 0.7877, 0.7979, 0.8111, 0.8311]
    l, σ = modellength(track, fill(75, 20), Ketcham2007FC()) .* 16.38
    @test l ≈ 12.610086362166317
    @test σ ≈  0.6062929557510826
    @test round.(track.r, sigdigits=4) ≈ [0.7176, 0.7212, 0.725, 0.7289, 0.7329, 0.737, 0.7413, 0.7457, 0.7503, 0.7552, 0.7603, 0.7657, 0.7716, 0.7779, 0.7848, 0.7925, 0.8014, 0.8119, 0.8254, 0.8454]

    l, σ = modellength(track, fill(100, 20), Ketcham1999FC()) .* 16.38
    @test l ≈ 10.874500968379554 
    @test σ ≈  0.6344191502094181
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4544, 0.55, 0.6184, 0.6876]
    l, σ = modellength(track, fill(100, 20), Ketcham2007FC()) .* 16.38
    @test l ≈ 10.982218818618923
    @test σ ≈  0.608187681868191
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4888, 0.6017, 0.6895]

    l, σ = modellength(track, reverse(1:20).*5, Ketcham1999FC()) .* 16.38
    @test l ≈ 14.363662652747362
    @test σ ≈  1.1822475087376938
    @test round.(track.r, sigdigits=4) ≈ [0.6124, 0.6806, 0.726, 0.7602, 0.7876, 0.8104, 0.8298, 0.8466, 0.8614, 0.8744, 0.886, 0.8963, 0.9056, 0.914, 0.9216, 0.9286, 0.935, 0.941, 0.9469, 0.9536]
    l, σ = modellength(track, reverse(1:20).*5, Ketcham2007FC()) .* 16.38
    @test l ≈ 14.558435101305054
    @test σ ≈  1.1716416735159647
    @test round.(track.r, sigdigits=4) ≈ [0.5945, 0.6823, 0.7347, 0.7723, 0.8014, 0.825, 0.8447, 0.8615, 0.8759, 0.8885, 0.8996, 0.9094, 0.9181, 0.926, 0.933, 0.9394, 0.9452, 0.9506, 0.956, 0.9619]

    l, σ =  modellength(track, reverse(1:20).*10, Ketcham1999FC()) .* 16.38 
    @test l ≈ 14.357128787814558
    @test σ ≈  1.1841227131850056
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6597, 0.7484, 0.8021, 0.8403, 0.8693, 0.8921, 0.9105, 0.9255, 0.9381, 0.9494]
    l, σ = modellength(track, reverse(1:20).*10, Ketcham2007FC()) .* 16.38
    @test l ≈ 14.538333058948192
    @test σ ≈  1.1899130108299063
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6564, 0.7592, 0.8162, 0.8549, 0.8834, 0.9053, 0.9225, 0.9364, 0.9479, 0.9581]

## -- Check the mean track length of a modelled Fish Canyon Apatite to be 15.35 +/- 0.06 um

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:28)), dpar=2.16)
    l, σ = modellength(track, fill(20., 28), Ketcham2007FC()) .* 16.38
    @test l ≈ 15.35 atol=0.06
    @test l ≈ 15.376461527029509 
    @test σ ≈  0.09967577151793756

    @test Thermochron.model_ll(track, fill(20., 28), Ketcham2007FC()) ≈ 1.0651088634382597

## --- Test c-axis equivalent model lengths

    @test lcmod(track) ≈ lcmod(15, 35) ≈ 15.405678663327869

    l = [13.9, 12.73, 13.98, 12.75, 13.77, 15.93, 13.99, 15.41, 13.23, 12.04, 13.01, 10.94]
    θ = [59.19, 78.45, 11.98, 68.73, 41.21, 35.05, 45.08, 41.42, 30.04, 53.9, 62.06, 70.21]
    @test lcmod.(l, θ) ≈ [14.980169300957634, 14.430422822748183, 14.085946809091633, 14.356417336926649, 14.570689505578352, 16.159195153315775, 14.812278052422833, 15.815567116293908, 13.865575929923889, 13.644277619942535, 14.435368624348033, 13.236076028425188] 

## --- 