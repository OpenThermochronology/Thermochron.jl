# Test fission track annealing

    @test FCKetcham1999 isa Thermochron.FanningCurvilinear{Float64}
    @test FCKetcham2007 isa Thermochron.SimplifiedCurvilinear{Float64}

## --- Ketcham et al. 1999
    am = FCKetcham1999
    @test am isa Thermochron.FanningCurvilinear{Float64}

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

## --- Ketcham et al. 2007

    am = FCKetcham2007
    @test am isa Thermochron.SimplifiedCurvilinear{Float64}

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

## --- Test "multikinetic" rmr0 model

    F = [1.75, 1.76, 1.64, 1.66, 1.64, 1.72, 1.72, 1.7, 1.66, 1.66]
    Cl = [0.01, 0.01, 0.0, 0.01, 0.0, 0.0, 0.01, 0.01, 0.0, 0.01]
    OH = [0.24, 0.23, 0.36, 0.33, 0.36, 0.28, 0.27, 0.29, 0.34, 0.33]
    @test rmr0model.(F, Cl, OH) ≈ [0.8573573076438294, 0.857484193068046, 0.8569770580132927, 0.856210256388646, 0.8569770580132927, 0.857991689484071, 0.8569759763207477, 0.8567211907298894, 0.8572313896305666, 0.856210256388646]

    rmr0 = rmr0model.(F, Cl, OH)
    rmr = Thermochron.reltracklength.(1:10, 95, FCKetcham2007)
    @test Thermochron.rlr.(rmr, rmr0) ≈ [0.7390142328562013, 0.6849516733686434, 0.6428585390459669, 0.6061834471512542, 0.5490365748013377, 0.44555008284518977, 0.33418938122916036, 0.0, 0.0, 0.0]

## --- Test fission track model ages

    apatite = ApatiteFT(agesteps=reverse(cntr(0:100)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294
    show(apatite)
    display(apatite)

    # Isothermal residence
    @test modelage(apatite, fill(0, 100), FCKetcham1999) ≈ 89.47899236366592
    @test modelage(apatite, fill(0, 100), FCKetcham2007) ≈ 91.18888272552469

    @test modelage(apatite, fill(50, 100), FCKetcham1999) ≈ 71.75968375010271
    @test modelage(apatite, fill(50, 100), FCKetcham2007) ≈ 74.20874464530966

    @test modelage(apatite, fill(75, 100), FCKetcham1999) ≈ 22.088815691960797
    @test modelage(apatite, fill(75, 100), FCKetcham2007) ≈ 21.6047849747377

    @test modelage(apatite, fill(100, 100), FCKetcham1999) ≈ 0.4682983384208703 
    @test modelage(apatite, fill(100, 100), FCKetcham2007) ≈ 0.42338708872671615

    # Linear cooling
    @test modelage(apatite, reverse(1:100), FCKetcham1999) ≈ 66.15470736807784 
    @test modelage(apatite, reverse(1:100), FCKetcham2007) ≈ 67.87471588034019

    apatite = ApatiteFT(agesteps=reverse(cntr(0:200)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294

    # As above but longer history
    @test modelage(apatite, reverse(1:200), FCKetcham1999) ≈ 66.15470736807784
    @test modelage(apatite, reverse(1:200), FCKetcham2007) ≈ 67.87471588034019

    @test modelage(apatite, reverse(1:200)./2, FCKetcham1999) ≈ 124.23025599587466
    @test modelage(apatite, reverse(1:200)./2, FCKetcham2007) ≈ 127.60575486065872

    @test modelage(apatite, reverse(1:200).*2, FCKetcham1999) ≈ 34.92192150860748
    @test modelage(apatite, reverse(1:200).*2, FCKetcham2007) ≈ 35.79582640096202 

    apatite = ApatiteFT(age=25, age_sigma=3, agesteps=reverse(cntr(0:28)), dpar=2.16)
    @test modelage(apatite, fill(20., 28), FCKetcham2007) ≈ 25.25247092840902
    @test Thermochron.model_ll(apatite, fill(20., 28), FCKetcham2007) ≈ -2.0210920201889886

## --- Test track lengths

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:20)), F=1.75, Cl=0.01, OH=0.24)
    show(track)
    display(track)

    l, σ = modellength(track, fill(75, 20), FCKetcham1999) .* 16.38
    @test l ≈ 12.432312726056672
    @test σ ≈  0.5736749779212151
    @test round.(track.r, sigdigits=4) ≈ [0.7101, 0.7134, 0.7167, 0.7202, 0.7238, 0.7275, 0.7314, 0.7354, 0.7396, 0.7441, 0.7488, 0.7539, 0.7593, 0.7653, 0.7718, 0.7792, 0.7877, 0.7979, 0.8111, 0.8311]
    l, σ = modellength(track, fill(75, 20), FCKetcham2007) .* 16.38
    @test l ≈ 12.610086362166317
    @test σ ≈  0.6062929557510826
    @test round.(track.r, sigdigits=4) ≈ [0.7176, 0.7212, 0.725, 0.7289, 0.7329, 0.737, 0.7413, 0.7457, 0.7503, 0.7552, 0.7603, 0.7657, 0.7716, 0.7779, 0.7848, 0.7925, 0.8014, 0.8119, 0.8254, 0.8454]

    l, σ = modellength(track, fill(100, 20), FCKetcham1999) .* 16.38
    @test l ≈ 10.874500968379554 
    @test σ ≈  0.6344191502094181
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4544, 0.55, 0.6184, 0.6876]
    l, σ = modellength(track, fill(100, 20), FCKetcham2007) .* 16.38
    @test l ≈ 10.982218818618923
    @test σ ≈  0.608187681868191
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4888, 0.6017, 0.6895]

    l, σ = modellength(track, reverse(1:20).*5, FCKetcham1999) .* 16.38
    @test l ≈ 14.363662652747362
    @test σ ≈  1.1822475087376938
    @test round.(track.r, sigdigits=4) ≈ [0.6124, 0.6806, 0.726, 0.7602, 0.7876, 0.8104, 0.8298, 0.8466, 0.8614, 0.8744, 0.886, 0.8963, 0.9056, 0.914, 0.9216, 0.9286, 0.935, 0.941, 0.9469, 0.9536]
    l, σ = modellength(track, reverse(1:20).*5, FCKetcham2007) .* 16.38
    @test l ≈ 14.558435101305054
    @test σ ≈  1.1716416735159647
    @test round.(track.r, sigdigits=4) ≈ [0.5945, 0.6823, 0.7347, 0.7723, 0.8014, 0.825, 0.8447, 0.8615, 0.8759, 0.8885, 0.8996, 0.9094, 0.9181, 0.926, 0.933, 0.9394, 0.9452, 0.9506, 0.956, 0.9619]

    l, σ =  modellength(track, reverse(1:20).*10, FCKetcham1999) .* 16.38 
    @test l ≈ 14.357128787814558
    @test σ ≈  1.1841227131850056
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6597, 0.7484, 0.8021, 0.8403, 0.8693, 0.8921, 0.9105, 0.9255, 0.9381, 0.9494]
    l, σ = modellength(track, reverse(1:20).*10, FCKetcham2007) .* 16.38
    @test l ≈ 14.538333058948192
    @test σ ≈  1.1899130108299063
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6564, 0.7592, 0.8162, 0.8549, 0.8834, 0.9053, 0.9225, 0.9364, 0.9479, 0.9581]

## -- Check the mean track length of a modelled Fish Canyon Apatite to be 15.35 +/- 0.06 um

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:28)), dpar=2.16)
    l, σ = modellength(track, fill(20., 28), FCKetcham2007) .* 16.38
    @test l ≈ 15.35 atol=0.06
    @test l ≈ 15.376461527029509 
    @test σ ≈  0.09967577151793756

    @test Thermochron.model_ll(track, fill(20., 28), FCKetcham2007) ≈ 1.0651088634382597

## --- Test c-axis equivalent model lengths

    @test lcmod(track) ≈ lcmod(15, 35) ≈ 15.405678663327869

    l = [13.9, 12.73, 13.98, 12.75, 13.77, 15.93, 13.99, 15.41, 13.23, 12.04, 13.01, 10.94]
    θ = [59.19, 78.45, 11.98, 68.73, 41.21, 35.05, 45.08, 41.42, 30.04, 53.9, 62.06, 70.21]
    @test lcmod.(l, θ) ≈ [14.980169300957634, 14.430422822748183, 14.085946809091633, 14.356417336926649, 14.570689505578352, 16.159195153315775, 14.812278052422833, 15.815567116293908, 13.865575929923889, 13.644277619942535, 14.435368624348033, 13.236076028425188] 

## --- 