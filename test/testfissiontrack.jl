# Test fission track annealing

    @test Ketcham1999FC() isa Thermochron.ApatiteAnnealingModel{Float64}
    @test Ketcham2007FC() isa Thermochron.ApatiteAnnealingModel{Float64}
    @test Yamada2007PC() isa Thermochron.ZirconAnnealingModel{Float64}
    @test Guenthner2013FC() isa Thermochron.ZirconAnnealingModel{Float64}

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

## --- Jones et al. 2021 Parallel Linear monazite

    am = Jones2021FA()
    @test am isa Thermochron.Jones2021FA{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.7579417895276482
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.735387946933017
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.532403363581336
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.30686493763502365
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.0

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0007407975532689652, 0.0013501847848567561, 0.0019181794368659146, 0.0024608598465454146, 0.0029854616491727447, 0.00349609239238555, 0.003995390180771109, 0.004485186955341078, 0.00496682573501881, 0.00544133127424477]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [1349.8964671079773, 2962.55745499635, 4691.948952755405, 6501.7924618751285, 8373.91430130319, 10297.210702556848, 12264.133860023401, 14269.193377499583, 16308.20252639552, 18377.855521005877]

    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.5158835790552965
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.470775893866034
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.064806727162672
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.0
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Yamada et al. 2007 Parallel Curvilinear zircon

    am = Yamada2007PC()
    @test am isa Thermochron.Yamada2007PC{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9964516173825159
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.995078688092618
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.940753376497744
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.587467392379713
    @test Thermochron.reltracklength(1, 300, am) ≈ 0.047145496695413736
    @test Thermochron.reltracklength(1, 500, am) ≈ 4.8685772129709805e-21

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0020577410195178228, 0.0041154820390356525, 0.006173223058553454, 0.00823096407807132, 0.010288705097588996, 0.01234644611710693, 0.014404187136625298, 0.016461928156142672, 0.018519669175660315, 0.02057741019517803]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [222.98141544261495, 445.9628308852307, 668.9442463278432, 891.9256617704631, 1114.907077213062, 1337.8884926556889, 1560.8699080982742, 1783.8513235409293, 2006.8327389835247, 2229.814154426128] 

    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9955645217281448
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9938483601157724
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.92594172062218 
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.48433424047464124
    @test Thermochron.reltrackdensity(1, 300, am) ≈ 0.0
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Guenther 2013 Simplified Curvilinear zircon

    am = Guenthner2013FC()
    @test am isa Thermochron.Guenthner2013FC{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9987988651472297
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.998476438713864
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.9886125414156559
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.9156810748820664
    @test Thermochron.reltracklength(1, 300, am) ≈ 0.5980310995551887
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.029152068246667134

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0028785665193461073, 0.005690055255432659, 0.008476770230793686, 0.011247518024086979, 0.014006452961030012, 0.016756010561736853, 0.01949779438472927, 0.022232940810438962, 0.02496229740141236, 0.02768652061547207]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [196.56951935291298, 397.32010926193027, 599.6796926328979, 803.0912917913217, 1007.2887605410821, 1212.114684837377, 1417.464714169357, 1623.2644860315188, 1829.4582748135415, 2036.0027420018664]

    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9984985814340372
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9980955483923302
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.9857656767695699
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.8946013436025829
    @test Thermochron.reltrackdensity(1, 300, am) ≈ 0.4975388744439859
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
    println()
    display(zircon)

    # Isothermal residence
    @test modelage(zircon, fill(0, 100), Yamada2007PC()) ≈ 99.03307472647863
    @test modelage(zircon, fill(50, 100), Yamada2007PC()) ≈ 95.5875134964696
    @test modelage(zircon, fill(75, 100), Yamada2007PC()) ≈ 91.44848076970709
    @test modelage(zircon, fill(100, 100), Yamada2007PC()) ≈ 84.4038173798344
    @test modelage(zircon, fill(0, 100), Guenthner2013FC()) ≈ 99.72322073315222
    @test modelage(zircon, fill(50, 100), Guenthner2013FC()) ≈ 99.06905675284062
    @test modelage(zircon, fill(75, 100), Guenthner2013FC()) ≈ 98.3387065429949
    @test modelage(zircon, fill(100, 100), Guenthner2013FC()) ≈ 97.08526701487509

    # Linear cooling
    @test modelage(zircon, reverse(1:100), Yamada2007PC()) ≈ 95.86292049620485
    @test modelage(zircon, reverse(1:100), Guenthner2013FC()) ≈ 99.1377750755283

    # As above but longer history
    zircon = ZirconFT(agesteps=reverse(cntr(0:200)))
    @test zircon isa ZirconFT{Float64}

    @test modelage(zircon, reverse(1:200), Yamada2007PC()) ≈ 158.88973345523144
    @test modelage(zircon, reverse(1:200)./2, Yamada2007PC()) ≈ 190.62094301177618 
    @test modelage(zircon, reverse(1:200).*2, Yamada2007PC()) ≈ 84.30136045990102
    @test modelage(zircon, reverse(1:200), Guenthner2013FC()) ≈ 191.5897369291866
    @test modelage(zircon, reverse(1:200)./2, Guenthner2013FC()) ≈ 198.06912276144533
    @test modelage(zircon, reverse(1:200).*2, Guenthner2013FC()) ≈ 132.1966344789464 

    # Fish Canyon Tuff zircon example
    zircon = ZirconFT(age=27, age_sigma=3, agesteps=reverse(cntr(0:28)))
    @test modelage(zircon, fill(20., 28), Yamada2007PC()) ≈ 27.601238386594986
    @test Thermochron.model_ll(zircon, fill(20., 28), Yamada2007PC()) ≈ -2.03763346617919
    @test modelage(zircon, fill(20., 28), Guenthner2013FC()) ≈ 27.89690820436307
    @test Thermochron.model_ll(zircon, fill(20., 28), Guenthner2013FC()) ≈ -2.06224217337577

## --- Test apatite fission track model ages

    apatite = ApatiteFT(agesteps=reverse(cntr(0:100)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294
    show(apatite)
    println()
    display(apatite)

    # Isothermal residence
    @test modelage(apatite, fill(0, 100), Ketcham1999FC()) ≈ 89.54730962518295
    @test modelage(apatite, fill(50, 100), Ketcham1999FC()) ≈ 71.90164016667141
    @test modelage(apatite, fill(75, 100), Ketcham1999FC()) ≈ 22.135848309417995
    @test modelage(apatite, fill(100, 100), Ketcham1999FC()) ≈ 0.46834257592280504
    @test modelage(apatite, fill(0, 100), Ketcham2007FC()) ≈ 91.24703624063753
    @test modelage(apatite, fill(50, 100), Ketcham2007FC()) ≈ 74.34211338104792
    @test modelage(apatite, fill(75, 100), Ketcham2007FC()) ≈ 21.641396305809224
    @test modelage(apatite, fill(100, 100), Ketcham2007FC()) ≈ 0.42342028286401207

    # Linear cooling
    @test modelage(apatite, reverse(1:100), Ketcham1999FC()) ≈ 66.2164833900794
    @test modelage(apatite, reverse(1:100), Ketcham2007FC()) ≈ 67.93014962856782

    # As above but longer history
    apatite = ApatiteFT(agesteps=reverse(cntr(0:200)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294

    @test modelage(apatite, reverse(1:200), Ketcham1999FC()) ≈ 66.2164833900794
    @test modelage(apatite, reverse(1:200)./2, Ketcham1999FC()) ≈ 124.45929282749134
    @test modelage(apatite, reverse(1:200).*2, Ketcham1999FC()) ≈ 34.938406635567084
    @test modelage(apatite, reverse(1:200), Ketcham2007FC()) ≈ 67.93014962856782
    @test modelage(apatite, reverse(1:200)./2, Ketcham2007FC()) ≈ 127.81186983354935
    @test modelage(apatite, reverse(1:200).*2, Ketcham2007FC()) ≈ 35.810590270883125

    # Fish Canyon Tuff apatite example
    apatite = ApatiteFT(age=25, age_sigma=3, agesteps=reverse(cntr(0:28)), dpar=2.16)
    @test modelage(apatite, fill(20., 28), Ketcham2007FC()) ≈ 25.25753183867445
    @test Thermochron.model_ll(apatite, fill(20., 28), Ketcham2007FC()) ≈ -2.021235413424507 

## -- Test monazite fission track

    monazite = MonaziteFT(agesteps=reverse(cntr(0:100)))
    @test monazite isa MonaziteFT{Float64}
    show(monazite)
    println()
    display(monazite)

    # Isothermal residence
    @test modelage(monazite, fill(0, 100), Jones2021FA()) ≈ 43.43881113607752
    @test modelage(monazite, fill(50, 100), Jones2021FA()) ≈ 19.288380906685823
    @test modelage(monazite, fill(75, 100), Jones2021FA()) ≈ 7.179152374521148
    @test modelage(monazite, fill(100, 100), Jones2021FA()) ≈ 0.18694393293501763

    # Linear cooling
    @test modelage(monazite, reverse(1:100), Jones2021FA()) ≈ 24.054872653095966

    monazite = MonaziteFT(age=10, age_sigma=3, agesteps=reverse(cntr(0:28)))
    @test modelage(monazite, fill(20., 28), Jones2021FA()) ≈ 10.263270122494962
    @test Thermochron.model_ll(monazite, fill(20., 28), Jones2021FA()) ≈ -2.0214014417282553

## --- Test track lengths

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:20)), F=1.75, Cl=0.01, OH=0.24)
    show(track)
    println()
    display(track)

    l, σ = modellength(track, fill(75, 20), Ketcham1999FC())
    @test l ≈ 12.431835901971635
    @test σ ≈  0.5806318171054747
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.7101, 0.7134, 0.7167, 0.7202, 0.7238, 0.7275, 0.7314, 0.7354, 0.7396, 0.7441, 0.7488, 0.7539, 0.7593, 0.7653, 0.7718, 0.7792, 0.7877, 0.7979, 0.8111, 0.8311]
    l, σ = modellength(track, fill(75, 20), Ketcham2007FC())
    @test l ≈ 12.609580169379788
    @test σ ≈  0.6128738590327928
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.7176, 0.7212, 0.725, 0.7289, 0.7329, 0.737, 0.7413, 0.7457, 0.7503, 0.7552, 0.7603, 0.7657, 0.7716, 0.7779, 0.7848, 0.7925, 0.8014, 0.8119, 0.8254, 0.8454]

    l, σ = modellength(track, fill(100, 20), Ketcham1999FC())
    @test l ≈ 10.874448549717606
    @test σ ≈  0.6408060359179175
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4544, 0.55, 0.6184, 0.6876]
    l, σ = modellength(track, fill(100, 20), Ketcham2007FC())
    @test l ≈ 10.98218090689383
    @test σ ≈  0.6148371520144172
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4888, 0.6017, 0.6895]

    l, σ = modellength(track, reverse(1:20).*5, Ketcham1999FC())
    @test l ≈ 14.362772443772794
    @test σ ≈  1.1860436444177835
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.6124, 0.6806, 0.726, 0.7602, 0.7876, 0.8104, 0.8298, 0.8466, 0.8614, 0.8744, 0.886, 0.8963, 0.9056, 0.914, 0.9216, 0.9286, 0.935, 0.941, 0.9469, 0.9536]
    l, σ = modellength(track, reverse(1:20).*5, Ketcham2007FC())
    @test l ≈ 14.55756747909137
    @test σ ≈  1.1754875315974673
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.5945, 0.6823, 0.7347, 0.7723, 0.8014, 0.825, 0.8447, 0.8615, 0.8759, 0.8885, 0.8996, 0.9094, 0.9181, 0.926, 0.933, 0.9394, 0.9452, 0.9506, 0.956, 0.9619] 

    l, σ =  modellength(track, reverse(1:20).*10, Ketcham1999FC()) 
    @test l ≈ 14.356679949551062
    @test σ ≈  1.187729515554466 
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6597, 0.7484, 0.8021, 0.8403, 0.8693, 0.8921, 0.9105, 0.9255, 0.9381, 0.9494]
    l, σ = modellength(track, reverse(1:20).*10, Ketcham2007FC())
    @test l ≈ 14.537888339640748
    @test σ ≈  1.193516585511275
    @test round.(track.r ./ 16.38, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6564, 0.7592, 0.8162, 0.8549, 0.8834, 0.9053, 0.9225, 0.9364, 0.9479, 0.9581]

## -- Check the mean track length of a modelled Fish Canyon Apatite to be 15.35 +/- 0.06 um

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:28)), dpar=2.16)
    l, σ = modellength(track, fill(20., 28), Ketcham2007FC(); trackhist=true)
    @test l ≈ 15.35 atol=0.06
    @test l ≈ 15.376345688670545
    @test σ ≈ 0.1342688478623522

    @test round.(track.ldist, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2570000000000001e-222, 7.635e-132, 4.523e-65, 3.2259999999999996e-22, 0.007669, 0.9922, 0.0001799, 0.0, 0.0, 0.0]
    @test Thermochron.model_ll(track, fill(20., 28), Ketcham2007FC()) ≈ 0.840405624475848

## --- Test c-axis equivalent model lengths

    @test lcmod(track) ≈ lcmod(15, 35) ≈ 15.405678663327869

    l = [13.9, 12.73, 13.98, 12.75, 13.77, 15.93, 13.99, 15.41, 13.23, 12.04, 13.01, 10.94]
    θ = [59.19, 78.45, 11.98, 68.73, 41.21, 35.05, 45.08, 41.42, 30.04, 53.9, 62.06, 70.21]
    @test lcmod.(l, θ) ≈ [14.980169300957634, 14.430422822748183, 14.085946809091633, 14.356417336926649, 14.570689505578352, 16.159195153315775, 14.812278052422833, 15.815567116293908, 13.865575929923889, 13.644277619942535, 14.435368624348033, 13.236076028425188] 

## --- 