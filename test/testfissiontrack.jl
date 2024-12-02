# Test fission track annealing

    @test FCKetcham1999 isa Thermochron.FanningCurvilinear{Float64}
    @test FCKetcham2007 isa Thermochron.SimplifiedCurvilinear{Float64}

## --- Ketcham et al. 1999
    am = FCKetcham1999
    @test am isa Thermochron.FanningCurvilinear{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9697793598061583
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.964698241959342
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.8757409060414316
    @test Thermochron.reltracklength(1, 200, am) ≈ 0.6865992067636948
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.0

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0017399827659966166, 0.003298712056416963, 0.004795641117812894, 0.0062537982811099796, 0.007683787046370734, 0.009091721758817036, 0.010481577045220772, 0.01185614029716061, 0.013217472138400533, 0.01456715633288802]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [660.2918270322567, 1394.7742722247667, 2160.135692291735, 2946.265864298046, 3748.209693768477, 4562.984978421849, 5388.588033345723, 6223.575180578797, 7066.854409034921, 7917.569525688446]
 
    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9516469756898532
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9435171871349471
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.8011854496662908
    @test Thermochron.reltrackdensity(1, 200, am) ≈ 0.3212180867210108
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

## --- C-axis equivalent model lengths

    l = [13.9, 12.73, 13.98, 12.75, 13.77, 15.93, 13.99, 15.41, 13.23, 12.04, 13.01, 10.94]
    θ = [59.19, 78.45, 11.98, 68.73, 41.21, 35.05, 45.08, 41.42, 30.04, 53.9, 62.06, 70.21]
    @test lcmodel.(l, θ) ≈ [14.980169300957634, 14.430422822748183, 14.085946809091633, 14.356417336926649, 14.570689505578352, 16.159195153315775, 14.812278052422833, 15.815567116293908, 13.865575929923889, 13.644277619942535, 14.435368624348033, 13.236076028425188] 

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

    # Isothermal residence
    @test modelage(apatite, fill(0, 100), FCKetcham1999) ≈ 89.48198019056863
    @test modelage(apatite, fill(0, 100), FCKetcham2007) ≈ 91.18888272552469

    @test modelage(apatite, fill(50, 100), FCKetcham1999) ≈ 71.77101682177967
    @test modelage(apatite, fill(50, 100), FCKetcham2007) ≈ 74.20874464530966

    @test modelage(apatite, fill(75, 100), FCKetcham1999) ≈ 22.14984699770526
    @test modelage(apatite, fill(75, 100), FCKetcham2007) ≈ 21.6047849747377

    @test modelage(apatite, fill(100, 100), FCKetcham1999) ≈ 0.4704161148027164
    @test modelage(apatite, fill(100, 100), FCKetcham2007) ≈ 0.42338708872671615

    # Linear cooling
    @test modelage(apatite, reverse(1:100), FCKetcham1999) ≈ 66.77980999111053
    @test modelage(apatite, reverse(1:100), FCKetcham2007) ≈ 68.4228509734102

    apatite = ApatiteFT(agesteps=reverse(cntr(0:200)), F=1.75, Cl=0.01, OH=0.24)
    @test apatite isa ApatiteFT{Float64}
    @test apatite.rmr0 ≈ 0.8573573076438294

    @test modelage(apatite, reverse(1:200), FCKetcham1999) ≈ 65.22193612822629
    @test modelage(apatite, reverse(1:200), FCKetcham2007) ≈ 67.08221305145062

    @test modelage(apatite, reverse(1:200)./2, FCKetcham1999) ≈ 125.34528518852127
    @test modelage(apatite, reverse(1:200)./2, FCKetcham2007) ≈ 128.59190431998582

    @test modelage(apatite, reverse(1:200).*2, FCKetcham1999) ≈ 33.39910745952449
    @test modelage(apatite, reverse(1:200).*2, FCKetcham2007) ≈ 34.52279548202092

## --- Test track lengths

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:20)), F=1.75, Cl=0.01, OH=0.24)
    l, σ = modellength(track, fill(75, 20), FCKetcham1999) .* 16.38
    @test l ≈ 12.434055894969841
    @test σ ≈  0.5731491659110551 
    @test round.(track.r, sigdigits=4) ≈ [0.7103, 0.7135, 0.7169, 0.7203, 0.7239, 0.7276, 0.7315, 0.7355, 0.7398, 0.7442, 0.749, 0.754, 0.7595, 0.7654, 0.7719, 0.7793, 0.7878, 0.798, 0.8112, 0.8311]
    l, σ = modellength(track, fill(75, 20), FCKetcham2007) .* 16.38
    @test l ≈ 12.610086362166317
    @test σ ≈  0.6062929557510826
    @test round.(track.r, sigdigits=4) ≈ [0.7176, 0.7212, 0.725, 0.7289, 0.7329, 0.737, 0.7413, 0.7457, 0.7503, 0.7552, 0.7603, 0.7657, 0.7716, 0.7779, 0.7848, 0.7925, 0.8014, 0.8119, 0.8254, 0.8454]

    l, σ = modellength(track, fill(100, 20), FCKetcham1999) .* 16.38
    @test l ≈ 10.875830406686777
    @test σ ≈  0.6360569439228593
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4559, 0.5506, 0.6187, 0.6878]
    l, σ = modellength(track, fill(100, 20), FCKetcham2007) .* 16.38
    @test l ≈ 10.982218818618923
    @test σ ≈  0.608187681868191
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4888, 0.6017, 0.6895]

    l, σ = modellength(track, reverse(1:20).*5, FCKetcham1999) .* 16.38
    @test l ≈ 14.353530832531991
    @test σ ≈  1.1780179817828795
    @test round.(track.r, sigdigits=4) ≈ [0.6251, 0.6876, 0.7304, 0.7632, 0.7897, 0.8118, 0.8307, 0.8472, 0.8616, 0.8744, 0.8858, 0.896, 0.9052, 0.9135, 0.9211, 0.928, 0.9345, 0.9406, 0.9467, 0.9536]
    l, σ = modellength(track, reverse(1:20).*5, FCKetcham2007) .* 16.38
    @test l ≈ 14.548247656072547
    @test σ ≈  1.1726032551206245
    @test round.(track.r, sigdigits=4) ≈ [0.6095, 0.6893, 0.7388, 0.7749, 0.8031, 0.8261, 0.8454, 0.8618, 0.8761, 0.8885, 0.8994, 0.9091, 0.9178, 0.9256, 0.9326, 0.9389, 0.9448, 0.9503, 0.9557, 0.9619]

    l, σ =  modellength(track, reverse(1:20).*10, FCKetcham1999) .* 16.38 
    @test l ≈ 14.346205934457554
    @test σ ≈  1.1850553007509892
    @test round.(track.r, sigdigits=4) ≈[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6595, 0.7476, 0.8012, 0.8392, 0.8682, 0.8911, 0.9095, 0.9247, 0.9375, 0.9494]
    l, σ = modellength(track, reverse(1:20).*10, FCKetcham2007) .* 16.38
    @test l ≈ 14.529573753525176
    @test σ ≈  1.190954235529948
    @test round.(track.r, sigdigits=4) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6558, 0.7584, 0.8153, 0.854, 0.8825, 0.9044, 0.9218, 0.9358, 0.9475, 0.9581]

## -- Check the mean track length of a modelled Fish Canyon Apatite to be 15.35 +/- 0.06 um

    track = ApatiteTrackLength(length=15, angle=35, agesteps=reverse(cntr(0:28)), dpar=2.16)
    l, σ = modellength(track, fill(20., 28), FCKetcham2007) .* 16.38
    @test l ≈ 15.35 atol=0.06
    @test l ≈ 15.376461527029509 
    @test σ ≈  0.09967577151793756
    
## --- 