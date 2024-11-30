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

## --- 