# Test fission track annealing

    @test FCKetcham1999 isa Thermochron.FanningCurvilinear{Float64}
    @test FCKetcham2007 isa Thermochron.FanningCurvilinear{Float64}

## --- Ketcham et al. 1999
    am = FCKetcham1999
    @test am isa Thermochron.FanningCurvilinear{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.9697793598061583
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.964698241959342
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.8757409060414316
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.0

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0017399827659966166, 0.003298712056416963, 0.004795641117812894, 0.0062537982811099796, 0.007683787046370734, 0.009091721758817036, 0.010481577045220772, 0.01185614029716061, 0.013217472138400533, 0.01456715633288802]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [660.2918270322567, 1394.7742722247667, 2160.135692291735, 2946.265864298046, 3748.209693768477, 4562.984978421849, 5388.588033345723, 6223.575180578797, 7066.854409034921, 7917.569525688446]
 
    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.9516469756898532
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.9435171871349471
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.8011854496662908
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- Ketcham et al. 2007

    am = FCKetcham2007
    @test am isa Thermochron.FanningCurvilinear{Float64}

    @test Thermochron.reltracklength(1, 0, am) ≈ 0.3043691833790162
    @test Thermochron.reltracklength(1, 10, am) ≈ 0.30292827763276786
    @test Thermochron.reltracklength(1, 100, am) ≈ 0.290329105174936
    @test Thermochron.reltracklength(1, 500, am) ≈ 0.23553543695185192

    @test Thermochron.equivalenttime.(1:10, 100, 100, am) ≈ 1:10
    @test Thermochron.equivalenttime.(1:10, 50, 100, am) ≈ [0.0015478065155001086, 0.002954611923287363, 0.00431269182539665, 0.005640066461673555, 0.006945066485586564, 0.008232508754276639, 0.009505531052511628, 0.010766337684274287, 0.012016560593707997, 0.013257455657955832]
    @test Thermochron.equivalenttime.(1:10, 150, 100, am) ≈ [646.5878531291968, 1354.896862856624, 2088.5372542290797, 2839.1277381635086, 3602.579481389471, 4376.439427402521, 5159.063447504911, 5949.269302029918, 6746.164254698982, 7549.049388886254]
    
    @test Thermochron.reltrackdensity(1, 0, am) ≈ 0.0
    @test Thermochron.reltrackdensity(1, 10, am) ≈ 0.0
    @test Thermochron.reltrackdensity(1, 100, am) ≈ 0.0
    @test Thermochron.reltrackdensity(1, 500, am) ≈ 0.0

## --- 