# Test geometry utility functions

    # Intersection of two spheres
    @test size(intersectionfraction.(1.,1.,0:0.2:2)) == (11,)
    @test intersectionfraction.(1.,1.,0:0.2:2) ≈ [1//1, 9//20, 2//5, 7//20, 3//10, 1//4, 1//5, 3//20, 1//10, 1//20, 0//1]

    # Intersection of spherical shells
    crystalRadius = 41.
    dr = 1.
    rEdges = collect(0. : dr : crystalRadius) # Edges of each radius element
    relVolumes = (rEdges[2:end].^3 - rEdges[1:end-1].^3)/rEdges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element
    ralpha = 16.69

    # Zero when outside range
    @test intersectiondensity(rEdges,relVolumes,ralpha,60.) == zeros(length(relVolumes))
    # Single shell when centered
    dInt = intersectiondensity(rEdges,relVolumes,ralpha,0.)
    @test dInt[round(Int,ralpha)] > 0
    @test all(dInt[(1:length(relVolumes)) .!= round(Int,ralpha)] .== 0)
        # Specific value tests
    @test round.(intersectiondensity(rEdges,relVolumes,ralpha,0.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.3586, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]
    @test round.(intersectiondensity(rEdges,relVolumes,ralpha,5.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.81946, 11.0061, 10.1916, 9.48929, 8.87752, 8.33983, 7.86353, 7.43869, 7.05739, 6.71327, 4.38494, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]
    @test round.(intersectiondensity(rEdges,relVolumes,ralpha,20.), sigdigits=6) ≈ [0, 0, 0, 7.03672, 7.61584, 6.2396, 5.28378, 4.58152, 4.04384, 3.61901, 3.27489, 2.99049, 2.75152, 2.5479, 2.37232, 2.21938, 2.08496, 1.96588, 1.85967, 1.76435, 1.67832, 1.60028, 1.52918, 1.46413, 1.40439, 1.34933, 1.29842, 1.25122, 1.20733, 1.16641, 1.12817, 1.09236, 1.05876, 1.02716, 0.997389, 0.969297, 0.647731, 0, 0, 0, 0.0]
    @test round.(intersectiondensity(rEdges,relVolumes,ralpha,40.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.508458, 0.702195, 0.674665, 0.649212, 0.62561, 0.603663, 0.583204, 0.564086, 0.546181, 0.529379, 0.513579, 0.498694, 0.484649, 0.471372, 0.458804, 0.446888, 0.435576, 0.424822]


# Test damage annealing function
    tCryst = 3000.0 # Time (in AMyr)
    dt = 100 # time step size in Myr
    dr = 1
    tSteps = collect(0+dt/2 : dt : tCryst-dt/2)

    TSteps = fill(250.0, length(tSteps))
    pr_known = [0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0 0.0; 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0 0.0; 0.349337 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0 0.0; 0.347022 0.349337 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0 0.0; 0.344796 0.347022 0.349337 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0 0.0; 0.342653 0.344796 0.347022 0.349337 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002 0.0; 0.340589 0.342653 0.344796 0.347022 0.349337 0.351749 0.354267 0.356898 0.359654 0.362546 0.365588 0.368795 0.372185 0.375779 0.379603 0.383687 0.388066 0.392785 0.397898 0.403474 0.409603 0.4164 0.424023 0.432692 0.44273 0.454631 0.469222 0.488036 0.514459 0.559002]

    pr, Teq = anneal(dt, tSteps, TSteps, ZRDAAM())
    @test isa(pr, AbstractMatrix)
    @test size(pr) == (30,30)
    @test round.(pr, sigdigits=6) ≈ pr_known
    @test isa(Teq, AbstractVector)
    @test length(Teq) == 30

    for i=1:10
        anneal!(pr, Teq, dt, tSteps, TSteps)
        @test round.(pr, sigdigits=6) ≈ pr_known
    end

    TSteps = collect(range(650, 0, length=length(tSteps)))
    pr_known = [0.000205838 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.00018167 0.000345328 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000175903 0.000306842 0.000575964 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000174306 0.00029782 0.000515098 0.000955418 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173858 0.000295383 0.000501104 0.000859803 0.00157683 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173735 0.000294718 0.00049742 0.000838273 0.00142761 0.00259 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173702 0.000294541 0.000496445 0.000832761 0.00139476 0.00235863 0.00423469 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173694 0.000294495 0.000496194 0.000831348 0.00138659 0.00230891 0.00387837 0.00689242 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173692 0.000294484 0.000496131 0.000830996 0.00138456 0.00229693 0.00380376 0.00634758 0.0111652 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496117 0.000830911 0.00138408 0.00229407 0.00378637 0.00623665 0.0103389 0.0179917 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496113 0.000830892 0.00138397 0.00229341 0.00378237 0.00621171 0.0101757 0.0167508 0.0288076 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830888 0.00138394 0.00229326 0.00378149 0.0062062 0.0101404 0.0165134 0.0269675 0.0457439 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229323 0.0037813 0.00620503 0.0101329 0.0164642 0.0266275 0.0430615 0.0718099 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378126 0.0062048 0.0101314 0.0164542 0.02656 0.0425841 0.0679927 0.110913 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378126 0.00620476 0.0101311 0.0164523 0.026547 0.0424937 0.0673402 0.105664 0.167408 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265447 0.0424772 0.0672227 0.104805 0.160534 0.244771 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265443 0.0424744 0.0672025 0.104658 0.15946 0.236352 0.34328 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671992 0.104635 0.159287 0.235101 0.333839 0.457702 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159261 0.234913 0.332512 0.448193 0.577183 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234886 0.332326 0.446936 0.568685 0.688768 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332302 0.446773 0.567636 0.682035 0.782538 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446753 0.56751 0.681264 0.777758 0.85451 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567497 0.681181 0.777256 0.851417 0.905927 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567496 0.681172 0.777207 0.851122 0.904069 0.940745 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851097 0.903911 0.939694 0.96343 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903899 0.939615 0.962863 0.977805 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962826 0.977511 0.986726 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977495 0.98658 0.992174 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977494 0.986574 0.992104 0.995453 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977494 0.986573 0.992102 0.995421 0.9974]

    pr, Teq = anneal(dt, tSteps, TSteps, ZRDAAM())

    @test isa(pr, AbstractMatrix)
    @test size(pr) == (30,30)
    @test round.(pr, sigdigits=6) ≈ pr_known
    @test isa(Teq, AbstractVector)
    @test length(Teq) == 30

    for i=1:10
        anneal!(pr, Teq, dt, tSteps, TSteps)
        @test round.(pr, sigdigits=6) ≈ pr_known
    end


# Test basic linear algebra
    dl = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    d = [1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 1.0]
    du = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    A = Tridiagonal(dl, d, du)
    y = [0.0; -0.900592; -2.72489; -4.54554; -6.36539; -8.18497; -10.0045; -11.8128; -13.5304; -15.2259; -16.8875; -18.5449; -20.1978; -21.804; -23.3226; -24.8305; -26.3291; -27.762; -29.1241; -30.2665; -31.0612; -31.7774; -32.1437; -32.1666; -31.9801; -31.5697; -31.0949; -30.565; -29.9796; -29.339; -28.6424; -27.8906; -27.0836; -26.2213; -25.3036; -24.3306; 0.0]

    sol_known = [-159.97020774647876, 159.97020774647876, 479.0100312394363, 795.3249647323938, 1107.0943582253512, 1412.4983617183084, 1709.7173952112655, 1996.9319287042226, 2272.33366219718, 2534.2049956901374, 2780.8504291830945, 3010.6083626760515, 3221.8213961690085, 3412.836629661965, 3582.0478631549227, 3727.93649664788, 3848.9946301408377, 3943.7236636337952, 4010.6906971267535, 4048.533630619712, 4056.110064112671, 4032.62529760563, 3977.3631310985884, 3889.9572645915464, 3770.3847980845044, 3618.832231577462, 3435.70996507042, 3221.492798563378, 2976.710632056336, 2701.948865549294, 2397.8480990422527, 2065.104932535211, 1704.4711660281687, 1316.7537995211264, 902.8151330140844, 463.5728665070422, 0.0]

    @test isapprox(A\y, sol_known, atol=0.00000000001)

# Test creating and allocating a Zircon

    crystalRadius = 29.26
    Uppm = 462.98
    Thppm = 177.76
    Zircon(crystalRadius,dr,Uppm,Thppm,dt,reverse(tSteps))
    @time "Allocating a zircon" zircon = Zircon(crystalRadius,dr,Uppm,Thppm,dt,reverse(tSteps))
    @test isa(zircon, Zircon)
    @test zircon.ageSteps == reverse(tSteps)
    @test zircon.r238U ≈ fill(1.1714561176470587e18, 29)
    @test zircon.r235U ≈ fill(8.608533548562195e15, 29)
    @test zircon.r232Th ≈ fill(4.614097931034483e17, 29)
    @test zircon.rEdges == 0:dr:crystalRadius
    @test zircon.rSteps == (zircon.rEdges[2:end] + zircon.rEdges[1:end-1])/2
    @test zircon.nrSteps == length(zircon.rSteps) + 2 # Implicit radius steps inside and outside modeled range
    @test round.(zircon.relVolumes, sigdigits=5) ≈ [4.1002e-5, 0.00028701, 0.00077904, 0.0015171, 0.0025011, 0.0037312, 0.0052073, 0.0069294, 0.0088975, 0.011112, 0.013572, 0.016278, 0.01923, 0.022428, 0.025872, 0.029563, 0.033499, 0.037681, 0.042109, 0.046783, 0.051704, 0.05687, 0.062282, 0.06794, 0.073845, 0.079995, 0.086391, 0.093034, 0.099922]
    @test round.(zircon.r238UHe, sigdigits=5) ≈ [9.0969e18, 8.9297e18, 8.862e18, 8.825e18, 8.7994e18, 8.7792e18, 8.7621e18, 8.7468e18, 8.7327e18, 8.7196e18, 8.7071e18, 8.6955e18, 8.6013e18, 8.367e18, 8.1542e18, 7.8435e18, 7.4452e18, 7.0531e18, 6.6323e18, 6.2396e18, 5.8706e18, 5.522e18, 5.1911e18, 4.8756e18, 4.5735e18, 4.2834e18, 4.0039e18, 3.7338e18, 3.4722e18]
    @test round.(zircon.r235UHe, sigdigits=5) ≈ [6.0636e16, 6.0327e16, 5.9517e16, 5.8385e16, 5.7712e16, 5.7253e16, 5.6412e16, 5.4727e16, 5.3357e16, 5.2219e16, 5.0705e16, 4.8967e16, 4.7197e16, 4.5608e16, 4.3853e16, 4.2003e16, 4.0247e16, 3.809e16, 3.6078e16, 3.4188e16, 3.24e16, 3.0705e16, 2.909e16, 2.7541e16, 2.6051e16, 2.4615e16, 2.3226e16, 2.1879e16, 2.0569e16]
    @test round.(zircon.r232ThHe, sigdigits=5) ≈ [2.5039e18, 2.5086e18, 2.5127e18, 2.5129e18, 2.5114e18, 2.5029e18, 2.4566e18, 2.4204e18, 2.3906e18, 2.3652e18, 2.3432e18, 2.3013e18, 2.2224e18, 2.1221e18, 2.0316e18, 1.9492e18, 1.8734e18, 1.8028e18, 1.7199e18, 1.6261e18, 1.5376e18, 1.4537e18, 1.3739e18, 1.2975e18, 1.224e18, 1.1533e18, 1.085e18, 1.0187e18, 9.544e17]

    alphaDeposition_known = [3.4645e17 3.4182e17 3.3873e17 3.3578e17 3.3394e17 3.3257e17 3.3037e17 3.2675e17 3.2377e17 3.2126e17 3.181e17 3.1445e17 3.085e17 2.9933e17 2.9043e17 2.7902e17 2.6566e17 2.5176e17 2.3735e17 2.2378e17 2.1102e17 1.9894e17 1.8747e17 1.765e17 1.66e17 1.559e17 1.4616e17 1.3673e17 1.2759e17; 3.3271e17 3.282e17 3.2526e17 3.2252e17 3.208e17 3.1952e17 3.1747e17 3.1415e17 3.114e17 3.0908e17 3.0619e17 3.0283e17 2.9723e17 2.8841e17 2.799e17 2.6891e17 2.56e17 2.4262e17 2.2871e17 2.1561e17 2.0329e17 1.9163e17 1.8056e17 1.6998e17 1.5984e17 1.5009e17 1.4069e17 1.316e17 1.2278e17; 3.1998e17 3.1558e17 3.128e17 3.1024e17 3.0863e17 3.0742e17 3.0551e17 3.0246e17 2.9993e17 2.9779e17 2.9513e17 2.9205e17 2.8675e17 2.7826e17 2.701e17 2.5951e17 2.4702e17 2.3412e17 2.2068e17 2.0802e17 1.9611e17 1.8484e17 1.7414e17 1.6392e17 1.5412e17 1.447e17 1.3562e17 1.2684e17 1.1832e17; 3.0818e17 3.0388e17 3.0123e17 2.9885e17 2.9734e17 2.9621e17 2.9442e17 2.916e17 2.8927e17 2.8729e17 2.8485e17 2.8201e17 2.7699e17 2.6882e17 2.6098e17 2.5077e17 2.3867e17 2.2621e17 2.132e17 2.0095e17 1.8942e17 1.7852e17 1.6817e17 1.5828e17 1.488e17 1.3969e17 1.3091e17 1.2241e17 1.1417e17; 2.9721e17 2.9301e17 2.9049e17 2.8827e17 2.8685e17 2.8578e17 2.8411e17 2.8151e17 2.7935e17 2.7752e17 2.7528e17 2.7266e17 2.679e17 2.6002e17 2.5248e17 2.4261e17 2.3088e17 2.1883e17 2.0623e17 1.9436e17 1.832e17 1.7264e17 1.6261e17 1.5303e17 1.4385e17 1.3502e17 1.2652e17 1.1829e17 1.1031e17; 2.8701e17 2.8291e17 2.805e17 2.7842e17 2.7709e17 2.7608e17 2.7451e17 2.7211e17 2.7011e17 2.6842e17 2.6635e17 2.6394e17 2.5942e17 2.518e17 2.4454e17 2.3499e17 2.236e17 2.1194e17 1.9972e17 1.8821e17 1.7738e17 1.6714e17 1.5741e17 1.4813e17 1.3923e17 1.3067e17 1.2242e17 1.1444e17 1.0671e17; 2.7751e17 2.735e17 2.712e17 2.6924e17 2.6799e17 2.6704e17 2.6556e17 2.6334e17 2.6149e17 2.5992e17 2.5801e17 2.5579e17 2.5148e17 2.4411e17 2.3711e17 2.2787e17 2.168e17 2.055e17 1.9364e17 1.8246e17 1.7195e17 1.6201e17 1.5256e17 1.4355e17 1.3491e17 1.266e17 1.186e17 1.1085e17 1.0335e17; 2.6865e17 2.6472e17 2.6252e17 2.6068e17 2.595e17 2.586e17 2.5721e17 2.5515e17 2.5344e17 2.5198e17 2.5022e17 2.4816e17 2.4405e17 2.3691e17 2.3015e17 2.2119e17 2.1042e17 1.9946e17 1.8794e17 1.7707e17 1.6686e17 1.572e17 1.4802e17 1.3926e17 1.3087e17 1.228e17 1.1502e17 1.075e17 1.002e17; 2.6036e17 2.5651e17 2.5441e17 2.5268e17 2.5157e17 2.5071e17 2.4939e17 2.4748e17 2.4589e17 2.4454e17 2.4291e17 2.41e17 2.3708e17 2.3016e17 2.2362e17 2.1492e17 2.0444e17 1.938e17 1.8259e17 1.7202e17 1.6208e17 1.5269e17 1.4376e17 1.3524e17 1.2708e17 1.1923e17 1.1167e17 1.0435e17 9.7261e16; 2.5261e17 2.4884e17 2.4682e17 2.4519e17 2.4413e17 2.4332e17 2.4207e17 2.403e17 2.3882e17 2.3756e17 2.3605e17 2.3428e17 2.3053e17 2.2381e17 2.1748e17 2.0903e17 1.9882e17 1.8848e17 1.7756e17 1.6727e17 1.576e17 1.4845e17 1.3976e17 1.3146e17 1.2352e17 1.1588e17 1.0852e17 1.014e17 9.4498e16; 2.4533e17 2.4164e17 2.397e17 2.3816e17 2.3716e17 2.3638e17 2.3521e17 2.3356e17 2.3218e17 2.31e17 2.296e17 2.2796e17 2.2436e17 2.1783e17 2.1169e17 2.0348e17 1.9353e17 1.8346e17 1.7283e17 1.628e17 1.5337e17 1.4446e17 1.3599e17 1.2791e17 1.2017e17 1.1273e17 1.0556e17 9.8624e16 9.1902e16; 2.385e17 2.3488e17 2.3301e17 2.3156e17 2.3061e17 2.2987e17 2.2875e17 2.2722e17 2.2593e17 2.2483e17 2.2353e17 2.22e17 2.1855e17 2.122e17 2.0624e17 1.9824e17 1.8853e17 1.7873e17 1.6836e17 1.5859e17 1.4939e17 1.407e17 1.3245e17 1.2457e17 1.1702e17 1.0977e17 1.0277e17 9.6011e16 8.9457e16; 2.3208e17 2.2853e17 2.2673e17 2.2535e17 2.2445e17 2.2374e17 2.2268e17 2.2124e17 2.2004e17 2.19e17 2.178e17 2.1638e17 2.1306e17 2.0687e17 2.0109e17 1.933e17 1.8381e17 1.7426e17 1.6414e17 1.546e17 1.4563e17 1.3715e17 1.291e17 1.2141e17 1.1404e17 1.0697e17 1.0014e17 9.3545e16 8.7152e16; 2.2602e17 2.2254e17 2.208e17 2.1949e17 2.1864e17 2.1796e17 2.1695e17 2.156e17 2.1448e17 2.1351e17 2.1238e17 2.1105e17 2.0786e17 2.0183e17 1.9621e17 1.8861e17 1.7935e17 1.7003e17 1.6015e17 1.5084e17 1.4207e17 1.3379e17 1.2593e17 1.1842e17 1.1123e17 1.0432e17 9.7657e16 9.1216e16 8.4974e16; 2.203e17 2.1689e17 2.1521e17 2.1397e17 2.1315e17 2.125e17 2.1153e17 2.1027e17 2.0922e17 2.083e17 2.0725e17 2.0601e17 2.0293e17 1.9706e17 1.9158e17 1.8417e17 1.7511e17 1.6602e17 1.5637e17 1.4726e17 1.387e17 1.3061e17 1.2292e17 1.1559e17 1.0856e17 1.0181e17 9.5303e16 8.901e16 8.2912e16; 2.149e17 2.1155e17 2.0992e17 2.0874e17 2.0796e17 2.0733e17 2.0641e17 2.0523e17 2.0423e17 2.0337e17 2.0239e17 2.0123e17 1.9826e17 1.9252e17 1.8719e17 1.7995e17 1.7109e17 1.6222e17 1.5278e17 1.4387e17 1.355e17 1.2759e17 1.2007e17 1.129e17 1.0603e17 9.9434e16 9.3071e16 8.6919e16 8.0958e16; 2.0978e17 2.0649e17 2.0491e17 2.0378e17 2.0303e17 2.0243e17 2.0155e17 2.0044e17 1.995e17 1.9869e17 1.9777e17 1.9668e17 1.9381e17 1.8821e17 1.83e17 1.7594e17 1.6727e17 1.5859e17 1.4936e17 1.4065e17 1.3246e17 1.2472e17 1.1736e17 1.1035e17 1.0363e17 9.7174e16 9.095e16 8.4933e16 7.9102e16; 2.0491e17 2.0169e17 2.0016e17 1.9908e17 1.9836e17 1.9778e17 1.9693e17 1.9589e17 1.95e17 1.9424e17 1.9338e17 1.9235e17 1.8957e17 1.8409e17 1.7902e17 1.7211e17 1.6362e17 1.5514e17 1.461e17 1.3757e17 1.2956e17 1.2198e17 1.1478e17 1.0792e17 1.0134e17 9.5024e16 8.8932e16 8.3043e16 7.7337e16; 2.0029e17 1.9712e17 1.9564e17 1.9461e17 1.9392e17 1.9336e17 1.9255e17 1.9156e17 1.9072e17 1.9e17 1.8919e17 1.8822e17 1.8552e17 1.8017e17 1.7522e17 1.6846e17 1.6014e17 1.5184e17 1.4299e17 1.3464e17 1.2679e17 1.1937e17 1.1232e17 1.056e17 9.9161e16 9.2973e16 8.7009e16 8.1242e16 7.5655e16; 1.9589e17 1.9278e17 1.9134e17 1.9035e17 1.8969e17 1.8915e17 1.8836e17 1.8743e17 1.8664e17 1.8595e17 1.8519e17 1.8428e17 1.8166e17 1.7642e17 1.7158e17 1.6496e17 1.5681e17 1.4869e17 1.4002e17 1.3184e17 1.2415e17 1.1688e17 1.0997e17 1.0339e17 9.7077e16 9.1015e16 8.5172e16 7.9523e16 7.4049e16; 1.917e17 1.8864e17 1.8724e17 1.8629e17 1.8565e17 1.8513e17 1.8437e17 1.8349e17 1.8274e17 1.8208e17 1.8136e17 1.805e17 1.7796e17 1.7283e17 1.681e17 1.6162e17 1.5363e17 1.4567e17 1.3717e17 1.2915e17 1.2162e17 1.1449e17 1.0772e17 1.0127e17 9.5084e16 8.9142e16 8.3415e16 7.7879e16 7.2515e16; 1.8769e17 1.8468e17 1.8332e17 1.8241e17 1.8179e17 1.8129e17 1.8056e17 1.7972e17 1.79e17 1.7838e17 1.777e17 1.7689e17 1.7441e17 1.6939e17 1.6475e17 1.5841e17 1.5057e17 1.4278e17 1.3445e17 1.2658e17 1.1919e17 1.122e17 1.0557e17 9.9236e16 9.3173e16 8.7348e16 8.1733e16 7.6305e16 7.1045e16; 1.8385e17 1.809e17 1.7957e17 1.7869e17 1.781e17 1.7761e17 1.769e17 1.761e17 1.7542e17 1.7483e17 1.7419e17 1.7341e17 1.71e17 1.6608e17 1.6155e17 1.5533e17 1.4764e17 1.4e17 1.3182e17 1.2411e17 1.1686e17 1.1001e17 1.035e17 9.7287e16 9.134e16 8.5626e16 8.0119e16 7.4795e16 6.9636e16; 1.8018e17 1.7727e17 1.7598e17 1.7513e17 1.7456e17 1.7409e17 1.734e17 1.7264e17 1.7199e17 1.7142e17 1.7081e17 1.7008e17 1.6773e17 1.629e17 1.5846e17 1.5236e17 1.4481e17 1.3732e17 1.293e17 1.2173e17 1.1462e17 1.0789e17 1.0151e17 9.5414e16 8.9579e16 8.3972e16 7.8568e16 7.3344e16 6.8283e16; 1.7665e17 1.7379e17 1.7253e17 1.7172e17 1.7116e17 1.707e17 1.7004e17 1.6931e17 1.6869e17 1.6814e17 1.6757e17 1.6686e17 1.6457e17 1.5984e17 1.5549e17 1.495e17 1.4209e17 1.3475e17 1.2688e17 1.1945e17 1.1246e17 1.0586e17 9.9591e16 9.3611e16 8.7884e16 8.2381e16 7.7076e16 7.1949e16 6.6981e16; 1.7326e17 1.7045e17 1.6922e17 1.6843e17 1.679e17 1.6745e17 1.668e17 1.661e17 1.6551e17 1.6499e17 1.6444e17 1.6377e17 1.6153e17 1.5688e17 1.5262e17 1.4675e17 1.3947e17 1.3226e17 1.2454e17 1.1724e17 1.1038e17 1.039e17 9.7745e16 9.1873e16 8.625e16 8.0847e16 7.5639e16 7.0605e16 6.5728e16; 1.6999e17 1.6723e17 1.6603e17 1.6527e17 1.6475e17 1.6431e17 1.6368e17 1.6302e17 1.6245e17 1.6195e17 1.6142e17 1.6078e17 1.5859e17 1.5403e17 1.4985e17 1.4409e17 1.3694e17 1.2986e17 1.2228e17 1.1511e17 1.0837e17 1.0201e17 9.5963e16 9.0196e16 8.4673e16 7.9367e16 7.4253e16 6.9309e16 6.4519e16; 1.6685e17 1.6413e17 1.6296e17 1.6222e17 1.6172e17 1.6129e17 1.6068e17 1.6004e17 1.5949e17 1.5901e17 1.5851e17 1.5789e17 1.5575e17 1.5127e17 1.4717e17 1.4151e17 1.3449e17 1.2754e17 1.2009e17 1.1305e17 1.0643e17 1.0018e17 9.4241e16 8.8576e16 8.315e16 7.7938e16 7.2913e16 6.8057e16 6.3352e16; 1.6381e17 1.6114e17 1.6e17 1.5928e17 1.5879e17 1.5837e17 1.5778e17 1.5716e17 1.5663e17 1.5617e17 1.5569e17 1.551e17 1.5301e17 1.4861e17 1.4458e17 1.3902e17 1.3212e17 1.253e17 1.1798e17 1.1106e17 1.0456e17 9.8411e16 9.2575e16 8.7008e16 8.1677e16 7.6555e16 7.1618e16 6.6847e16 6.2223e16; 1.6088e17 1.5825e17 1.5713e17 1.5644e17 1.5596e17 1.5555e17 1.5497e17 1.5438e17 1.5387e17 1.5342e17 1.5296e17 1.5239e17 1.5034e17 1.4602e17 1.4207e17 1.3661e17 1.2982e17 1.2312e17 1.1593e17 1.0913e17 1.0274e17 9.6696e16 9.0961e16 8.549e16 8.025e16 7.5216e16 7.0364e16 6.5675e16 6.1131e16]
    @test round.(zircon.alphaDeposition, sigdigits=5) ≈ alphaDeposition_known

    alphaDamage_known = [3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17 3.5402e17; 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17 3.4023e17; 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17 3.2745e17; 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17 3.1559e17; 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17 3.0457e17; 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17 2.9431e17; 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17 2.8474e17; 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17 2.7581e17; 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17 2.6745e17; 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17 2.5963e17; 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17 2.5228e17; 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17 2.4538e17; 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17 2.3888e17; 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17 2.3275e17; 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17 2.2696e17; 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17 2.2148e17; 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17 2.1628e17; 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17 2.1135e17; 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17 2.0665e17; 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17 2.0218e17; 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17 1.9791e17; 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17 1.9383e17; 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17 1.8992e17; 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17 1.8617e17; 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17 1.8257e17; 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17 1.791e17; 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17 1.7577e17; 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17 1.7255e17; 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17 1.6944e17; 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17 1.6644e17]
    @test round.(zircon.alphaDamage, sigdigits=5) ≈ alphaDamage_known


# Test whole integrated age program
    diffusionparams = (;
        DzEa = 165.0, # kJ/mol
        DzD0 = 193188.0, # cm^2/sec
        DN17Ea = 71.0, # kJ/mol
        DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
    )
    HeAgeSpherical(zircon,TSteps,pr,diffusionparams)
    @time "Running HeAgeSpherical" age = HeAgeSpherical(zircon,TSteps,pr,diffusionparams)
    @test age ≈ 520.0297717798045
    # Re-run to ensure internal state does not change
    for i=1:10
        @test HeAgeSpherical(zircon,TSteps,pr,diffusionparams) ≈ 520.0297717798045
    end

    crystalRadius = 35.
    Uppm = 1107.
    Thppm = 351.
    zircon = Zircon(crystalRadius,dr,Uppm,Thppm,dt,reverse(tSteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        @test HeAgeSpherical(zircon,TSteps,pr,diffusionparams) ≈ 309.7600561440283
    end

    crystalRadius = 135.
    Uppm = 1738.
    Thppm = 1171.
    zircon = Zircon(crystalRadius,dr,Uppm,Thppm,dt,reverse(tSteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        @test HeAgeSpherical(zircon,TSteps,pr,diffusionparams) ≈ 16.02209841621174
    end

    crystalRadius = 135.
    Uppm = 500.
    Thppm = 400.
    zircon = Zircon(crystalRadius,dr,Uppm,Thppm,dt,reverse(tSteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        @test HeAgeSpherical(zircon,TSteps,pr,diffusionparams) ≈ 777.5627957477788
    end
