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
    TSteps = collect(range(650, 0, length=length(tSteps)))

    pr_known = [0.000205838 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.00018167 0.000345328 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000175903 0.000306842 0.000575964 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000174306 0.00029782 0.000515098 0.000955418 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173858 0.000295383 0.000501104 0.000859803 0.00157683 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173735 0.000294718 0.00049742 0.000838273 0.00142761 0.00259 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173702 0.000294541 0.000496445 0.000832761 0.00139476 0.00235863 0.00423469 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173694 0.000294495 0.000496194 0.000831348 0.00138659 0.00230891 0.00387837 0.00689242 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173692 0.000294484 0.000496131 0.000830996 0.00138456 0.00229693 0.00380376 0.00634758 0.0111652 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496117 0.000830911 0.00138408 0.00229407 0.00378637 0.00623665 0.0103389 0.0179917 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496113 0.000830892 0.00138397 0.00229341 0.00378237 0.00621171 0.0101757 0.0167508 0.0288076 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830888 0.00138394 0.00229326 0.00378149 0.0062062 0.0101404 0.0165134 0.0269675 0.0457439 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229323 0.0037813 0.00620503 0.0101329 0.0164642 0.0266275 0.0430615 0.0718099 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378126 0.0062048 0.0101314 0.0164542 0.02656 0.0425841 0.0679927 0.110913 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378126 0.00620476 0.0101311 0.0164523 0.026547 0.0424937 0.0673402 0.105664 0.167408 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265447 0.0424772 0.0672227 0.104805 0.160534 0.244771 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265443 0.0424744 0.0672025 0.104658 0.15946 0.236352 0.34328 0 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671992 0.104635 0.159287 0.235101 0.333839 0.457702 0 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159261 0.234913 0.332512 0.448193 0.577183 0 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234886 0.332326 0.446936 0.568685 0.688768 0 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332302 0.446773 0.567636 0.682035 0.782538 0 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446753 0.56751 0.681264 0.777758 0.85451 0 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567497 0.681181 0.777256 0.851417 0.905927 0 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567496 0.681172 0.777207 0.851122 0.904069 0.940745 0 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851097 0.903911 0.939694 0.96343 0 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903899 0.939615 0.962863 0.977805 0 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962826 0.977511 0.986726 0 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977495 0.98658 0.992174 0 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977494 0.986574 0.992104 0.995453 0.0; 0.000173691 0.000294481 0.000496112 0.000830887 0.00138394 0.00229322 0.00378125 0.00620475 0.010131 0.0164519 0.0265442 0.0424739 0.0671987 0.104631 0.159257 0.234883 0.332299 0.446751 0.567495 0.681172 0.777203 0.851095 0.903898 0.939609 0.962824 0.977494 0.986573 0.992102 0.995421 0.9974]

    pr = DamageAnnealing(dt,tSteps,TSteps)

    @test size(pr) == (30,30)
    @test round.(pr, sigdigits=6) ≈ pr_known

# Test basic linear algebra
    dl = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    d = [1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 1.0]
    du = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    A = Tridiagonal(dl, d, du)
    y = [0.0; -0.900592; -2.72489; -4.54554; -6.36539; -8.18497; -10.0045; -11.8128; -13.5304; -15.2259; -16.8875; -18.5449; -20.1978; -21.804; -23.3226; -24.8305; -26.3291; -27.762; -29.1241; -30.2665; -31.0612; -31.7774; -32.1437; -32.1666; -31.9801; -31.5697; -31.0949; -30.565; -29.9796; -29.339; -28.6424; -27.8906; -27.0836; -26.2213; -25.3036; -24.3306; 0.0]

    sol_known = [-159.97020774647876, 159.97020774647876, 479.0100312394363, 795.3249647323938, 1107.0943582253512, 1412.4983617183084, 1709.7173952112655, 1996.9319287042226, 2272.33366219718, 2534.2049956901374, 2780.8504291830945, 3010.6083626760515, 3221.8213961690085, 3412.836629661965, 3582.0478631549227, 3727.93649664788, 3848.9946301408377, 3943.7236636337952, 4010.6906971267535, 4048.533630619712, 4056.110064112671, 4032.62529760563, 3977.3631310985884, 3889.9572645915464, 3770.3847980845044, 3618.832231577462, 3435.70996507042, 3221.492798563378, 2976.710632056336, 2701.948865549294, 2397.8480990422527, 2065.104932535211, 1704.4711660281687, 1316.7537995211264, 902.8151330140844, 463.5728665070422, 0.0]

    @test isapprox(A\y, sol_known, atol=0.00000000001)

# Test whole integrated age program
    diffusionparams = (;
        DzEa = 165.0, # kJ/mol
        DzD0 = 193188.0, # cm^2/sec
        DN17Ea = 71.0, # kJ/mol
        DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
    )

    crystalRadius = 29.26
    Uppm = 462.98
    Thppm = 177.76
    @test round(ZrnHeAgeSpherical(dt,reverse(tSteps),TSteps,pr,crystalRadius,dr,Uppm,Thppm,diffusionparams), sigdigits=5) ≈ 520.03

    crystalRadius = 35.
    Uppm = 1107.
    Thppm = 351.
    @test round(ZrnHeAgeSpherical(dt,reverse(tSteps),TSteps,pr,crystalRadius,dr,Uppm,Thppm,diffusionparams), sigdigits=5) ≈ 309.76
