## --- Test helium utilities

    @test Thermochron.calc_He(1000., 100, 100/137.818, 30) ≈ 151.8663003836161
    @test Thermochron.calc_He(1000., 100, 100/137.818, 30, 20) ≈ 151.99628115817814
    @test Thermochron.calc_dHedt(1000., 100, 100/137.818, 30) ≈ 0.16764721917364458
    @test Thermochron.calc_dHedt(1000., 100, 100/137.818, 30, 20) ≈ 0.16777762141120173
    @test Thermochron.newton_he_age(151.8663003836161, 100, 100/137.818, 30) ≈ 1000
    @test Thermochron.newton_he_age(151.99628115817814, 100, 100/137.818, 30, 20) ≈ 1000

## --- Test diff utilities

    for _ in 1:4
        local x = rand(100)
        local d = diff(x)
        @test Thermochron.maxdiff(x) === maximum(d)
        @test Thermochron.mindiff(x) === minimum(d)
        @test Thermochron.maxabsdiff(x) === maximum(abs.(d))
    end

    @test Thermochron.isdistinct(1:10, 5, 0.5)
    @test !Thermochron.isdistinct(1:10, 5, 1.5)
    @test Thermochron.isdistinct(collect(1:10), 5, 0.5, 10)
    @test !Thermochron.isdistinct(collect(1:10), 5, 1.5, 10)

## --- Test geometry utility functions

    # Intersection of two spheres
    @test size(Thermochron.intersectionfraction.(1.,1.,0:0.2:2)) == (11,)
    @test Thermochron.intersectionfraction.(1.,1.,0:0.2:2) ≈ [1//1, 9//20, 2//5, 7//20, 3//10, 1//4, 1//5, 3//20, 1//10, 1//20, 0//1]

    # Intersection of spherical shells
    crystalradius = 41.
    dr = 1.
    redges = (0. : dr : crystalradius) # Edges of each radius element
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element
    ralpha = 16.69

    # Zero when outside range
    @test Thermochron.intersectiondensity(redges,relvolumes,ralpha,60.) == zeros(length(relvolumes))
    # Single shell when centered
    dInt = Thermochron.intersectiondensity(redges,relvolumes,ralpha,0.)
    @test dInt[round(Int,ralpha)] > 0
    @test all(dInt[(1:length(relvolumes)) .!= round(Int,ralpha)] .== 0)
        # Specific value tests
    @test round.(Thermochron.intersectiondensity(redges,relvolumes,ralpha,0.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.3586, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]
    @test round.(Thermochron.intersectiondensity(redges,relvolumes,ralpha,5.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.81946, 11.0061, 10.1916, 9.48929, 8.87752, 8.33983, 7.86353, 7.43869, 7.05739, 6.71327, 4.38494, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]
    @test round.(Thermochron.intersectiondensity(redges,relvolumes,ralpha,20.), sigdigits=6) ≈ [0, 0, 0, 7.03672, 7.61584, 6.2396, 5.28378, 4.58152, 4.04384, 3.61901, 3.27489, 2.99049, 2.75152, 2.5479, 2.37232, 2.21938, 2.08496, 1.96588, 1.85967, 1.76435, 1.67832, 1.60028, 1.52918, 1.46413, 1.40439, 1.34933, 1.29842, 1.25122, 1.20733, 1.16641, 1.12817, 1.09236, 1.05876, 1.02716, 0.997389, 0.969297, 0.647731, 0, 0, 0, 0.0]
    @test round.(Thermochron.intersectiondensity(redges,relvolumes,ralpha,40.), sigdigits=6) ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.508458, 0.702195, 0.674665, 0.649212, 0.62561, 0.603663, 0.583204, 0.564086, 0.546181, 0.529379, 0.513579, 0.498694, 0.484649, 0.471372, 0.458804, 0.446888, 0.435576, 0.424822]

## --- Test basic linear algebra
    dl = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    d = [1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 1.0]
    du = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    A = Tridiagonal(dl, d, du)
    y = [0.0; -0.900592; -2.72489; -4.54554; -6.36539; -8.18497; -10.0045; -11.8128; -13.5304; -15.2259; -16.8875; -18.5449; -20.1978; -21.804; -23.3226; -24.8305; -26.3291; -27.762; -29.1241; -30.2665; -31.0612; -31.7774; -32.1437; -32.1666; -31.9801; -31.5697; -31.0949; -30.565; -29.9796; -29.339; -28.6424; -27.8906; -27.0836; -26.2213; -25.3036; -24.3306; 0.0]

    sol_known = [-159.97020774647876, 159.97020774647876, 479.0100312394363, 795.3249647323938, 1107.0943582253512, 1412.4983617183084, 1709.7173952112655, 1996.9319287042226, 2272.33366219718, 2534.2049956901374, 2780.8504291830945, 3010.6083626760515, 3221.8213961690085, 3412.836629661965, 3582.0478631549227, 3727.93649664788, 3848.9946301408377, 3943.7236636337952, 4010.6906971267535, 4048.533630619712, 4056.110064112671, 4032.62529760563, 3977.3631310985884, 3889.9572645915464, 3770.3847980845044, 3618.832231577462, 3435.70996507042, 3221.492798563378, 2976.710632056336, 2701.948865549294, 2397.8480990422527, 2065.104932535211, 1704.4711660281687, 1316.7537995211264, 902.8151330140844, 463.5728665070422, 0.0]

    @test isapprox(A\y, sol_known, atol=0.00000000001)

## ---