# Test damage annealing function for RDAAM

    dm = RDAAM(rmr0=0.0, kappa=1.0) # RDAAM, without rmr0 correction
    tCryst = 3000.0 # Time (in AMyr)
    dt = 100 # time step size in Myr
    dr = 1
    tsteps = collect(0+dt/2 : dt : tCryst-dt/2)

    Tsteps = fill(250.0, length(tsteps))
    pr_known = [0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0 0; 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0 0; 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0 0; 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0 0; 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0 0; 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0 0; 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0 0; 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0 0; 0.193029 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0 0; 0.192268 0.193029 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0 0; 0.191537 0.192268 0.193029 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0 0; 0.190834 0.191537 0.192268 0.193029 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174 0; 0.190157 0.190834 0.191537 0.192268 0.193029 0.193822 0.194651 0.195517 0.196426 0.197381 0.198386 0.199447 0.200571 0.201764 0.203037 0.204398 0.205861 0.207442 0.209161 0.211042 0.213118 0.215431 0.218041 0.221029 0.224518 0.228699 0.233897 0.240729 0.250609 0.268174]
    
    pr, Teq = anneal(dt, tsteps, Tsteps, dm)

    @test isa(pr, AbstractMatrix)
    @test size(pr) == (30,30)
    @test round.(pr, sigdigits=6) ≈ pr_known
    @test isa(Teq, AbstractVector)
    @test length(Teq) == 30

    for i=1:10
        anneal!(pr, Teq, dt, tsteps, Tsteps, dm)
        @test round.(pr, sigdigits=6) ≈ pr_known
    end

    Tsteps = collect(range(650, 0, length=length(tsteps)))
    pr_known = [0.00089077 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000877932 0.00125167 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876788 0.00123386 0.00175621 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876687 0.00123228 0.00173155 0.00246073 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123214 0.00172939 0.00242672 0.00344343 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.0017292 0.00242376 0.00339665 0.00481261 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242351 0.00339263 0.00474848 0.00671799 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339229 0.00474304 0.00663038 0.00936595 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339227 0.00474259 0.00662306 0.00924671 0.0130397 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662247 0.00923693 0.0128782 0.0181258 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923614 0.0128652 0.017908 0.0251467 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178908 0.0248548 0.0347997 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178895 0.0248323 0.0344116 0.0479982 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248306 0.0343825 0.0474872 0.0659054 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343804 0.0474499 0.0652409 0.0899424 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474473 0.0651938 0.0890922 0.121737 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651907 0.089034 0.120672 0.16297 0 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.0890302 0.120602 0.161671 0.215064 0 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120598 0.161589 0.213534 0.278713 0 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213442 0.276987 0.353318 0 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276888 0.351469 0.436558 0 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351369 0.434693 0.524391 0 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351365 0.434599 0.522631 0.61168 0 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522548 0.610132 0.693315 0 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610064 0.692046 0.765336 0 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610062 0.691996 0.764365 0.825582 0 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610062 0.691994 0.76433 0.824886 0.873691 0 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610062 0.691994 0.764329 0.824863 0.873219 0.910624 0 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610062 0.691994 0.764329 0.824863 0.873206 0.910322 0.938065 0; 0.000876678 0.00123213 0.00172919 0.00242349 0.00339226 0.00474255 0.00662242 0.00923608 0.0128641 0.0178894 0.0248305 0.0343802 0.0474471 0.0651905 0.08903 0.120597 0.161584 0.213437 0.276883 0.351364 0.434595 0.522544 0.610062 0.691994 0.764329 0.824863 0.873206 0.910315 0.937882 0.957907]
    
    pr, Teq = anneal(dt, tsteps, Tsteps, dm)

    @test isa(pr, AbstractMatrix)
    @test size(pr) == (30,30)
    @test round.(pr, sigdigits=6) ≈ pr_known
    @test isa(Teq, AbstractVector)
    @test length(Teq) == 30

    for i=1:10
        anneal!(pr, Teq, dt, tsteps, Tsteps, dm)
        @test round.(pr, sigdigits=6) ≈ pr_known
    end


# Test creating and allocating an ApatiteHe

    crystalradius = 29.26
    Uppm = 33.0
    Thppm = 24.2
    ApatiteHe(crystalradius,dr,Uppm,Thppm,dt,reverse(tsteps))
    @time "Allocating an apatite" apatite = ApatiteHe(crystalradius,dr,Uppm,Thppm,dt,reverse(tsteps))
    @test isa(apatite, ApatiteHe)
    @test apatite.agesteps == reverse(tsteps)
    @test apatite.r238U ≈ fill(8.349831932773109e16, 29)
    @test apatite.r235U ≈ fill(6.13593691093681e14, 29)
    @test apatite.r232Th ≈ fill(6.281568965517242e16, 29)
    @test apatite.redges == 0:dr:crystalradius
    @test apatite.rsteps == (apatite.redges[2:end] + apatite.redges[1:end-1])/2
    @test apatite.nrsteps == length(apatite.rsteps) + 2 # Implicit radius steps inside and outside modeled range
    @test round.(apatite.relvolumes, sigdigits=5) ≈ [4.1002e-5, 0.00028701, 0.00077904, 0.0015171, 0.0025011, 0.0037312, 0.0052073, 0.0069294, 0.0088975, 0.011112, 0.013572, 0.016278, 0.01923, 0.022428, 0.025872, 0.029563, 0.033499, 0.037681, 0.042109, 0.046783, 0.051704, 0.05687, 0.062282, 0.06794, 0.073845, 0.079995, 0.086391, 0.093034, 0.099922]
    @test round.(apatite.r238UHe, sigdigits=5) ≈ [6.0329e17, 5.8732e17, 5.8564e17, 5.8516e17, 5.8592e17, 5.9207e17, 5.9279e17, 5.8884e17, 5.8546e17, 5.7858e17, 5.6654e17, 5.5386e17, 5.4263e17, 5.2168e17, 4.9705e17, 4.7442e17, 4.4831e17, 4.2407e17, 4.0151e17, 3.8034e17, 3.6032e17, 3.4132e17, 3.2321e17, 3.0586e17, 2.8919e17, 2.7311e17, 2.5756e17, 2.4248e17, 2.2782e17]
    @test round.(apatite.r235UHe, sigdigits=5) ≈ [3.5158e15, 3.6585e15, 3.6284e15, 3.5583e15, 3.4902e15, 3.4403e15, 3.4004e15, 3.3132e15, 3.2163e15, 3.1108e15, 3.0195e15, 2.9387e15, 2.8288e15, 2.7254e15, 2.631e15, 2.5037e15, 2.3853e15, 2.2747e15, 2.1705e15, 2.0719e15, 1.9782e15, 1.8884e15, 1.8022e15, 1.7189e15, 1.6384e15, 1.5604e15, 1.4845e15, 1.4104e15, 1.3379e15]
    @test round.(apatite.r232ThHe, sigdigits=5) ≈ [3.0641e17, 3.1219e17, 2.9997e17, 2.9417e17, 2.9066e17, 2.9001e17, 2.9181e17, 2.9286e17, 2.9035e17, 2.8328e17, 2.7161e17, 2.6069e17, 2.5099e17, 2.4218e17, 2.3414e17, 2.267e17, 2.1929e17, 2.0857e17, 1.9851e17, 1.8902e17, 1.7999e17, 1.7135e17, 1.6307e17, 1.5511e17, 1.4743e17, 1.3998e17, 1.3275e17, 1.2572e17, 1.1884e17]

    alphadeposition_known = [2.2867e16 2.2765e16 2.26e16 2.2429e16 2.2305e16 2.2362e16 2.2319e16 2.2071e16 2.1799e16 2.14e16 2.0874e16 2.0356e16 1.9827e16 1.9077e16 1.8257e16 1.7431e16 1.6535e16 1.5681e16 1.4883e16 1.4132e16 1.3421e16 1.2745e16 1.2098e16 1.1478e16 1.088e16 1.0303e16 9.7439e15 9.2006e15 8.6716e15; 2.2037e16 2.1917e16 2.1758e16 2.1599e16 2.1487e16 2.155e16 2.1513e16 2.1281e16 2.1027e16 2.0649e16 2.0143e16 1.9643e16 1.9138e16 1.8413e16 1.7619e16 1.6823e16 1.5957e16 1.5131e16 1.4359e16 1.3634e16 1.2946e16 1.2292e16 1.1667e16 1.1068e16 1.049e16 9.9325e15 9.3923e15 8.8675e15 8.3564e15; 2.1266e16 2.1131e16 2.0977e16 2.0829e16 2.0727e16 2.0795e16 2.0764e16 2.0547e16 2.0309e16 1.995e16 1.9463e16 1.898e16 1.8496e16 1.7796e16 1.7025e16 1.6257e16 1.542e16 1.462e16 1.3873e16 1.317e16 1.2505e16 1.1872e16 1.1267e16 1.0687e16 1.0128e16 9.5884e15 9.0658e15 8.5581e15 8.0638e15; 2.055e16 2.04e16 2.0251e16 2.0114e16 2.002e16 2.0094e16 2.0067e16 1.9864e16 1.9641e16 1.9299e16 1.8829e16 1.8363e16 1.7898e16 1.722e16 1.6472e16 1.573e16 1.4919e16 1.4143e16 1.3419e16 1.2738e16 1.2094e16 1.148e16 1.0894e16 1.0332e16 9.7908e15 9.2682e15 8.762e15 8.2703e15 7.7916e15; 1.9882e16 1.972e16 1.9576e16 1.9448e16 1.9363e16 1.944e16 1.9419e16 1.9227e16 1.9018e16 1.8692e16 1.8238e16 1.7787e16 1.7341e16 1.6683e16 1.5956e16 1.5238e16 1.4452e16 1.3699e16 1.2996e16 1.2336e16 1.171e16 1.1115e16 1.0547e16 1.0002e16 9.4767e15 8.9699e15 8.4791e15 8.0023e15 7.5382e15; 1.9259e16 1.9086e16 1.8947e16 1.8827e16 1.875e16 1.8831e16 1.8813e16 1.8633e16 1.8436e16 1.8125e16 1.7686e16 1.7249e16 1.6819e16 1.6181e16 1.5474e16 1.4779e16 1.4015e16 1.3284e16 1.2601e16 1.196e16 1.1352e16 1.0774e16 1.0223e16 9.6932e15 9.1835e15 8.6915e15 8.215e15 7.7522e15 7.3017e15; 1.8678e16 1.8495e16 1.836e16 1.8247e16 1.8177e16 1.8261e16 1.8247e16 1.8078e16 1.7892e16 1.7594e16 1.7169e16 1.6745e16 1.6331e16 1.5711e16 1.5022e16 1.4349e16 1.3606e16 1.2895e16 1.2232e16 1.1608e16 1.1017e16 1.0456e16 9.9194e15 9.4048e15 8.9095e15 8.4313e15 7.9683e15 7.5186e15 7.0809e15; 1.8133e16 1.7943e16 1.7811e16 1.7705e16 1.7641e16 1.7728e16 1.7718e16 1.7558e16 1.7382e16 1.7096e16 1.6685e16 1.6273e16 1.5873e16 1.527e16 1.4599e16 1.3945e16 1.3223e16 1.2531e16 1.1885e16 1.1278e16 1.0704e16 1.0157e16 9.6353e15 9.1347e15 8.6528e15 8.1877e15 7.7373e15 7.2999e15 6.8742e15; 1.7623e16 1.7425e16 1.7297e16 1.7198e16 1.7139e16 1.7228e16 1.7221e16 1.707e16 1.6903e16 1.6629e16 1.623e16 1.5829e16 1.5443e16 1.4856e16 1.4201e16 1.3566e16 1.2863e16 1.2189e16 1.156e16 1.0969e16 1.0409e16 9.8768e15 9.3688e15 8.8813e15 8.4121e15 7.9592e15 7.5207e15 7.0949e15 6.6805e15; 1.7143e16 1.694e16 1.6815e16 1.6722e16 1.6668e16 1.6759e16 1.6754e16 1.6611e16 1.6453e16 1.619e16 1.5802e16 1.5412e16 1.5038e16 1.4467e16 1.3827e16 1.321e16 1.2525e16 1.1868e16 1.1254e16 1.0678e16 1.0132e16 9.6134e15 9.1183e15 8.6431e15 8.1859e15 7.7446e15 7.3173e15 6.9024e15 6.4986e15; 1.6692e16 1.6485e16 1.6363e16 1.6274e16 1.6225e16 1.6317e16 1.6315e16 1.6179e16 1.603e16 1.5776e16 1.5398e16 1.5019e16 1.4656e16 1.4099e16 1.3475e16 1.2874e16 1.2206e16 1.1564e16 1.0966e16 1.0404e16 9.8716e15 9.3654e15 8.8824e15 8.419e15 7.973e15 7.5426e15 7.1259e15 6.7213e15 6.3275e15; 1.6268e16 1.6056e16 1.5937e16 1.5853e16 1.5808e16 1.5901e16 1.5901e16 1.5772e16 1.563e16 1.5385e16 1.5018e16 1.4647e16 1.4296e16 1.3752e16 1.3142e16 1.2557e16 1.1905e16 1.1278e16 1.0694e16 1.0145e16 9.6256e15 9.1315e15 8.66e15 8.2076e15 7.7723e15 7.3522e15 6.9455e15 6.5507e15 6.1664e15; 1.5867e16 1.5652e16 1.5535e16 1.5455e16 1.5414e16 1.5508e16 1.5511e16 1.5388e16 1.5252e16 1.5016e16 1.4658e16 1.4296e16 1.3955e16 1.3424e16 1.2827e16 1.2257e16 1.162e16 1.1008e16 1.0437e16 9.9008e15 9.3931e15 8.9104e15 8.4498e15 8.0079e15 7.5828e15 7.1724e15 6.7752e15 6.3896e15 6.0143e15; 1.5487e16 1.527e16 1.5156e16 1.508e16 1.5042e16 1.5137e16 1.5141e16 1.5024e16 1.4894e16 1.4666e16 1.4317e16 1.3964e16 1.3632e16 1.3113e16 1.2529e16 1.1972e16 1.135e16 1.0752e16 1.0194e16 9.6694e15 9.173e15 8.7012e15 8.2509e15 7.819e15 7.4034e15 7.0024e15 6.6141e15 6.2373e15 5.8705e15; 1.5128e16 1.4909e16 1.4797e16 1.4725e16 1.469e16 1.4785e16 1.4791e16 1.4679e16 1.4555e16 1.4334e16 1.3993e16 1.3648e16 1.3325e16 1.2818e16 1.2246e16 1.1702e16 1.1094e16 1.0509e16 9.9628e15 9.4498e15 8.9643e15 8.5028e15 8.0624e15 7.6399e15 7.2335e15 6.8412e15 6.4615e15 6.093e15 5.7343e15; 1.4788e16 1.4567e16 1.4457e16 1.4388e16 1.4356e16 1.4452e16 1.4459e16 1.4352e16 1.4233e16 1.4019e16 1.3685e16 1.3348e16 1.3033e16 1.2537e16 1.1977e16 1.1446e16 1.0851e16 1.0278e16 9.7434e15 9.2413e15 8.7661e15 8.3143e15 7.8833e15 7.4699e15 7.0721e15 6.6883e15 6.3167e15 5.9561e15 5.6051e15; 1.4464e16 1.4242e16 1.4134e16 1.4068e16 1.4039e16 1.4135e16 1.4143e16 1.404e16 1.3927e16 1.3718e16 1.3393e16 1.3062e16 1.2756e16 1.227e16 1.1721e16 1.1202e16 1.0619e16 1.0058e16 9.5345e15 9.0428e15 8.5774e15 8.135e15 7.7129e15 7.3081e15 6.9186e15 6.5428e15 6.179e15 5.8259e15 5.4823e15; 1.4155e16 1.3933e16 1.3827e16 1.3764e16 1.3737e16 1.3833e16 1.3842e16 1.3744e16 1.3634e16 1.3432e16 1.3113e16 1.279e16 1.249e16 1.2015e16 1.1477e16 1.0969e16 1.0398e16 9.8481e15 9.3353e15 8.8535e15 8.3976e15 7.9641e15 7.5506e15 7.154e15 6.7724e15 6.4043e15 6.0479e15 5.702e15 5.3654e15; 1.3861e16 1.3639e16 1.3535e16 1.3474e16 1.3449e16 1.3545e16 1.3555e16 1.3461e16 1.3355e16 1.3159e16 1.2846e16 1.253e16 1.2237e16 1.1771e16 1.1243e16 1.0746e16 1.0187e16 9.6478e15 9.1451e15 8.6728e15 8.2259e15 7.801e15 7.3957e15 7.0069e15 6.633e15 6.2721e15 5.9228e15 5.5838e15 5.254e15; 1.358e16 1.3358e16 1.3256e16 1.3198e16 1.3174e16 1.3269e16 1.3281e16 1.319e16 1.3088e16 1.2897e16 1.2591e16 1.2281e16 1.1995e16 1.1538e16 1.102e16 1.0533e16 9.985e15 9.4562e15 8.9632e15 8.5e15 8.0617e15 7.645e15 7.2475e15 6.8664e15 6.4996e15 6.1458e15 5.8033e15 5.4709e15 5.1475e15; 1.3312e16 1.309e16 1.2989e16 1.2933e16 1.2912e16 1.3006e16 1.3018e16 1.2931e16 1.2833e16 1.2646e16 1.2347e16 1.2042e16 1.1762e16 1.1314e16 1.0806e16 1.0329e16 9.7915e15 9.2726e15 8.7889e15 8.3344e15 7.9044e15 7.4956e15 7.1057e15 6.7318e15 6.372e15 6.0249e15 5.689e15 5.3629e15 5.0457e15; 1.3054e16 1.2834e16 1.2734e16 1.268e16 1.266e16 1.2754e16 1.2767e16 1.2683e16 1.2588e16 1.2406e16 1.2112e16 1.1813e16 1.1539e16 1.1099e16 1.0601e16 1.0133e16 9.6057e15 9.0963e15 8.6216e15 8.1756e15 7.7535e15 7.3523e15 6.9697e15 6.6027e15 6.2497e15 5.909e15 5.5794e15 5.2594e15 4.9481e15; 1.2807e16 1.2588e16 1.249e16 1.2437e16 1.2418e16 1.2512e16 1.2525e16 1.2444e16 1.2352e16 1.2174e16 1.1886e16 1.1593e16 1.1325e16 1.0893e16 1.0403e16 9.9444e15 9.4271e15 8.927e15 8.4608e15 8.023e15 7.6086e15 7.2147e15 6.839e15 6.4788e15 6.1322e15 5.7978e15 5.4741e15 5.1601e15 4.8545e15; 1.2569e16 1.2352e16 1.2255e16 1.2204e16 1.2186e16 1.2279e16 1.2293e16 1.2215e16 1.2126e16 1.1952e16 1.1669e16 1.1381e16 1.1118e16 1.0694e16 1.0213e16 9.763e15 9.2551e15 8.764e15 8.3062e15 7.8761e15 7.4691e15 7.0823e15 6.7133e15 6.3596e15 6.0192e15 5.6908e15 5.373e15 5.0646e15 4.7645e15; 1.2341e16 1.2124e16 1.2029e16 1.198e16 1.1963e16 1.2055e16 1.207e16 1.1994e16 1.1907e16 1.1738e16 1.146e16 1.1176e16 1.0919e16 1.0503e16 1.003e16 9.5882e15 9.0894e15 8.6069e15 8.1571e15 7.7346e15 7.3348e15 6.9548e15 6.5923e15 6.2447e15 5.9104e15 5.5878e15 5.2756e15 4.9727e15 4.6779e15; 1.212e16 1.1906e16 1.1812e16 1.1764e16 1.1748e16 1.184e16 1.1854e16 1.1781e16 1.1697e16 1.1531e16 1.1257e16 1.0979e16 1.0726e16 1.0318e16 9.8528e15 9.4194e15 8.9295e15 8.4552e15 8.0133e15 7.5981e15 7.2052e15 6.8317e15 6.4756e15 6.134e15 5.8055e15 5.4885e15 5.1818e15 4.8841e15 4.5944e15; 1.1908e16 1.1695e16 1.1602e16 1.1555e16 1.1541e16 1.1631e16 1.1646e16 1.1575e16 1.1493e16 1.1331e16 1.1062e16 1.0788e16 1.0541e16 1.0139e16 9.6819e15 9.2562e15 8.7749e15 8.3087e15 7.8743e15 7.4662e15 7.08e15 6.7129e15 6.3628e15 6.0271e15 5.7042e15 5.3927e15 5.0911e15 4.7986e15 4.5139e15; 1.1702e16 1.1492e16 1.14e16 1.1354e16 1.134e16 1.143e16 1.1445e16 1.1376e16 1.1297e16 1.1137e16 1.0873e16 1.0604e16 1.0361e16 9.9658e15 9.5164e15 9.0984e15 8.6254e15 8.167e15 7.7398e15 7.3386e15 6.9589e15 6.598e15 6.2538e15 5.9238e15 5.6063e15 5.3e15 5.0036e15 4.7159e15 4.4361e15; 1.1503e16 1.1295e16 1.1204e16 1.116e16 1.1146e16 1.1235e16 1.1251e16 1.1183e16 1.1106e16 1.095e16 1.069e16 1.0425e16 1.0186e16 9.7981e15 9.3562e15 8.9454e15 8.4805e15 8.0297e15 7.6096e15 7.215e15 6.8416e15 6.4867e15 6.1483e15 5.8237e15 5.5115e15 5.2103e15 4.9188e15 4.636e15 4.3608e15; 1.1311e16 1.1104e16 1.1015e16 1.0971e16 1.0959e16 1.1046e16 1.1062e16 1.0997e16 1.0922e16 1.0768e16 1.0513e16 1.0252e16 1.0017e16 9.6355e15 9.2008e15 8.7971e15 8.34e15 7.8966e15 7.4834e15 7.0953e15 6.728e15 6.3789e15 6.046e15 5.7267e15 5.4197e15 5.1234e15 4.8367e15 4.5585e15 4.2879e15]
    @test round.(apatite.alphadeposition, sigdigits=5) ≈ alphadeposition_known
    # println( round.(apatite.alphadeposition, sigdigits=5))

    alphadamage_known = [2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16 2.6258e16; 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16 2.527e16; 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16 2.4354e16; 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16 2.3504e16; 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16 2.2713e16; 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16 2.1977e16; 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16 2.129e16; 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16 2.0649e16; 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16 2.0048e16; 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16 1.9486e16; 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16 1.8958e16; 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16 1.8461e16; 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16 1.7993e16; 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16 1.7551e16; 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16 1.7134e16; 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16 1.6738e16; 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16 1.6363e16; 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16 1.6007e16; 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16 1.5667e16; 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16 1.5344e16; 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16 1.5035e16; 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16 1.4739e16; 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16 1.4456e16; 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16 1.4184e16; 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16 1.3923e16; 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16 1.3672e16; 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16 1.343e16; 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16 1.3196e16; 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16 1.297e16; 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16 1.2751e16]
    @test round.(apatite.alphadamage, sigdigits=5) ≈ alphadamage_known
    # println( round.(apatite.alphadamage, sigdigits=5))


# Test integrated age program for ApatiteHe

    HeAgeSpherical(apatite,Tsteps,pr,dm) # to not time compilation
    @time "Running HeAgeSpherical" age = HeAgeSpherical(apatite,Tsteps,pr,dm)
    @test age ≈ 125.2393082743598 
    dm = RDAAM() # Standard RDAAM
    pr, Teq = anneal(dt, tsteps, Tsteps, dm)
    # Re-run to ensure internal state does not change
    for i=1:10
        # @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 125.2393082743598  # with no rmr0 correction
        @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 125.2331472062392  # with rmr0 = 0.83
    end

    crystalradius = 35.
    Uppm = 110.7
    Thppm = 35.1
    apatite = ApatiteHe(crystalradius,dr,Uppm,Thppm,dt,reverse(tsteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        # @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 150.6059246264699  # with no rmr0 correction
        @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 150.3786208747766  # with rmr0 = 0.83
    end

    crystalradius = 135.
    Uppm = 173.8
    Thppm = 117.1
    apatite = ApatiteHe(crystalradius,dr,Uppm,Thppm,dt,reverse(tsteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        # @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 266.3624679808305  # with no rmr0 correction
        @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 263.9370169558578  # with rmr0 = 0.83
    end

    crystalradius = 135.
    Uppm = 50.0
    Thppm = 40.0
    apatite = ApatiteHe(crystalradius,dr,Uppm,Thppm,dt,reverse(tsteps))
    # Re-run to ensure internal state does not change
    for i=1:10
        # @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 263.9393563618151  # with no rmr0 correction
        @test HeAgeSpherical(apatite,Tsteps,pr,dm) ≈ 263.87540959714664  # with rmr0 = 0.83
    end
