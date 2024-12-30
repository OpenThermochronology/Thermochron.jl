## Test creating and allocating an GenericAr

    tCryst = 3000.0     # [Myr] Crystallization age
    dt = 100            # [Myr] time step size
    r = 29.26           # [μm] radius
    dr = 1              # [μm] radius step
    K40 = 14.0*1e4*κ40K # [ppm] K-40
    D0 = exp(4.094)     # [cm^2/s] frequency factor
    Ea = 49.22*4.184    # [kJ/mol] activation energy

    tsteps = (0+dt/2 : dt : tCryst-dt/2)
    Tsteps = collect(range(650, 0, length=length(tsteps)))

    GenericAr(r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    @time "Allocating a mineral" mineral = GenericAr(r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    @test isa(mineral, GenericAr)
    show(mineral)
    display(mineral)

    @test mineral.agesteps == reverse(tsteps)
    @test mineral.r40K ≈ fill(2.4623236248306448e17, 29)
    @test mineral.redges == 0:dr:r
    @test mineral.rsteps == (mineral.redges[2:end] + mineral.redges[1:end-1])/2
    @test mineral.nrsteps == length(mineral.rsteps) + 2 # Implicit radius steps inside and outside modeled range

    argondeposition_known = [7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15 7.3528e15; 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15 6.9555e15; 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15 6.5797e15; 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15 6.2242e15; 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15 5.8879e15; 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15 5.5698e15; 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15 5.2689e15; 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15 4.9842e15; 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15 4.7149e15; 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15 4.4601e15; 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15 4.2192e15; 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15 3.9912e15; 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15 3.7755e15; 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15 3.5715e15; 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15 3.3786e15; 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15 3.196e15; 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15 3.0233e15; 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15 2.86e15; 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15 2.7055e15; 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15 2.5593e15; 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15 2.421e15; 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15 2.2902e15; 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15 2.1665e15; 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15 2.0494e15; 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15 1.9387e15; 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15 1.8339e15; 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15 1.7348e15; 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15 1.6411e15; 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15 1.5524e15; 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15 1.4686e15]
    @test round.(mineral.argondeposition, sigdigits=5) ≈ argondeposition_known
    # println( round.(mineral.argondeposition, sigdigits=5))

## --- Test integrated age program for GenericAr

    modelage(mineral,Tsteps) # to not time compilation
    @time "Running modelage" age = modelage(mineral,Tsteps)
    @test age ≈ 904.2471019329474 
    # Re-run to ensure internal state does not change
    for _ in 1:4
        @test modelage(mineral,Tsteps) ≈  904.2471019329474 
    end

    r = 35.
    mineral = GenericAr(r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    # Re-run to ensure internal state does not change
    for _ in 1:4
        @test modelage(mineral,Tsteps) ≈ 917.8400803564615
    end

    r = 135.
    mineral = GenericAr(r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    # Re-run to ensure internal state does not change
    for _ in 1:4
        @test modelage(mineral,Tsteps) ≈ 1023.7342920730781
    end

    r = 1350.
    mineral = GenericAr(r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    # Re-run to ensure internal state does not change
    for _ in 1:4
        @test modelage(mineral,Tsteps) ≈ 1235.131761596538
    end

## --- As above but check calculated age as well

    r = 35.
    mineral = GenericAr(age=915, age_sigma=15, r=r,dr=dr,K40=K40,D0=D0,Ea=Ea,agesteps=reverse(tsteps))
    @test modelage(mineral,Tsteps) ≈ 917.8400803564615

## --- Test log likelihood

    @test Thermochron.model_ll(mineral,Tsteps) ≈ -3.6449133041539015

## --- End of file