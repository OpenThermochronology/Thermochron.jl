
simannealmodel = (
    σModel = 10., # Model uncertainty [Ma]
    σAnnealing = 100., # Initial uncertainty [Ma]
    λAnnealing = 10 / 10^5, # lambda [1/n]
)

@test simannealsigma(1, 10; simannealmodel) ≈ 110.44365174144839
@test simannealsigma(10^5, 10; simannealmodel) ≈ 14.145346247832224
