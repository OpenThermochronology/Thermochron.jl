
simannealparams = (
    10., # Model uncertainty [Ma]
    100., # Initial uncertainty [Ma]
    10 / 10^5, # lambda [1/n]
)

@test simannealsigma(1, 10, params=simannealparams) ≈ 110.44365174144839
@test simannealsigma(10^5, 10, params=simannealparams) ≈ 14.145346247832224
