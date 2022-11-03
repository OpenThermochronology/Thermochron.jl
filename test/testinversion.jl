


simannealparams = (;
    ModelUncertainty = 10.0, # Ma
    InitialUncertainty = 100.0, # [Ma]
    lambda = 10 / 10^5, # [1/n]
)

@test SimAnnealSigma(1, 10, simannealparams) ≈ 110.44365174144839
@test SimAnnealSigma(10^5, 10, simannealparams) ≈ 14.145346247832224
