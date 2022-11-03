
lambda = 10 ./ 10^5

@test SimAnnealSigma(1, lambda, 10, 100, 10) ≈ 110.44365174144839
@test SimAnnealSigma(10^5, lambda, 10, 100, 10) ≈ 14.145346247832224
