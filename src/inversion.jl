
    """
    ```julia
    SimAnnealSigma(n, AnalyticalSigma, params)
    ```
    To avoid getting stuck in local optima, decrease uncertainty slowly by simulated annealing:
    Calculate annealing sigma which declines from InitialUncertainty to ModelUncertainty with
    a decay constant of lambda. Returns:

        sigma = sqrt(AnalyticalSigma^2 + AnnealingSigma^2)

    """
    function SimAnnealSigma(n, AnalyticalSigma, params)
        AnnealingSigma = params.InitialUncertainty*exp(-params.lambda*n) + params.ModelUncertainty
        sigma = sqrt(AnalyticalSigma^2 + AnnealingSigma^2)
        return sigma
    end
    export SimAnnealSigma
