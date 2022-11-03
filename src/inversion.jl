
    """
    ```julia
    simannealsigma(n, sigma_analytical; params=(25.0, 35.0, 10/10^5))
    ```
    To avoid getting stuck in local optima, decrease uncertainty slowly by
    simulated annealing. Parameters are specified as a tuple `params` of the
    form (σₘ, σᵢ, λ), where annealing uncertainty declines from `σᵢ+σₘ` to `σₘ`
    with a decay constant of λ.

    Returns the annealing uncertainty added in quadrature with analytical
    uncertainty, or in other words

        sigma_annealing = σᵢ*exp(-λ*n) + σₘ
        sigma = sqrt(sigma_analytical^2 + sigma_annealing^2)

    """
    function simannealsigma(n, sigma_analytical; params=(25.0, 35.0, 10/10^5))
        σₘ, σᵢ, λ = params
        sigma_annealing = σᵢ*exp(-λ*n) + σₘ
        return sqrt(sigma_analytical^2 + sigma_annealing^2)
    end
    export simannealsigma
