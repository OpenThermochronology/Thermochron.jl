## -- Argon age functions

# The amount of ingrown argon since time t
calc_Ar(t, K40) = K40*BR40K*(exp(λ40K*t)-1)
# First time derivative of the amount of ingrown argon since time t
calc_dArdt(t, K40) = K40*BR40K*λ40K*exp(λ40K*t)
# Use Newton's method to solve for Ar age
function newton_ar_age(Ar::T, K40; iterations::Int=16) where {T<:Number}
    Tf = float(T)
    argonage = one(Tf)
    for _ in 1:iterations
        ∂Ar∂t = calc_dArdt(argonage, K40) # Calculate derivative
        argonage += (Ar - calc_Ar(argonage, K40))/∂Ar∂t # Move towards zero (calc_Ar(argonage) == μAr)
    end
    return max(argonage, zero(Tf))
end
# Apply to an ArgonSample
function newton_age(mineral::ArgonSample)
    # Daughter concentrations, in atoms/gram
    μAr = final_diffusant(mineral)
    # Parent concentrations in atoms/gram
    μ40K = meanparent(mineral)
    # Numerically solve for raw Ar age of the grain (i.e., as measured)
    return newton_ar_age(μAr, μ40K)
end
meanparent(mineral::PlanarAr) = nanmean(mineral.r40K)
meanparent(mineral::SphericalAr) = nanmean(mineral.r40K, mineral.relvolumes)

## --- End of File