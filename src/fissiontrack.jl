## --- Fission track models
SEC_MYR = 1E6*365.25*24*3600

abstract type FissionTrackModel{T} end

struct FanningCurvilinear{T<:AbstractFloat} <: FissionTrackModel{T} 
    C0::T
    C1::T
    C2::T
    C3::T
    alpha::T
    beta::T
end

l₀ = 16.38 # [um] initial track length
σl₀ = 0.09 # [um] initial track length uncertainty


g(r, α ,β) = (((1-r^β)/β)^α-1)/α
r(g, α, β) = (1-(g*α+1)^(1/α)*β)^(1/β)

# From Ketcham, 1999
DR1999 = FanningCurvilinear(-106.1806145, 2.196512263, -155.9001311, -9.786405925, -0.48078, -6.3626)

# From Ketcham, 2007
DR = FanningCurvilinear(-62.8742, 1.3060, -85.4861, -8.3589, -0.3900, -9.1435)
RN = FanningCurvilinear(-41.2567, 0.8622, -65.8430, -7.9770, -0.3180, -3.1249)
B3 = FanningCurvilinear(-19.0257, 0.4135, -25.8627, -7.1209, -0.2187, -19.3036)
HS = FanningCurvilinear(-27.6610, 0.5880, -28.7268, -7.1267, -0.3720, -14.5778)

function equivalenttime(t, T, Teq, fc::FanningCurvilinear)
    exp(fc.C2 + (log(t*SEC_MYR)-fc.C2)*(log(1/(Teq+273.15))-fc.C3)/(log(1/(T+273.15))-fc.C3))/SEC_MYR
end

function reltracklength(t, T, fc::FanningCurvilinear{E}) where {E}
    # Relative track length 𝑟 (=l/l₀) following the 
    # Fanning Curvilinear equations of Ketcham et al. 1999
    g = fc.C0 + fc.C1*(log(t*SEC_MYR)-fc.C2)/(log(1/(T+273.15))-fc.C3)
    (1-max(g*fc.alpha+1, zero(E))^(1/fc.alpha)*fc.beta)^(1/fc.beta)
end

function reltrackdensity(r)
    # Ketchem et al., 2000 relationship between 
    # relative track density 𝑟 and relative track length ρ
    if r < 0.765
       9.205*r^2 - 9.157*r + 2.269
    else
        1.6*r-0.6
    end
end

## --- End of File