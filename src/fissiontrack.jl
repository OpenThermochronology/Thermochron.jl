## --- Fission track models

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


DR = FanningCurvilinear(-62.8742, 1.3060, -85.4861, -8.3589, -0.3900, -9.1435)
RN = FanningCurvilinear(-41.2567, 0.8622, -65.8430, -7.9770, -0.3180, -3.1249)
B3 = FanningCurvilinear(-19.0257, 0.4135, -25.8627, -7.1209, -0.2187, -19.3036)
HS = FanningCurvilinear(-27.6610, 0.5880, -28.7268, -7.1267, -0.3720, -14.5778)



## --- End of File