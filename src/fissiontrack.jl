## --- Fission track models
SEC_MYR = 1E6*365.25*24*3600

abstract type AnnealingModel{T} end

struct FanningCurvilinear{T<:AbstractFloat} <: AnnealingModel{T} 
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

using LsqFit: curve_fit

ellipse(x, lc) = @. sqrt(abs((1 - x^2/lc^2)*( 1.632*lc - 10.879)^2))
alrline(x, θalr) = @. (0.1035*θalr - 2.250) + x * tan(deg2rad(θalr))


"""
```julia
lcmodel(l, θ)
```
Calculate the model c-axis equivalent length ("lc,mod") following the approach
of Donelick et al. 1999 (doi: 10.2138/am-1999-0902) given a measured "confined"
fission track length `l` [microns] and angle from the c-axis `θ` [degrees].
"""
function lcmodel(l, θ)
    x = l*cos(deg2rad(θ))
    y = l*sin(deg2rad(θ))
    
    fobj = curve_fit(ellipse, Float64[x], Float64[y], Float64[len])
    lc_mod = only(fobj.param)
    # Return if we're above the minimum length for ALR
    lc_mod > 12.96 && return lc_mod

    # Return if we're below the minimum angle for ALR
    θalr = 0.304 * exp(0.439*lc_mod)
    θ < θalr && return lc_mod
 
    # Otherwise, fit to a linear segment for ALR
    fobj = curve_fit(alrline, [x], [y], Float64[θalr])

    lc_mod = log(only(fobj.param)/0.304)/0.439
    return lc_mod
end

"""
```julia
rmr0(F, Cl, OH, Mn=0, Fe=0, others=0)
```
Calculate rmr0 as a function of composition (specified in terms of
atoms per fomula unit, or APFU) for "multikinetic" apatite fission 
track thermochronology.

Implements the equation
```
rmr0 = (-0.0495 -0.0348F +0.3528|Cl - 1| +0.0701|OH - 1| 
        -0.8592Mn -1.2252Fe -0.1721Others)^0.1433
```
(equation 11) from Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
"""
function rmr0(F, Cl, OH, Mn=0, Fe=0, others=0)
    sum((F, Cl, OH)) ≈ 2 || error("F, Cl, and OH should sum to 2")
    h = - 0.0348F + 0.3528abs(Cl - 1) + 0.0701abs(OH - 1) 
        - 0.8592Mn - 1.2252Fe - 0.1721others -0.0495
    return h^0.1433
end


## --- 

abstract type FissionTrackLength{T} end

struct ApatiteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T       # [um]
    angle::T        # [degrees]
    Dpar::T         # [μm]
    F::T            # [APFU]
    Cl::T           # [APFU]
    OH::T           # [APFU]
end

abstract type FissionTrackAge{T} end

struct ApatiteFT{T<:AbstractFloat} <: FissionTrackAge{T}
    FTage::T        # [Ma]
    FTage_sigma::T  # [Ma]
    Dpar::T         # [μm]
    F::T            # [APFU]
    Cl::T           # [APFU]
    OH::T           # [APFU]
    UPbage::T       # [Ma]
    UPbage_sigma::T # [Ma]
end


## --- End of File