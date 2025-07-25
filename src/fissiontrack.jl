## --- Concrete types for fission track annealing models

"""
```julia
Yamada2007PC(
    c0p = -63.37     # Yamada et al. 2007 zircon
    c1p = 0.212      # Yamada et al. 2007 zircon
    bp = 43.00       # Yamada et al. 2007 zircon
)
```
Parallel Curvilinear zircon annealing model of Yamada, 2007 
(doi: 10.1016/j.chemgeo.2006.09.002)
"""
Base.@kwdef struct Yamada2007PC{T<:AbstractFloat} <: ZirconAnnealingModel{T} 
    c0p::T = -63.37 - 0.212log(3600) # Yamada et al. 2007 zircon, converted from log(h) to log(s)
    c1p::T = 0.212      # Yamada et al. 2007 zircon
    bp::T = 43.00       # Yamada et al. 2007 zircon
end

"""
```julia
Guenthner2013FC(
    C0 = 6.24534         # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C1 = -0.11977        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C2 = -314.93688      # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C3 = -14.2868        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    alpha = -0.057206897 # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
)
```
Fanning Curvilinear zircon annealing model with simplified Box-Cox transform 
from Guenthner et al. 2013 (doi: 10.2475/03.2013.01)
"""
Base.@kwdef struct Guenthner2013FC{T<:AbstractFloat} <: FanningCurvilinearZircon{T} 
    C0::T = 6.24534         # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C1::T = -0.11977        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C2::T = -314.93688      # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C3::T = -14.2868        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    alpha::T = -0.057206897 # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
end

"""
```julia
Ketcham1999FC(
    C0::T = -19.844     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C1::T = 0.38951     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C2::T = -51.253     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C3::T = -7.6423     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    alpha::T = -0.12327 # Box-Cox transform parameter
    beta::T = -11.988   # Box-Cox transform parameter
)
```
Fanning Curvilinear apatite annealing model from Ketcham, 1999 (doi: 10.2138/am-1999-0903)
"""
Base.@kwdef struct Ketcham1999FC{T<:AbstractFloat} <: FanningCurvilinearApatite{T} 
    C0::T = -19.844     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C1::T = 0.38951     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C2::T = -51.253     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C3::T = -7.6423     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    alpha::T = -0.12327 # Box-Cox transform parameter
    beta::T = -11.988   # Box-Cox transform parameter
end

"""
```julia
Ketcham2007FC(
    C0 = 0.39528     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C1 = 0.01073     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C2 = -65.12969   # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C3 = -7.91715    # "Simultaneous fit" from Ketcham et al. 2007 apatite
    alpha = 0.04672  # "Simultaneous fit" from Ketcham et al. 2007 apatite
)
```
Fanning Curvilinear apatite annealing model with simplified Box-Cox transform 
from Ketcham, 2007 (doi: 10.2138/am.2007.2281)
"""
Base.@kwdef struct Ketcham2007FC{T<:AbstractFloat} <: FanningCurvilinearApatite{T} 
    C0::T = 0.39528     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C1::T = 0.01073     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C2::T = -65.12969   # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C3::T = -7.91715    # "Simultaneous fit" from Ketcham et al. 2007 apatite
    alpha::T = 0.04672  # "Simultaneous fit" from Ketcham et al. 2007 apatite
end

"""
```julia
Jones2021FA(
    C0 = 1.374           # Annealing parameter
    C1 = -4.192812e-5    # Annealing parameter
    C2 = -22.70885029    # Annealing parameter
    C3 = 0.0             # Annealing parameter
)
```
Parallel Arrhenius monazite annealing model modified from
Jones et al., 2021 (doi: 10.5194/gchron-3-89-2021)
"""
Base.@kwdef struct Jones2021FA{T<:AbstractFloat} <: MonaziteAnnealingModel{T} 
    C0::T = 1.374           # Annealing parameter
    C1::T = -4.192812e-5    # Annealing parameter
    C2::T = -22.70885029    # Annealing parameter
    C3::T = 0.0             # Annealing parameter
end

## --- Fission track functions

# Fanning Curvilinear
function equivalenttime(t::Number, T::Number, Teq::Number, am::FanningCurvilinear)
    exp(am.C2 + (log(t*SEC_MYR)-am.C2)*(log(1/(Teq+273.15))-am.C3)/(log(1/(T+273.15))-am.C3))/SEC_MYR
end
# Fanning Linear (Fanning Arrhenius)
function equivalenttime(t::Number, T::Number, Teq::Number, am::Jones2021FA)
    exp(am.C2 + (log(t*SEC_MYR)-am.C2)*((1/(Teq+273.15))-am.C3)/((1/(T+273.15))-am.C3))/SEC_MYR
end
# Parallel Curvilinear
function equivalenttime(t::Number, T::Number, Teq::Number, am::Yamada2007PC)
    exp(am.bp*log(1/(Teq+273.15)) + (log(t*SEC_MYR) - am.bp*log(1/(T+273.15))))/SEC_MYR
end

"""
```julia
reltracklength(t, T, am::AnnealingModel) 
```
Calculate the relative track length `r` (equal to `l/l₀`) expected after 
isothermal heating at `T` C for `t` Myr and annealing parameters `am`. 
    
Possible annealing model types and the references for t he equations 
which they respetively implement include 
  `Ketcham1999FC`       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `Ketcham2007FC`       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
  `Yamada2007PC`        Parallel Curvilinear zircon model of Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002)

See also: `reltrackdensity`.
"""
# Fanning Curvilinear, Box-Cox
function reltracklength(t::Number, T::Number, am::Ketcham1999FC{F}) where {F<:AbstractFloat}
    g = am.C0 + am.C1*(log(t*SEC_MYR)-am.C2)/(log(1/(T+273.15))-am.C3)
    (1-am.beta * max(g*am.alpha+1, zero(F))^(1/am.alpha))^(1/am.beta)
end
# Fanning Curvilinear, simplified Box-Cox
function reltracklength(t::Number, T::Number, am::Union{Ketcham2007FC, Guenthner2013FC})
    g = am.C0 + am.C1*(log(t*SEC_MYR)-am.C2)/(log(1/(T+273.15))-am.C3)
    r = 1/(g^(1/am.alpha) + 1)
end
# Fanning Linear (Fanning Arrhenius), no Box-Cox
function reltracklength(t::Number, T::Number, am::Jones2021FA{F}) where {F<:AbstractFloat}
    g = am.C0 + am.C1*(log(t*SEC_MYR)-am.C2)/((1/(T+273.15))-am.C3)
    r = max(min(g, one(F)), zero(F))
end
# Parallel Curvilinear, no Box-Cox
function reltracklength(t::Number, T::Number, am::Yamada2007PC)
    r = exp(-exp(am.c0p + am.c1p*(log(t*SEC_MYR) - am.bp*log(1/(T+273.15)))))
end


"""
```julia
reltrackdensity(t, T, am::AnnealingModel)
```
Calculate the relative track density `ρ` corresponding to a given 
relative track length `r` 

Follows the relations of Ketcham et al. (2000), equations 7a and 7b 
for apatite and Tagami et al. (1999) for zircon

See also: `reltracklength`.
"""
reltrackdensity(t::Number, T::Number, am::ZirconAnnealingModel) = reltrackdensityzrn(reltracklength(t, T, am))
reltrackdensity(t::Number, T::Number, am::MonaziteAnnealingModel) = reltrackdensitymnz(reltracklength(t, T, am))
reltrackdensity(t::Number, T::Number, am::ApatiteAnnealingModel) = reltrackdensityap(reltracklength(t, T, am))
function reltrackdensityzrn(r::T) where {T<:Number}
    Tf = float(T)
    if r < 0.2
        zero(Tf)
    elseif r < 1
        Tf(1.25)*(r-Tf(0.2))
    else
        one(Tf)
    end
end
function reltrackdensitymnz(r::T) where {T<:Number}
    Tf = float(T)
    if r < 0.5
        zero(Tf)
    elseif r < 1
        Tf(2)*(r-Tf(0.5))
    else
        one(Tf)
    end
end
function reltrackdensityap(r::T) where {T<:Number}
    Tf = float(T)
    if r < 0.5274435106696789
        zero(Tf)
    elseif r < 0.765
        Tf(9.205)*r^2 - Tf(9.157)*r + Tf(2.269)
    elseif r < 1
        Tf(1.6)*r-Tf(0.6)
    else
        one(Tf)
    end
end

ellipse(x, lc) = @. sqrt(abs((1 - x^2/lc^2)*( 1.632*lc - 10.879)^2))
alrline(x, θalr) = @. (0.1035*θalr - 2.250) + x * tan(deg2rad(θalr))

"""
```julia
lcmod(l, θ)
```
Calculate the model c-axis equivalent length ("lc,mod") given a measured
"confined" fission track length `l` [microns] and angle from the c-axis 
`θ` [degrees] following the approach of Donelick et al. 1999 
(doi: 10.2138/am-1999-0902) 
"""
lcmod(x::ApatiteTrackLengthOriented) = x.lcmod
function lcmod(l, θ)
    x = l*cos(deg2rad(θ))
    y = l*sin(deg2rad(θ))
    
    fobj = curve_fit(ellipse, Float64[x], Float64[y], Float64[l])
    lcm = only(fobj.param)
    # Return if we're above the minimum length for ALR
    lcm > 12.96 && return lcm

    # Return if we're below the minimum angle for ALR
    θalr = 0.304 * exp(0.439*lcm)
    θ < θalr && return lcm
 
    # Otherwise, fit to a linear segment for ALR
    fobj = curve_fit(alrline, [x], [y], Float64[θalr])

    lcm = log(only(fobj.param)/0.304)/0.439
    return lcm
end
export lcmod

"""
```julia
rmr0model(F, Cl, OH, Mn=0, Fe=0, others=0)
```
Calculate rmr0 as a function of composition (specified in terms of
atoms per fomula unit, or APFU) for "multikinetic" apatite fission 
track thermochronology.

Implements equation 11 from Ketcham et al. 2007 
(doi: 10.2138/am.2007.2281)
```
rmr0 = (-0.0495 -0.0348F +0.3528|Cl - 1| +0.0701|OH - 1| 
        -0.8592Mn -1.2252Fe -0.1721Others)^0.1433
```
"""
function rmr0model(F, Cl, OH, Mn=0, Fe=0, others=0)
    F+Cl+OH ≈ 2 || error("F, Cl, and OH should sum to 2")
    h = (- 0.0348F + 0.3528abs(Cl - 1) + 0.0701abs(OH - 1) +
         - 0.8592Mn - 1.2252Fe - 0.1721others -0.0495      )
    return h^0.1433
end
export rmr0model

"""
```julia
rmr0fromcl(Cl)
```
Calculate `rmr0` as a function of chlorine content `Cl` [APFU] for 
"multikinetic" apatite fission track following the relation (Fig. 7a) 
of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
```
rmr0 = 1 - exp(2.107(1 - abs(Cl - 1)) - 1.834)
```
"""
rmr0fromcl(Cl) = 1 - exp(2.107*(1 - abs(Cl - 1)) - 1.834)
export rmr0fromcl

"""
```julia
rmr0fromdpar(dpar)
```
Calculate `rmr0` as a function of `dpar` for "multikinetic" apatite 
fission track following the relation (Fig. 7b) of Ketcham et al. 1999
(doi: 10.2138/am-1999-0903)
```
rmr0 = 1 - exp(0.647(dpar-1.75) - 1.834)
```
"""
rmr0fromdpar(dpar) = 1 - exp(0.647(dpar-1.75) - 1.834)
export rmr0fromdpar


# dpar = [1.5851, 1.6458, 1.7066, 1.8646, 1.9497, 1.8767, 1.8403, 2.0469, 2.1441, 2.3385, 2.4236, 2.4358, 2.9948, 4.5625, 4.9757]
# l0 = [15.7784, 16.2913, 16.0571, 15.9679, 16.2969, 16.2077, 16.2077, 16.1463, 16.1408, 16.4195, 16.3805, 16.4474, 16.3805, 16.9882, 16.9268]
# l0m = 15.63 .+ 0.283.*dpar
# nanstd(l0 - l0m)
"""
```julia
l0 = apatitel0fromdpar(dpar)
```
Calculate `l0` as a function of `dpar` for "multikinetic" apatite 
fission track following the relation (equation 1) of Carlson et al. 1999
(doi: 10.2138/am-1999-0901), that is

``l_{0,m} = 15.63 + 0.283*D_{par}``

The results may be used along with a constant uncertainty of 0.1367 based on the scatter around the 
appropriate line in Figure 1 of Carlson et al. 1999.
"""
apatitel0fromdpar(dpar) = 15.63 + 0.283*dpar
export apatitel0fromdpar

# dpar = [1.5894, 1.6402, 1.7012, 1.8333, 1.8638, 1.8841, 1.9451, 2.0569, 2.1484, 2.3415, 2.4329, 2.4431, 2.9919, 4.5772, 4.9837]
# l0 = [16.189, 16.5693, 16.4301, 16.5183, 16.3606, 16.5229, 16.5507, 16.4487, 16.4394, 16.6481, 16.6342, 16.6759, 16.6064, 17.0701, 17.0377]
# l0mod = 16.10 .+ 0.205.*dpar
"""
```julia
l0 = apatitel0modfromdpar(dpar)
```
Calculate c-axis-projected `l0` as a function of `dpar` for "multikinetic" apatite 
fission track following the relation (equation 1) of Carlson et al. 1999
(doi: 10.2138/am-1999-0901), that is

``l_{0,c,mod} = 16.10 + 0.205*D_{par}``

The results may be used along with a constant uncertainty of 0.1311 based on the scatter around the 
appropriate line in Figure 1 of Carlson et al. 1999.
"""
apatitel0modfromdpar(dpar) = 16.10 + 0.205*dpar
export apatitel0modfromdpar

## --- Track length conversion for "multikinetic" fission track

function rlr(rmr::T, rmr0::T, kappa=1.04-rmr0) where {T<:AbstractFloat}
    if rmr < rmr0
        zero(T)
    elseif rmr0 == 0
        rmr
    else
        ((rmr-rmr0)/(1-rmr0))^min(kappa, one(T))
    end
end

## --- Calculate fission track ages (accounting for decay) from normalized counts (i.e. counts/(counts-per-Ma at t0))

calc_ft(t) = (exp(λ238U*t)-1)/λ238U
calc_dftdt(t) = exp(λ238U*t)
function newton_ft_age(ftobs::T; iterations::Int=16) where {T<:Number}
    Tf = float(T)
    ftage = one(Tf)
    for _ in 1:iterations
        ∂ft∂t = calc_dftdt(ftage) # Calculate derivative
        ftage += (ftobs - calc_ft(ftage))/∂ft∂t # Move towards zero
    end
    return max(ftage, zero(Tf))
end

"""
```julia
modelage(mineral::ZirconFT, Tsteps, am::ZirconAnnealingModel)
modelage(mineral::MonaziteFT, Tsteps, am::MonaziteAnnealingModel)
modelage(mineral::ApatiteFT, Tsteps, am::ApatiteAnnealingModel)
```
Calculate the precdicted fission track age of an apatite that has experienced a given 
t-T path (specified by `mineral.tsteps` for time and `Tsteps` for temperature, at a
time resolution of `step(mineral.tsteps)`) and given annealing model parameters `am`.

Possible annealing model types and the references for the equations 
which they respetively implement include 
  `Ketcham1999FC`       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `Ketcham2007FC`       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
  `Yamada2007PC`        Parallel Curvilinear zircon model of Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002)
  `Guenthner2013FC`     Fanning Curvilinear zircon model of Guenthner et al. 2013 (doi: 10.2475/03.2013.01)
  `Jones2021FA`         Fanning Arrhenius (Fanning Linear) model adapted from Jones et al. 2021 (doi: 10.5194/gchron-3-89-2021)

"""
function modelage(zircon::ZirconFT{T}, Tsteps::AbstractVector, am::ZirconAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = zircon.agesteps::FloatRange
    @assert issorted(zircon.tsteps)
    @assert eachindex(agesteps) == eachindex(zircon.tsteps) == eachindex(Tsteps)
    ΔT = zircon.offset::T
    teq = dt = step(zircon.tsteps)
    r = reltracklength(teq, Tsteps[end], am)
    ftobs = dt * reltrackdensityzrn(r) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r = reltracklength(teq, Tsteps[i]+ΔT, am)
        ftobs += dt * reltrackdensityzrn(r) * exp(λ238U * agesteps[i])
    end
    return newton_ft_age(ftobs)
end
function modelage(monazite::MonaziteFT{T}, Tsteps::AbstractVector, am::MonaziteAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = monazite.agesteps::FloatRange
    @assert issorted(monazite.tsteps)
    @assert eachindex(agesteps) == eachindex(monazite.tsteps) == eachindex(Tsteps)
    teq = dt = step(monazite.tsteps)
    ΔT = monazite.offset::T
    r = reltracklength(teq, Tsteps[end], am)
    ftobs = dt * reltrackdensitymnz(r) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r = reltracklength(teq, Tsteps[i]+ΔT, am)
        ftobs += dt * reltrackdensitymnz(r) * exp(λ238U * agesteps[i])
    end
    return newton_ft_age(ftobs)
end
function modelage(apatite::ApatiteFT{T}, Tsteps::AbstractVector, am::ApatiteAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = apatite.agesteps::FloatRange
    @assert issorted(apatite.tsteps)
    @assert eachindex(agesteps) == eachindex(apatite.tsteps) == eachindex(Tsteps)
    rmr0 = apatite.rmr0::T
    ΔT = apatite.offset::T
    teq = dt = step(apatite.tsteps)
    r = rlr(reltracklength(teq, Tsteps[end], am), rmr0)
    ftobs = dt * reltrackdensityap(r) * exp(λ238U * agesteps[end]) 
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r = rlr(reltracklength(teq, Tsteps[i]+ΔT, am), rmr0)
        ftobs += dt * reltrackdensityap(r) * exp(λ238U * agesteps[i])
    end
    return newton_ft_age(ftobs)
end

function model_ll(mineral::FissionTrackSample, Tsteps::AbstractVector, am::AnnealingModel)
    age = modelage(mineral, Tsteps, am)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end


"""
```julia
modellength(track::ApatiteTrackLengthOriented, Tsteps, am::ApatiteAnnealingModel)
```
Calculate the predicted mean and standard deviation of the distribution of fission  
track lengths of an apatite that has experienced a given t-T path (specified by 
`track.tsteps` for time and `Tsteps` for temperature, at a time resolution of 
`step(mineral.tsteps)`) and given annealing model parameters `am`.

Possible annealing model types and the references for the equations 
which they respetively implement include 
  `Ketcham1999FC`       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `Ketcham2007FC`       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
"""
function modellength(track::ApatiteTrackLengthOriented{T}, Tsteps::AbstractVector, am::ApatiteAnnealingModel{T}; trackhist::Bool=false) where {T <: AbstractFloat}
    agesteps = track.agesteps
    tsteps = track.tsteps
    ΔT = track.offset::T
    rmr0 = track.rmr0
    r = track.r
    pr = track.pr
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps) == eachindex(pr) == eachindex(r)
    teq = dt = step(tsteps)
    r[end] = rlr(reltracklength(teq, Tsteps[end]+ΔT, am), rmr0)
    pr[end] = reltrackdensityap(r[end]) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r[i] = rlr(reltracklength(teq, Tsteps[i]+ΔT, am), rmr0)
        pr[i] = reltrackdensityap(r[i]) * exp(λ238U * agesteps[i])
    end
    r .*= track.l0 # Convert from reduced length to length
    μ, σ = nanmean(r, pr), nanstd(r, pr)
    σ = sqrt(σ^2 + track.l0_sigma^2)
    if trackhist
        binlikelihoods!(track, σ)
    end
    return μ, σ
end
function modellength(track::MonaziteTrackLength{T}, Tsteps::AbstractVector, am::MonaziteAnnealingModel{T}; trackhist::Bool=false) where {T <: AbstractFloat}
    agesteps = track.agesteps
    tsteps = track.tsteps
    ΔT = track.offset::T
    r = track.r
    pr = track.pr
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps) == eachindex(pr) == eachindex(r)
    teq = dt = step(tsteps)
    r[end] = reltracklength(teq, Tsteps[end]+ΔT, am)
    pr[end] = reltrackdensitymnz(r[end]) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r[i] = reltracklength(teq, Tsteps[i]+ΔT, am)
        pr[i] = reltrackdensitymnz(r[i]) * exp(λ238U * agesteps[i])
    end
    r .*= track.l0 # Convert from reduced length to length
    μ, σ = nanmean(r, pr), nanstd(r, pr)
    σ = sqrt(σ^2 + track.l0_sigma^2)
    if trackhist
        binlikelihoods!(track, σ)
    end
    return μ, σ
end
function modellength(track::ZirconTrackLength{T}, Tsteps::AbstractVector, am::ZirconAnnealingModel{T}; trackhist::Bool=false) where {T <: AbstractFloat}
    agesteps = track.agesteps
    tsteps = track.tsteps
    ΔT = track.offset::T
    r = track.r
    pr = track.pr
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps) == eachindex(pr) == eachindex(r)
    teq = dt = step(tsteps)
    r[end] = reltracklength(teq, Tsteps[end]+ΔT, am)
    pr[end] = reltrackdensityzrn(r[end]) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1]+ΔT, Tsteps[i]+ΔT, am) + dt
        r[i] = reltracklength(teq, Tsteps[i]+ΔT, am)
        pr[i] = reltrackdensityzrn(r[i]) * exp(λ238U * agesteps[i])
    end
    r .*= track.l0 # Convert from reduced length to length
    μ, σ = nanmean(r, pr), nanstd(r, pr)
    σ = sqrt(σ^2 + track.l0_sigma^2)
    if trackhist
        binlikelihoods!(track, σ)
    end
    return μ, σ
end

function binlikelihoods!(track::FissionTrackLength{T}, bandwidth::T) where {T<:AbstractFloat}
    fill!(track.ldist, zero(T))
    if bandwidth > 0
        kernel = Normal{T}(zero(T), bandwidth)
        @assert eachindex(track.ldist) == 1:length(track.ledges)-1
        @assert eachindex(track.ledges) == 1:length(track.ledges)
        @inbounds for i in eachindex(track.pr)
            if (track.pr[i] > 0) && (track.r[i] > 0)
                lastcdf = cdf(kernel, first(track.ledges) - track.r[i])
                for li in eachindex(track.ldist)
                    nextcdf = cdf(kernel, track.ledges[li + 1] - track.r[i])
                    track.ldist[li] += (nextcdf - lastcdf) * track.pr[i]
                    lastcdf = nextcdf
                end
            end
        end
        track.ldist ./= nansum(track.ldist)*step(track.ledges)  # Normalize
    end
    return track
end

function model_ll(track::FissionTrackLength, Tsteps::AbstractVector, am::AnnealingModel)
    l,σ = modellength(track, Tsteps, am)
    return model_ll(track, σ)
end

function model_ll(track::FissionTrackLength{T}, σ::T) where {T<:AbstractFloat}
    ll = typemin(T)
    Σpr = nansum(track.pr)
    if σ > 0 && Σpr > 0
        kernel = Normal{T}(val(track), σ)
        @inbounds for i in eachindex(track.pr)
            if (track.pr[i] > 0) && (track.r[i] > 0)
                lpr = log(track.pr[i])
                ll = logaddexp(ll, logpdf(kernel, track.r[i])+lpr)
            end
        end
        ll -= log(Σpr)
    end
    return ll
end

## --- End of File