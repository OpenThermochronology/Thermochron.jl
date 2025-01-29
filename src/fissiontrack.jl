## --- Fission track functions

function equivalenttime(t::Number, T::Number, Teq::Number, fc::Union{SimplifiedCurvilinear,FanningCurvilinear})
    exp(fc.C2 + (log(t*SEC_MYR)-fc.C2)*(log(1/(Teq+273.15))-fc.C3)/(log(1/(T+273.15))-fc.C3))/SEC_MYR
end
function equivalenttime(t::Number, T::Number, Teq::Number, am::ParallelCurvilinear)
    exp(am.bp*log(1/(Teq+273.15)) + (log(t*SEC_MYR) - am.bp*log(1/(T+273.15))))/SEC_MYR
end

"""
```julia
reltracklength(t, T, am::AnnealingModel) 
```
Calculate the relative track length `r` (equal to `l/l₀`) expected after 
isothermal heating at `T` C for `t` Myr and annealing parameters `am`. 
    
Possible annealing model types and the references for the equations 
which they respetively implement include 
  `FanningCurvilinear`      Ketcham et al. 1999 apatite (doi: 10.2138/am-1999-0903)
  `SimplifiedCurvilinear`   Ketcham et al. 2007 apatite (doi: 10.2138/am.2007.2281)
  `ParallelCurvilinear`     Yamada et al. 2005 zircon (doi: 10.1016/j.chemgeo.2006.09.002)

See also: `reltrackdensity`.
"""
function reltracklength(t::Number, T::Number, fc::FanningCurvilinear{Tf}) where {Tf}
    g = fc.C0 + fc.C1*(log(t*SEC_MYR)-fc.C2)/(log(1/(T+273.15))-fc.C3)
    (1-max(g*fc.alpha+1, zero(Tf))^(1/fc.alpha)*fc.beta)^(1/fc.beta)
end
function reltracklength(t::Number, T::Number, fc::SimplifiedCurvilinear)
    g = fc.C0 + fc.C1*(log(t*SEC_MYR)-fc.C2)/(log(1/(T+273.15))-fc.C3)
    r = 1/(g^(1/fc.alpha) + 1)
end
function reltracklength(t::Number, T::Number, am::ParallelCurvilinear)
    r = exp(-exp(am.c0p + am.c1p*(log(t*SEC_MYR) - am.bp*log(1/(T+273.15)))))
end


"""
```julia
reltrackdensity(r)
reltrackdensity(t, T, am::AnnealingModel)
```
Calculate the relative track density `ρ` corresponding to a given relative
track length `r` following the approach of Ketcham et al. 2000, 
equations 7a and 7b.

See also: `reltracklength`.
"""
function reltrackdensity(r::T) where T<:Number
    Tf = float(T)
    if r < 0.5274435106696789
        zero(Tf)
    elseif r < 0.765
        9.205*r^2 - 9.157*r + 2.269
    elseif r < 1
        1.6*r-0.6
    else
        one(Tf)
    end
end
reltrackdensity(t::Number, T::Number, am::AnnealingModel) = reltrackdensity(reltracklength(t, T, am))


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
lcmod(x::FissionTrackLength) = x.lcmod
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
    h = - 0.0348F + 0.3528abs(Cl - 1) + 0.0701abs(OH - 1) 
        - 0.8592Mn - 1.2252Fe - 0.1721others -0.0495
    return h^0.1433
end
export rmr0model

"""
```julia
rmr0fromdpar(dpar)
```
Calculate `rmr0` as a function of `dpar` for "multikinetic" apatite 
fission track following the relation (Fig. 7) of Ketcham et al. 1999
(doi: 10.2138/am-1999-0903)
```
rmr0 = 1 - exp(0.647(dpar-1.75) - 1.834)
```
"""
rmr0fromdpar(dpar) = 1 - exp(0.647(dpar-1.75) - 1.834)
export rmr0fromdpar

## --- 

function rlr(rmr::T, rmr0::T, kappa=1.04-rmr0) where {T<:AbstractFloat}
    (max(rmr-rmr0, zero(T))/(1-rmr0))^kappa
end

## ---

"""
```julia
modelage(mineral::ApatiteFT, Tsteps, am::AnnealingModel)
```
Calculate the precdicted fission track age of an apatite that has experienced a given 
t-T path (specified by `mineral.tsteps` for time and `Tsteps` for temperature, at a
time resolution of `step(mineral.tsteps)`) and given annealing model parameters `am`.

Possible annealing model types and the references for the equations 
which they respetively implement include 
  `FanningCurvilinear`      Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `SimplifiedCurvilinear`   Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
"""
function modelage(apatite::ApatiteFT{T}, Tsteps::AbstractVector, am::AnnealingModel{T}) where {T <: AbstractFloat}
    tsteps = apatite.tsteps
    rmr0 = apatite.rmr0
    @assert issorted(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps)
    teq = dt = step(tsteps)
    r = rlr(reltracklength(teq, Tsteps[end], am), rmr0)
    ftage = dt * reltrackdensity(r)
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1], Tsteps[i], am) + dt
        r = rlr(reltracklength(teq, Tsteps[i], am), rmr0)
        ftage += dt * reltrackdensity(r)
    end
    return ftage
end

function model_ll(mineral::FissionTrackSample, Tsteps::AbstractVector, am::AnnealingModel)
    age = modelage(mineral, Tsteps, am)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end


"""
```julia
modellength(track::ApatiteTrackLength, Tsteps, am::AnnealingModel)
```
Calculate the predicted mean and standard deviation of the distribution of fission  
track lengths of an apatite that has experienced a given t-T path (specified by 
`track.tsteps` for time and `Tsteps` for temperature, at a time resolution of 
`step(mineral.tsteps)`) and given annealing model parameters `am`.

Possible annealing model types and the references for the equations 
which they respetively implement include 
  `FanningCurvilinear`      Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `SimplifiedCurvilinear`   Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
"""
function modellength(track::ApatiteTrackLength{T}, Tsteps::AbstractVector, am::AnnealingModel{T}) where {T <: AbstractFloat}
    tsteps = track.tsteps
    rmr0 = track.rmr0
    r = track.r
    pr = track.pr
    @assert issorted(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == eachindex(r)
    teq = dt = step(tsteps)
    r[end] = rlr(reltracklength(teq, Tsteps[end], am), rmr0)
    pr[end] = reltrackdensity(r[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1], Tsteps[i], am) + dt
        r[i] = rlr(reltracklength(teq, Tsteps[i], am), rmr0)
        pr[i] = reltrackdensity(r[i])
    end
    return nanmean(r, pr), nanstd(r, pr)
end

function model_ll(track::ApatiteTrackLength, Tsteps::AbstractVector, am::AnnealingModel)
    l,σ = modellength(track, Tsteps, am) .* am.l0
    lc = lcmod(track)
    δ = l - lc
    σ² = σ^2 + am.l0_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end

## --- End of File