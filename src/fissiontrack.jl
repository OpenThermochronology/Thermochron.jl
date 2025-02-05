## --- Fission track functions

function equivalenttime(t::Number, T::Number, Teq::Number, fc::FanningCurvilinearApatite)
    exp(fc.C2 + (log(t*SEC_MYR)-fc.C2)*(log(1/(Teq+273.15))-fc.C3)/(log(1/(T+273.15))-fc.C3))/SEC_MYR
end
function equivalenttime(t::Number, T::Number, Teq::Number, am::Yamada2007PC)
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
  `Ketcham1999FC`       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `Ketcham2007FC`       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)
  `Yamada2007PC`        Parallel Curvilinear zircon model of Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002)

See also: `reltrackdensity`.
"""
function reltracklength(t::Number, T::Number, fc::Ketcham1999FC{Tf}) where {Tf}
    g = fc.C0 + fc.C1*(log(t*SEC_MYR)-fc.C2)/(log(1/(T+273.15))-fc.C3)
    (1-max(g*fc.alpha+1, zero(Tf))^(1/fc.alpha)*fc.beta)^(1/fc.beta)
end
function reltracklength(t::Number, T::Number, fc::Ketcham2007FC)
    g = fc.C0 + fc.C1*(log(t*SEC_MYR)-fc.C2)/(log(1/(T+273.15))-fc.C3)
    r = 1/(g^(1/fc.alpha) + 1)
end
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
reltrackdensity(t::Number, T::Number, am::ApatiteAnnealingModel) = reltrackdensityap(reltracklength(t, T, am))
reltrackdensity(t::Number, T::Number, am::ZirconAnnealingModel) = reltrackdensityzrn(reltracklength(t, T, am))
function reltrackdensityap(r::T) where T<:Number
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
function reltrackdensityzrn(r::T) where T<:Number
    Tf = float(T)
    if r < 0.2
        zero(Tf)
    elseif r < 1
        Tf(1.25)*(r-Tf(0.2))
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
modelage(mineral::ZirconFT, Tsteps, am::ZirconAnnealingModel)
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
"""
function modelage(zircon::ZirconFT{T}, Tsteps::AbstractVector, am::ZirconAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = zircon.agesteps
    tsteps = zircon.tsteps
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps)
    teq = dt = step(tsteps)
    r = reltracklength(teq, Tsteps[end], am)
    ΣUw = Uw = exp(λ238U * agesteps[end]) # To correct for U decay
    ftage = dt * reltrackdensityzrn(r) * Uw
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1], Tsteps[i], am) + dt
        r = reltracklength(teq, Tsteps[i], am)
        ΣUw += Uw = exp(λ238U * agesteps[i])
        ftage += dt * reltrackdensityzrn(r) * Uw
    end
    μUw = ΣUw/length(tsteps)
    return ftage/μUw
end
function modelage(apatite::ApatiteFT{T}, Tsteps::AbstractVector, am::ApatiteAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = apatite.agesteps
    tsteps = apatite.tsteps
    rmr0 = apatite.rmr0
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps)
    teq = dt = step(tsteps)
    r = rlr(reltracklength(teq, Tsteps[end], am), rmr0)
    ΣUw = Uw = exp(λ238U * agesteps[end]) # To correct for U decay
    ftage = dt * reltrackdensityap(r) * Uw
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1], Tsteps[i], am) + dt
        r = rlr(reltracklength(teq, Tsteps[i], am), rmr0)
        ΣUw += Uw = exp(λ238U * agesteps[i])
        ftage += dt * reltrackdensityap(r) * Uw
    end
    μUw = ΣUw/length(tsteps)
    return ftage/μUw
end

function model_ll(mineral::FissionTrackSample, Tsteps::AbstractVector, am::AnnealingModel)
    age = modelage(mineral, Tsteps, am)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end


"""
```julia
modellength(track::ApatiteTrackLength, Tsteps, am::ApatiteAnnealingModel)
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
function modellength(track::ApatiteTrackLength{T}, Tsteps::AbstractVector, am::ApatiteAnnealingModel{T}) where {T <: AbstractFloat}
    agesteps = track.agesteps
    tsteps = track.tsteps
    rmr0 = track.rmr0
    r = track.r
    pr = track.pr
    @assert issorted(tsteps)
    @assert eachindex(agesteps) == eachindex(tsteps) == eachindex(Tsteps) == eachindex(pr) == eachindex(r)
    teq = dt = step(tsteps)
    r[end] = rlr(reltracklength(teq, Tsteps[end], am), rmr0)
    pr[end] = reltrackdensityap(r[end]) * exp(λ238U * agesteps[end])
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)),1)
        teq = equivalenttime(teq, Tsteps[i+1], Tsteps[i], am) + dt
        r[i] = rlr(reltracklength(teq, Tsteps[i], am), rmr0)
        pr[i] = reltrackdensityap(r[i]) * exp(λ238U * agesteps[i])
    end
    r .*= am.l0 # Convert from reduced length to length
    μ, σ = nanmean(r, pr), nanstd(r, pr)
    h = (4*σ^5/(3 * sum(pr)))^(1/5) # Silverman's rule for kernel bandwidth
    binlikelihoods!(track, h)
    return μ, σ
end

function binlikelihoods!(track::ApatiteTrackLength{T}, bandwidth::T) where {T<:AbstractFloat}
    if bandwidth > 0
        kernel = Normal(zero(T), bandwidth)
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
    else
        fill!(track.ldist, zero(T))
    end
    return track
end

function model_ll(track::ApatiteTrackLength, Tsteps::AbstractVector, am::ApatiteAnnealingModel)
    l,σ = modellength(track, Tsteps, am)
    return model_ll(track)
end

function model_ll(track::ApatiteTrackLength{T}) where {T}
    lc = lcmod(track)
    σ = nanstd(track.r, track.pr)
    h = (4*σ^5/(3 * sum(track.pr)))^(1/5) # Silverman's rule
    ll = typemin(T)
    if h > 0
        kernel = Normal(lc, h)
        @inbounds for i in eachindex(track.pr)
            if (track.pr[i] > 0) && (track.r[i] > 0)
                lpr = log(track.pr[i])
                ll = logaddexp(ll, logpdf(kernel, track.r[i])+lpr)
            end
        end
    end
    return ll - log(nansum(track.pr))
end

## --- End of File