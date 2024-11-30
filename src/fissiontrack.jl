## --- Fission track models

l₀ = 16.38 # [um] initial track length
σl₀ = 0.09 # [um] initial track length uncertainty

# "Simultaneous fit" Fanning Curvilinear models 
const FCKetcham1999 = FanningCurvilinear(-19.84402202, 0.3895104539, -51.25312954, -7.642358713, -0.12327, -11.988)
export FCKetcham1999

const FCKetcham2007 = SimplifiedCurvilinear(0.39528, 0.01073, -65.12969, -7.91715, 0.04672)
export FCKetcham2007

function equivalenttime(t::Number, T::Number, Teq::Number, fc::Union{SimplifiedCurvilinear,FanningCurvilinear})
    exp(fc.C2 + (log(t*SEC_MYR)-fc.C2)*(log(1/(Teq+273.15))-fc.C3)/(log(1/(T+273.15))-fc.C3))/SEC_MYR
end

"""
```julia
reltracklength(t, T, am::AnnealingModel) 
```
Calculate the relative track length `r` (equal to `l/l₀`) expected after 
isothermal heating at `T` C for `t` Myr and annealing parameters `am`. 
    
Possible annealing model types and the references for the equations 
which they respetively implement include 
  `FanningCurvilinear`      Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)
  `SimplifiedCurvilinear`   Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)

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
reltrackdensity(t::Number, T::Number, am::AnnealingModel) = reltrackdensity(reltracklength(t, T, am))
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


ellipse(x, lc) = @. sqrt(abs((1 - x^2/lc^2)*( 1.632*lc - 10.879)^2))
alrline(x, θalr) = @. (0.1035*θalr - 2.250) + x * tan(deg2rad(θalr))

"""
```julia
lcmodel(l, θ)
```
Calculate the model c-axis equivalent length ("lc,mod") given a measured
"confined" fission track length `l` [microns] and angle from the c-axis 
`θ` [degrees] following the approach of Donelick et al. 1999 
(doi: 10.2138/am-1999-0902) 
"""
function lcmodel(l, θ)
    x = l*cos(deg2rad(θ))
    y = l*sin(deg2rad(θ))
    
    fobj = curve_fit(ellipse, Float64[x], Float64[y], Float64[l])
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
export lcmodel

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
    sum((F, Cl, OH)) ≈ 2 || error("F, Cl, and OH should sum to 2")
    h = - 0.0348F + 0.3528abs(Cl - 1) + 0.0701abs(OH - 1) 
        - 0.8592Mn - 1.2252Fe - 0.1721others -0.0495
    return h^0.1433
end
export rmr0model

## --- 

function rlr(rmr::T, rmr0::T, kappa=1.04-rmr0) where {T<:AbstractFloat}
    (max(rmr-rmr0, zero(T))/(1-rmr0))^kappa
end

## ---

function modelage(apatite::ApatiteFT{Tf}, Tsteps, am::AnnealingModel{Tf}) where {Tf <: AbstractFloat}
    tsteps = apatite.tsteps
    @assert issorted(tsteps)
    agesteps = apatite.agesteps
    @assert eachindex(tsteps) == eachindex(agesteps) == eachindex(Tsteps)
    Teq = mean(Tsteps)
    teq = ftage = zero(Tf)
    @inbounds for i in Iterators.drop(reverse(eachindex(Tsteps)), 1)
        teq += equivalenttime(dt, Tsteps[i], Teq, am)
        ftage += agesteps[i] * reltrackdensity(teq, Teq, am)
    end
    return ftage
end

## --- End of File