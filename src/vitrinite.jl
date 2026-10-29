## --- Vitrinite reflectance kinetic model

# note that VR predictions for rapid thermal excursions are sensitive to agesteps density in that interval

"""
```julia
NielsenBasinRo(;
    logA = 60.9856,                 # [log(1/Myr)] Rate constant
    Ea = (34:2:72) .* 4184.0,       # [kJ/mol] activation energy
    stoich = [0.0185,..,0.0096],    # [unitless] fraction of each component
    Ro0 = 0.2104,                   # [unitless] preexponential factor
)
```
The basin%Ro model of Nielsen et al. 2015 (doi: 10.1111/bre.12160)
"""
struct NielsenBasinRo{T<:AbstractFloat,N} <: VitriniteReflectanceModel{T}
    logA::T
    Ea::NTuple{N,T}
    stoich::NTuple{N,T}
    Ro0::T
end
function NielsenBasinRo(;
        logA::T = 60.9856,                      # [log(1/Myr)] Rate constant
        Ea::AbstractVector{T} = (34:2:72) .* 4184.0, # [kJ/mol] activation energy
        stoich::AbstractVector{T} = [0.0185,0.0143,0.0569,0.0478,0.0497,0.0344,0.0344,0.0322,0.0282,0.0062,0.1155,0.1041,0.1023,0.0760,0.0593,0.0512,0.0477,0.0086,0.0246,0.0096],  # [unitless] fraction of each component
        Ro0::T = 0.2104,                        # [unitless] preexponential factor
    ) where {T}
    @assert length(Ea) == length(stoich) "Ea and stoich must be same length"
    Ea_tuple = (T.(Ea)...,) 
    stoich_tuple = (T.(stoich ./ sum(stoich))...,) # Normalize
    NielsenBasinRo(T(logA), Ea_tuple, stoich_tuple, T(Ro0))
end

## --- Forward model

"""
```julia
modelage(vro::VitriniteRo, Tsteps, model::VitriniteReflectanceModel)
```
Calculate the predicted vitrinite reflectance (%Ro) of a sample that has
experienced a given t-T path (specified by `vro.tsteps` for time and
`Tsteps` for temperature), given kinetic model parameters `model`.
"""
function modelage(vro::VitriniteRo{T}, Tsteps::AbstractVector, model::NielsenBasinRo{T}, rp::RegionalParameters{T}=RegionalParameters{T}()) where {T<:AbstractFloat}
    tsteps = tsteps_geol(vro)
    @assert issorted(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps)
    xi = vro.xi::Vector{T}
    @assert eachindex(xi) == eachindex(model.stoich) == eachindex(model.Ea)
    logA = model.logA::T
    ΔT = T(273.15) + temperatureoffset(vro, rp)
    copyto!(xi, model.stoich)
    @inbounds for i in eachindex(tsteps)
        TK = Tsteps[i] + ΔT
        dt = abs(step_at(tsteps, i))
        for j in eachindex(xi)
            k = exp(logA - model.Ea[j]/(8.3143 * TK))
            xi[j] *= exp(-k*dt)
        end
    end
    X = sum(xi)
    return X < T(1e-10) ? model.Ro0 * exp(T(37)) : model.Ro0 * exp(3.7 * (-log(X)))
end

model_ll(vro::VitriniteRo, Tsteps::AbstractVector, model::VitriniteReflectanceModel, rp::RegionalParameters=RegionalParameters()) =
    norm_ll(modelage(vro, Tsteps, model, rp), stdev(vro), value(vro))

## --- End of File
