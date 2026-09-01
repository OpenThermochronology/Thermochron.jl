## --- Vitrinite reflectance kinetic model

# note that VR predictions for rapid thermal excursions are sensitive to agesteps density in that interval

struct NielsenBasinRo{T<:AbstractFloat} <: VitriniteReflectanceModel{T}
    logA::T
    Ea::Vector{T}
    stoich::Vector{T}
    Ro0::T
end
function NielsenBasinRo(T::Type{<:AbstractFloat}=Float64;
        logA::Number = 60.9856,
        Ea::AbstractVector = [34.,36.,38.,40.,42.,44.,46.,48.,50.,52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.] .* 4184.0,
        stoich::AbstractVector = normalize_stoich([0.0185,0.0143,0.0569,0.0478,0.0497,0.0344,0.0344,0.0322,0.0282,0.0062,0.1155,0.1041,0.1023,0.0760,0.0593,0.0512,0.0477,0.0086,0.0246,0.0096]),
        Ro0::Number = 0.2104,
    )
    NielsenBasinRo(T(logA), T.(Ea), T.(stoich), T(Ro0))
end
normalize_stoich(x) = x ./ sum(x)

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
    ΔT = T(273.15) + temperatureoffset(vro, rp)
    @assert issorted(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps)
    A = exp(model.logA) / SEC_MYR
    xi = copy(model.stoich)
    @inbounds for i in eachindex(tsteps, Tsteps)
        Tavg = Tsteps[i] + ΔT
        dt = abs(step_at(tsteps, i)) * SEC_MYR   # dt is a duration -- must be positive
        for j in eachindex(xi, model.Ea)
            k = A * exp(-model.Ea[j] / (8.3143 * Tavg))
            xi[j] *= exp(-k*dt)
        end
    end
    X = sum(xi)
    return model.Ro0 * exp(3.7 * (-log(max(X, T(1e-10)))))
end

model_ll(vro::VitriniteRo, Tsteps::AbstractVector, model::VitriniteReflectanceModel, rp::RegionalParameters=RegionalParameters()) =
    norm_ll(modelage(vro, Tsteps, model, rp), stdev(vro), value(vro))

## --- End of File
