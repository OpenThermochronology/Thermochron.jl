## --- Check validity and consistency of temporal discretization

function checktimediscretization(::Type{T}, agesteps, tsteps=nothing) where {T<:AbstractFloat}
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps = ((first(agesteps) - step_at(agesteps,1)/2) .- agesteps))
    isnothing(agesteps) && (agesteps = (last(tsteps) + step_at(tsteps, lastindex(tsteps))/2) .- tsteps)
    @assert issorted(tsteps, lt=<=) "`tsteps` must be in strictly increasing order"
    @assert first(tsteps) >= 0 "all `tsteps` must be positive"
    @assert issorted(agesteps, lt=<=, rev=true) "`agesteps` must be in strictly decreasing order"
    @assert last(agesteps) >= 0 "all `agesteps` must be positive"
    @assert tsteps ≈ (first(agesteps) - step_at(agesteps,1)/2) .- agesteps "`tsteps` and `agesteps` must represent the same chronology"
    return applyeltype(T, agesteps), applyeltype(T, tsteps)
end

# Ensure a specific element type
applyeltype(::Type{T}, x::AbstractArray{T}) where {T} = x
applyeltype(::Type{T}, x::AbstractArray) where {T} = T.(x)
applyeltype(::Type{T}, x::OrdinalRange) where {T} = range(T(first(x)), T(last(x)), length(x))

## --- Parse imported datasets as Chronometer objects

"""
```julia
chronometers([T=Float64], data, model)
```
Construct a vector of `Chronometer` objects given a dataset `data`
and model parameters `model`.
"""
chronometers(ds, model; kwargs...) = chronometers(Float64, ds, model; kwargs...)
function chronometers(T::Type{<:AbstractFloat}, ds, model;
        zirconvolumeweighting = :cylindrical,
        apatitevolumeweighting = :cylindrical,
    )
    # Spatial discretization
    dr = haskey(model, :dr) ? model.dr : one(T)

    # Temporal discretization
    tsteps = haskey(model, :tsteps) ? model.tsteps : nothing
    agesteps = haskey(model, :agesteps) ? model.agesteps : nothing
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    haskey(ds, :mineral) || @error "dataset must contain a column labeled `mineral`"
    mineral = ds.mineral
    crystage = if haskey(ds, :crystallization_age_Ma)        
        ds.crystallization_age_Ma
    elseif haskey(ds, :crystAge) # Legacy option
        ds.crystAge
    else
        @error "dataset must contain a column labeled `crystallization age [Ma]`"
    end
    @assert eachindex(mineral) == eachindex(crystage)

    # Default damage models for each mineral
    zdm = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
    adm = (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
    zftm = (haskey(model, :zftm) ? model.zftm : Yamada2007PC())::ZirconAnnealingModel{T}
    mftm = (haskey(model, :mftm) ? model.mftm : Jones2021FA())::MonaziteAnnealingModel{T}
    aftm = (haskey(model, :aftm) ? model.aftm : Ketcham2007FC())::ApatiteAnnealingModel{T}
    uaftm = (haskey(model, :uaftm) ? model.aftm : Ketcham1999FC(:unoriented))::ApatiteAnnealingModel{T}

    chrons = Chronometer[]
    damodels = Model[]
    for i in eachindex(mineral)
        first_index = findclosest(crystage[i], agesteps)
        mineral = lowercase(string(ds.mineral[i]))

        if mineral == "zircon"
            # Zircon helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ZirconHe(T;
                    age = ds.raw_He_age_Ma[i], 
                    age_sigma = ds.raw_He_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = ds.halfwidth_um[i], 
                    dr = dr, 
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                    volumeweighting = zirconvolumeweighting,
                )
                push!(chrons, c)
                push!(damodels, zdm)
            end
            # Zircon fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) &&  (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = ZirconFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end
            # Zircon fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                c = ZirconTrackLength(T;
                    length = ds.track_length_um[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN,
                    l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end

        elseif mineral == "monazite"
            # Monazite fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) && (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = MonaziteFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
            # Monazite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                c = MonaziteTrackLength(T;
                    length = ds.track_length_um[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN,
                    l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
            
        elseif mineral == "apatite"
            # Apatite helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ApatiteHe(T;
                    age = ds.raw_He_age_Ma[i], 
                    age_sigma = ds.raw_He_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = ds.halfwidth_um[i], 
                    dr = dr, 
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                    volumeweighting = apatitevolumeweighting,
                )
                push!(chrons, c)
                push!(damodels, adm)
            end
            # Apatite fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) && (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = ApatiteFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                    F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                    Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                    OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                    rmr0 =haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, aftm)
            end
            # Apatite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                if haskey(ds, :track_angle_degrees) && !isnan(ds.track_angle_degrees[i])
                    c = ApatiteTrackLengthOriented(T;
                        length = ds.track_length_um[i], 
                        angle = ds.track_angle_degrees[i],
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN,
                        l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN,
                        dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                        F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                        Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                        OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                        rmr0 = haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                        agesteps = agesteps[first_index:end],
                    )
                    push!(chrons, c)
                    push!(damodels, aftm)
                else
                    c = ApatiteTrackLength(T;
                        length = ds.track_length_um[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN,
                        l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN,
                        dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                        F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                        Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                        OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                        rmr0 = haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                        agesteps = agesteps[first_index:end],
                    )
                    push!(chrons, c)
                    push!(damodels, uaftm)
                end
            end
        elseif haskey(ds, :D0_cm_2_s) && haskey(ds, :Ea_kJ_mol) && (0 < ds.D0_cm_2_s[i]) && (0 < ds.Ea_kJ_mol[i])
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : "spherical"
            if (geometry == "slab") || (geometry == "planar")
                # Planar slab helium
                if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                    c = PlanarHe(T;
                        age = ds.raw_He_age_Ma[i], 
                        age_sigma = ds.raw_He_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        stoppingpower = alphastoppingpower(ds.mineral[i]),
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                        Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                        Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                        U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                        Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                        Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
                # Planar slab argon
                if haskey(ds, :raw_Ar_age_Ma) && haskey(ds, :raw_Ar_age_sigma_Ma) && (0 < ds.raw_Ar_age_sigma_Ma[i]/ds.raw_Ar_age_Ma[i])
                    c = PlanarAr(T;
                        age = ds.raw_Ar_age_Ma[i], 
                        age_sigma = ds.raw_Ar_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
            else
                (geometry === "sphere") || (geometry === "spherical") || @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                # Spherical helium
                if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                    c = SphericalHe(T;
                        age = ds.raw_He_age_Ma[i], 
                        age_sigma = ds.raw_He_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        stoppingpower = alphastoppingpower(ds.mineral[i]),
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                        Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                        Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                        U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                        Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                        Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
                # Spherical argon
                if haskey(ds, :raw_Ar_age_Ma) && haskey(ds, :raw_Ar_age_sigma_Ma) && (0 < ds.raw_Ar_age_sigma_Ma[i]/ds.raw_Ar_age_Ma[i])
                    c = SphericalAr(T;
                        age = ds.raw_Ar_age_Ma[i], 
                        age_sigma = ds.raw_Ar_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
            end
        elseif haskey(ds, :mdd_file) && !isempty(ds.mdd_file[i])
            mdds = importdataset(ds.mdd_file[i], importas=:Tuple)
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : "spherical"
            r = (haskey(ds, :halfwidth_um) && !isnan(ds.halfwidth_um[i])) ? ds.halfwidth_um[i] : 100
            if (geometry == "slab") || (geometry == "planar")
                c = MultipleDomain(T, PlanarAr;
                    age = mdds.age_Ma,
                    age_sigma = mdds.age_sigma_Ma,
                    fraction_experimental = mdds.fraction_degassed,
                    tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
                    Tsteps_experimental = mdds.temperature_C,
                    fit = mdds.fit,
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = r,
                    dr = dr, 
                    volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    agesteps = agesteps[first_index:end],
                )
            else
                (geometry === "spherical") || @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                c = MultipleDomain(T, SphericalAr;
                    age = mdds.age_Ma,
                    age_sigma = mdds.age_sigma_Ma,
                    fraction_experimental = mdds.fraction_degassed,
                    tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
                    Tsteps_experimental = mdds.temperature_C,
                    fit = mdds.fit,
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = r,
                    dr = dr, 
                    volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    agesteps = agesteps[first_index:end],
                )
            end
            tdomains = .!isnan.(mdds.lnD0_a_2)
            dm = MDDiffusivity(
                D0 = (T.(exp.(mdds.lnD0_a_2[tdomains]).*(r/10000)^2)...,),
                D0_logsigma = (T.(haskey(mdds, :lnD0_a_2_sigma) ? mdds.lnD0_a_2_sigma[tdomains] : fill(log(2)/2, count(tdomains)))...,),
                Ea = (T.(mdds.Ea_kJ_mol[tdomains])...,),
                Ea_logsigma = (T.(haskey(mdds, :Ea_logsigma) ? mdds.Ea_logsigma[tdomains] : fill(log(2)/4, count(tdomains)))...,),
            )
            push!(chrons, c)
            push!(damodels, dm)
        end
    end

    isempty(chrons) && @error "No chronometers found"
    return unionize(chrons), unionize(damodels)
end

## --- Other sample processing related to parsed Chronometer objects

"""
```julia
function empiricaluncertainty!(σcalc, chrons, C::Type{<:HeliumSample}; 
    fraction::Number = 1/sqrt(2), 
    sigma_eU::Number = (C<:ZirconHe) ? 100.0 : (C<:ApatiteHe) ? 10.0 : 25.0,
    sigma_offset::Number = 10.0,
)
```
Given a vector of chronometers `chrons`, update the uncertainties in `σcalc`
to reflect the uncertainty implied by the observed (empirical) scatter in the
helium ages for samples with similar eU and temperature offset (i.e., elevation).

By default, half the empirical variance (`1/sqrt(2)` of the standard deviation)
is interpreted as external uncertainty which should be reflected in `σcalc`.

Similarity in eU and offset is assessed on the basis of Gaussian kernels with
bandwidth equal to `sigma_eU` for eU and `sigma_offset` for temperature offset.
"""
function empiricaluncertainty!(σcalc::AbstractVector{T}, chrons::AbstractArray{<:Chronometer{T}}, ::Type{C}; 
        fraction::Number = 1/sqrt(2), 
        sigma_eU::Number = (C<:ZirconHe) ? 100.0 : (C<:ApatiteHe) ? 10.0 : 25.0,
        sigma_offset::Number = 10.0,
    ) where {T<:AbstractFloat, C<:HeliumSample}
    @assert eachindex(σcalc) == eachindex(chrons)
    @assert 0 <= fraction <= 1
    inds = findall(cᵢ->cᵢ isa C, chrons)
    chrons_C = chrons[inds]
    eU_C = eU.(chrons_C)
    ages_C = value.(chrons_C)
    ΔT_C = temperatureoffset.(chrons_C)
    for i ∈ inds
        nearest_eU = minimum(j->abs(eU(chrons[j]) - eU(chrons[i])), setdiff(inds, i))
        eU_kernel = Normal(eU(chrons[i]), max(sigma_eU, nearest_eU/2))
        ΔT_kernel = Normal(temperatureoffset(chrons[i]), sigma_offset)
        W = pdf.(eU_kernel, eU_C) .* pdf.(ΔT_kernel, ΔT_C)
        # Assume some fraction of weighted variance is from unknown external uncertainty
        σₑ = nanstd(ages_C, W) * fraction   # External uncertainty (est)
        σcalc[i] = max(σₑ, σcalc[i])
    end
    return σcalc
end

## --- End of File