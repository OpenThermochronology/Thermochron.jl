## --- Check validity and consistency of temporal discretization

function checktimediscretization(::Type{T}, agesteps, tsteps=nothing) where {T<:AbstractFloat}
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps = ((first(agesteps) - step_at(agesteps,1)/2) .- agesteps))
    isnothing(agesteps) && (agesteps = (last(tsteps) + step_at(tsteps, lastindex(tsteps))/2) .- tsteps)
    @assert issorted(tsteps, lt=<=) "`tsteps` must be in strictly increasing order"
    @assert first(tsteps) >= 0 "all `tsteps` must be positive"
    @assert issorted(agesteps, lt=<=, rev=true) "`agesteps` must be in strictly decreasing order"
    @assert last(agesteps) >= 0 "all `agesteps` must be positive"
    @assert eachindex(agesteps) == eachindex(tsteps) "`tsteps` and `agesteps` must have equivalent indices"
    @assert tsteps ≈ (first(agesteps) - step_at(agesteps,1)/2) .- agesteps "`tsteps` and `agesteps` must represent the same chronology"
    return applyeltype(T, agesteps), applyeltype(T, tsteps)
end

# Ensure a specific element type
applyeltype(::Type{T}, x::AbstractArray{T}) where {T} = x
applyeltype(::Type{T}, x::AbstractArray) where {T} = T.(x)
applyeltype(::Type{T}, x::OrdinalRange) where {T} = range(T(first(x)), T(last(x)), length(x))

## --- Parse some input parameters into desired forms

function parseaftparams(;
        l0::Number = NaN,
        l0_sigma::Number = NaN,
        dpar::Number = NaN,
        F::Number = NaN,
        Cl::Number = NaN,
        OH::Number = NaN,
        rmr0::Number = NaN,
        oriented::Bool = false,
    )
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            apatite_rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(Cl)
            apatite_rmr0fromcl(Cl)
        elseif !isnan(dpar)
            apatite_rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    if isnan(dpar)
        # Estimate dpar using the relation of Ketcham et al. 1999 (Fig. 7b)
        dpar = apatite_dparfromrmr0(rmr0)
    end
    if isnan(l0)
        # Use the relation of Carlson et al. 1999 for c-axis-projected vs unoriented tracks (equation 1)
        l0 = oriented ? apatite_l0modfromdpar(dpar) : apatite_l0fromdpar(dpar)
    end
    if isnan(l0_sigma)
        # Scatter around the fit of Carlson et al. 1999 for c-axis-projected vs unoriented tracks
        l0_sigma = oriented ? 0.1311 : 0.1367
    end
    return l0, l0_sigma, rmr0
end

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
    uaftm = (haskey(model, :uaftm) ? model.uaftm : Ketcham1999FC(:unoriented))::ApatiteAnnealingModel{T}

    # Dictionaries to store reused `r` and `pr` vectors for fission track length chronometers
    # These will be indexed by hash, such that identical tracks can reuse the same `r` and `pr`
    rdict = Dict{UInt64,Vector{T}}()
    prdict = Dict{UInt64,Vector{T}}()
    calcdict = Dict{UInt64,Vector{T}}()

    # Create and fill vectors of chronometers
    chrons = Chronometer[]
    damodels = Model[]
    for i in eachindex(mineral)
        first_index = findclosest(crystage[i], agesteps)
        sample_agesteps = agesteps[first_index:end]
        mineral = lowercase(string(ds.mineral[i]))
        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0
        name = (haskey(ds, :grain_name) && ds.grain_name[i]==ds.grain_name[i]) ? string(ds.grain_name[i]) : ""
        notes = (haskey(ds, :notes) && ds.notes[i]==ds.notes[i]) ? string(ds.notes[i]) : ""
        He_step_heating_file = haskey(ds, :He_step_heating_file) ? ds.He_step_heating_file[i] : ""
        Ar_step_heating_file = haskey(ds, :Ar_step_heating_file) ? ds.Ar_step_heating_file[i] : haskey(ds, :mdd_file) ? ds.mdd_file[i] : ""

        if mineral === "zircon"
            # Zircon helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ZirconHe(T;
                    age = ds.raw_He_age_Ma[i],
                    age_sigma = ds.raw_He_age_sigma_Ma[i],
                    r = ds.halfwidth_um[i],
                    dr, offset, name, notes,
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                    agesteps = sample_agesteps,
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
                    offset, name, notes,
                    agesteps = sample_agesteps,
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end
            # Zircon fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN
                l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN
                h = hash((ZirconTrackLength, offset, sample_agesteps, l0, l0_sigma))
                haskey(rdict, h) || (rdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                haskey(prdict, h) || (prdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                haskey(calcdict, h) || (calcdict[h] = zeros(T, 2))
                c = ZirconTrackLength(T;
                    length = ds.track_length_um[i],
                    offset, l0, l0_sigma, name, notes,
                    agesteps = sample_agesteps,
                    r = rdict[h],
                    pr = prdict[h],
                    calc = calcdict[h],
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end
        elseif mineral === "monazite"
            # Monazite fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) && (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = MonaziteFT(T;
                    age = ds.FT_age_Ma[i],
                    age_sigma = ds.FT_age_sigma_Ma[i],
                    offset, name, notes,
                    agesteps = sample_agesteps,
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
            # Monazite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN
                l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN
                h = hash((MonaziteTrackLength, offset, sample_agesteps, l0, l0_sigma))
                haskey(rdict, h) || (rdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                haskey(prdict, h) || (prdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                haskey(calcdict, h) || (calcdict[h] = zeros(T, 2))
                c = MonaziteTrackLength(T;
                    length = ds.track_length_um[i],
                    offset, l0, l0_sigma, name, notes,
                    agesteps = sample_agesteps,
                    r = rdict[h],
                    pr = prdict[h],
                    calc = calcdict[h],
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
        elseif mineral === "apatite"
            # Apatite helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ApatiteHe(T;
                    age = ds.raw_He_age_Ma[i],
                    age_sigma = ds.raw_He_age_sigma_Ma[i],
                    r = ds.halfwidth_um[i],
                    dr, offset, name, notes,
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                    agesteps = sample_agesteps,
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
                    offset, name, notes,
                    dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                    F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                    Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                    OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                    rmr0 =haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                    agesteps = sample_agesteps,
                )
                push!(chrons, c)
                push!(damodels, aftm)
            end
            # Apatite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                l0 = haskey(ds, :l0_um) ? ds.l0_um[i] : NaN
                l0_sigma = haskey(ds, :l0_sigma_um) ? ds.l0_sigma_um[i] : NaN
                dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN
                F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN
                Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN
                OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN
                rmr0 = haskey(ds, :rmr0) ? ds.rmr0[i] : NaN
                if haskey(ds, :track_angle_degrees) && !isnan(ds.track_angle_degrees[i])
                    l0, l0_sigma, rmr0 = parseaftparams(;l0, l0_sigma, dpar, F, Cl, OH, rmr0, oriented=true)
                    h = hash((ApatiteTrackLengthOriented, offset, sample_agesteps, l0, l0_sigma, rmr0))
                    haskey(rdict, h) || (rdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                    haskey(prdict, h) || (prdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                    haskey(calcdict, h) || (calcdict[h] = zeros(T, 2))
                    c = ApatiteTrackLengthOriented(T;
                        length = ds.track_length_um[i],
                        angle = ds.track_angle_degrees[i],
                        offset, l0, l0_sigma, dpar, F, Cl, OH, rmr0, name, notes,
                        agesteps = sample_agesteps,
                        r = rdict[h],
                        pr = prdict[h],
                        calc = calcdict[h],
                    )
                    push!(chrons, c)
                    push!(damodels, aftm)
                else
                    l0, l0_sigma, rmr0 = parseaftparams(;l0, l0_sigma, dpar, F, Cl, OH, rmr0, oriented=false)
                    h = hash((ApatiteTrackLength, offset, sample_agesteps, l0, l0_sigma, rmr0))
                    haskey(rdict, h) || (rdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                    haskey(prdict, h) || (prdict[h] = zeros(T, lastindex(agesteps)-first_index+1))
                    haskey(calcdict, h) || (calcdict[h] = zeros(T, 2))
                    c = ApatiteTrackLength(T;
                        length = ds.track_length_um[i],
                        offset, l0, l0_sigma, dpar, F, Cl, OH, rmr0, name, notes,
                        agesteps = sample_agesteps,
                        r = rdict[h],
                        pr = prdict[h],
                        calc = calcdict[h],
                    )
                    push!(chrons, c)
                    push!(damodels, uaftm)
                end
            end
        end
        if haskey(ds, :D0_cm_2_s) && haskey(ds, :Ea_kJ_mol) && (0 < ds.D0_cm_2_s[i]) && (0 < ds.Ea_kJ_mol[i])
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : ""
            dm = Diffusivity(
                D0 = T(ds.D0_cm_2_s[i]),
                D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                Ea = T(ds.Ea_kJ_mol[i]),
                Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
            )
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                DomainType = if (geometry == "slab") || (geometry == "planar")
                    PlanarHe
                elseif (geometry === "spherical")
                    SphericalHe
                else
                    @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                    SphericalHe
                end
                c = DomainType(T;
                    age = ds.raw_He_age_Ma[i],
                    age_sigma = ds.raw_He_age_sigma_Ma[i],
                    stoppingpower = alphastoppingpower(ds.mineral[i]),
                    r = ds.halfwidth_um[i],
                    dr, offset, name, notes,
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                    agesteps = sample_agesteps,
                )
                push!(chrons, c)
                push!(damodels, dm)
            elseif haskey(ds, :raw_Ar_age_Ma) && haskey(ds, :raw_Ar_age_sigma_Ma) && (0 < ds.raw_Ar_age_sigma_Ma[i]/ds.raw_Ar_age_Ma[i])
                DomainType = if (geometry == "slab") || (geometry == "planar")
                    PlanarAr
                elseif (geometry === "spherical")
                    SphericalAr
                else
                    @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                    SphericalAr
                end
                # Planar slab argon
                c = DomainType(T;
                    age = ds.raw_Ar_age_Ma[i],
                    age_sigma = ds.raw_Ar_age_sigma_Ma[i],
                    r = ds.halfwidth_um[i],
                    dr, offset, name, notes,
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    K40_matrix = (haskey(ds, :K40_matrix_ppm) && !isnan(ds.K40_matrix_ppm[i])) ? ds.K40_matrix_ppm[i] : 0,
                    grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                    agesteps = sample_agesteps,
                )
                push!(chrons, c)
                push!(damodels, dm)
            end
        end
        if !isempty(Ar_step_heating_file) && Ar_step_heating_file==Ar_step_heating_file
            dds = importdataset(Ar_step_heating_file, importas=:Tuple)
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : ""
            DomainType = if (geometry == "slab") || (geometry == "planar")
                PlanarAr
            elseif (geometry === "spherical")
                SphericalAr
            else
                @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                SphericalAr
            end
            r = (haskey(ds, :halfwidth_um) && !isnan(ds.halfwidth_um[i])) ? ds.halfwidth_um[i] : 100
            fraction_experimental = dds.fraction_degassed
            fraction_experimental_sigma = haskey(dds, :fraction_experimental_sigma) ? dds.fraction_experimental_sigma : fill(0.005, size(fraction_experimental))
            tsteps_experimental = issorted(dds.time_s, lt=<=) ? dds.time_s : cumsum(dds.time_s)
            if haskey(dds, :lnD0_a_2) && count(!isnan, dds.lnD0_a_2) > 0
                c = MultipleDomain(T, DomainType;
                        step_age = dds.age_Ma,
                        step_age_sigma = dds.age_sigma_Ma,
                        fraction_experimental,
                        fraction_experimental_sigma,
                        tsteps_experimental,
                        Tsteps_experimental = dds.temperature_C,
                        fit = dds.fit,
                        offset, r, dr, name, notes,
                        volume_fraction = dds.volume_fraction[.!isnan.(dds.volume_fraction)],
                        K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                        K40_matrix = (haskey(ds, :K40_matrix_ppm) && !isnan(ds.K40_matrix_ppm[i])) ? ds.K40_matrix_ppm[i] : 0,
                        grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                        agesteps = sample_agesteps,
                )
                tdomains = .!isnan.(dds.lnD0_a_2)
                dm = MDDiffusivity(
                    D0 = (T.(exp.(dds.lnD0_a_2[tdomains]).*(r/10000)^2)...,),
                    D0_logsigma = (T.(haskey(dds, :lnD0_a_2_sigma) ? dds.lnD0_a_2_sigma[tdomains] : fill(log(2)/2, count(tdomains)))...,),
                    Ea = (T.(dds.Ea_kJ_mol[tdomains])...,),
                    Ea_logsigma = (T.(haskey(dds, :Ea_logsigma) ? dds.Ea_logsigma[tdomains] : fill(log(2)/4, count(tdomains)))...,),
                )
                push!(chrons, c)
                push!(damodels, dm)
            elseif haskey(ds, :D0_cm_2_s) && haskey(ds, :Ea_kJ_mol)
                c = SingleDomain(T, DomainType;
                    step_age = dds.age_Ma,
                    step_age_sigma = dds.age_sigma_Ma,
                    fraction_experimental,
                    fraction_experimental_sigma,
                    tsteps_experimental,
                    Tsteps_experimental = dds.temperature_C,
                    fit = dds.fit,
                    offset, r, dr, name, notes,
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    K40_matrix = (haskey(ds, :K40_matrix_ppm) && !isnan(ds.K40_matrix_ppm[i])) ? ds.K40_matrix_ppm[i] : 0,
                    grainsize_matrix = (haskey(ds, :grainsize_matrix_mm) && !isnan(ds.grainsize_matrix_mm[i])) ? ds.grainsize_matrix_mm[i] : 1,
                    agesteps = sample_agesteps,
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
        if !isempty(He_step_heating_file) && He_step_heating_file==He_step_heating_file
            dds = importdataset(He_step_heating_file, importas=:Tuple)
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : ""
            DomainType = if mineral === "apatite"
                ApatiteHe
            elseif mineral === "zircon"
                ZirconHe
            elseif (geometry === "slab") || (geometry === "planar")
                PlanarHe
            elseif (geometry === "spherical")
                SphericalHe
            else
                @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                SphericalHe
            end
            r = (haskey(ds, :halfwidth_um) && !isnan(ds.halfwidth_um[i])) ? ds.halfwidth_um[i] : 100
            Rstep = dds.He_4_He_3
            Rstep_sigma = dds.He_4_He_3_sigma
            Rbulk = nanmean(dds.He_4_He_3, dds.He_3)
            fraction_experimental = cumsum(dds.He_3)
            total_He_3 = last(fraction_experimental)
            fraction_experimental ./= total_He_3 # Rescale from 0-1
            fraction_experimental_sigma = haskey(dds, :He_3_sigma) ? dds.He_3_sigma./total_He_3 : fill(0.005, size(fraction_experimental))
            tsteps_experimental = issorted(dds.time_s, lt=<=) ? dds.time_s : cumsum(dds.time_s)
            if haskey(dds, :lnD0_a_2) && count(!isnan, dds.lnD0_a_2) > 0
                @warn "Multiple domain He not implemented, skipping"
            else
                c = SingleDomain(Float64, DomainType;
                    step_age = Rstep./Rbulk,
                    step_age_sigma = Rstep_sigma./Rbulk,
                    fraction_experimental,
                    fraction_experimental_sigma,
                    tsteps_experimental,
                    Tsteps_experimental = dds.temperature_C,
                    fit = dds.fit,
                    r, dr, name, notes,
                    age = ds.raw_He_age_Ma[i],
                    age_sigma = ds.raw_He_age_sigma_Ma[i],
                    U238 = ds.U238_ppm[i],
                    Th232 = ds.Th232_ppm[i],
                    Sm147 = ds.Sm147_ppm[i],
                    agesteps,
                )
                dm = if mineral === "apatite"
                    SDDiffusivity(
                        model = adm,
                        scale = 1.0 + rand()/1e6,   # Twiddle to establish uniqueness 
                        scale_logsigma = log(2),
                    )
                elseif mineral === "zircon"
                    SDDiffusivity(
                        model = zdm,
                        scale = 1.0 + rand()/1e6,   # Twiddle to establish uniqueness 
                        scale_logsigma = log(2),
                    )
                else # Custom diffusivity
                    Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                end
                push!(chrons, c)
                push!(damodels, dm)
            end
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