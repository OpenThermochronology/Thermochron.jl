module PlotsExt

    using Thermochron
    using Plots

    # Plot chronometers in general
    Plots.plot(y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), y, args...; framestyle, kwargs...)
    Plots.plot(x, y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)
    Plots.plot(x::Vector{<:Chronometer}, y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.value.(y); yerror=2*Thermochron.stdev.(y), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x, y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, x, Thermochron.value.(y); yerror=2*Thermochron.stdev.(y), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:Chronometer}, y; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.value.(x), y; xerror=2*Thermochron.stdev.(x), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:Chronometer}, y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.value.(x), Thermochron.value.(y); xerror=2*Thermochron.stdev.(x), yerror=2*Thermochron.stdev.(y), framestyle, kwargs...)
    end

    # Age-eU plots
    Thermochron.ageeuplot(x::Vector{<:Chronometer}, args...; framestyle=:box, xlabel="eU [ppm]", ylabel="Age [Ma]", kwargs...) = ageeuplot!(plot(), x, args...; framestyle, xlabel, ylabel, kwargs...)
    Thermochron.ageeuplot!(hdl::Plots.Plot, x::Vector{<:Chronometer}; seriestype=:scatter, mscolor=:auto, kwargs...) = plot!(hdl, Thermochron.eU.(x), Thermochron.value.(x); yerror=2*Thermochron.stdev.(x), seriestype, mscolor, kwargs...)
    Thermochron.ageeuplot!(hdl::Plots.Subplot, x::Vector{<:Chronometer}; seriestype=:scatter, mscolor=:auto, kwargs...) = plot!(hdl, Thermochron.eU.(x), Thermochron.value.(x); yerror=2*Thermochron.stdev.(x), seriestype, mscolor, kwargs...)

    # Error boxes for Ar-Ar age spectra
    Thermochron.errorbox(xc::AbstractVector, y::AbstractVector, t::BitVector=trues(length(y)); kwargs...) = errorbox!(plot(), xc, y, t; kwargs...)
    function Thermochron.errorbox!(h::Union{Plots.Plot, Plots.Subplot}, xc::AbstractVector, y::AbstractVector, t::BitVector=trues(length(y)); yerror::AbstractVector=zeros(size(x)), startvalue=0, framestyle=:box, label="", kwargs...)
        @assert eachindex(y) == eachindex(yerror) == eachindex(t)
        @assert (eachindex(xc) == eachindex(y)) || eachindex(xc)==firstindex(y):lastindex(y)+1
        length(xc) == length(y) && (xc = [startvalue; xc])
        xl, yl = zeros(5), zeros(5)
        labelled = false
        for i in eachindex(y)
            if t[i]
                xl .= (xc[i], xc[i], xc[i+1], xc[i+1], xc[i])
                yl .= (y[i]-yerror[i], y[i]+yerror[i], y[i]+yerror[i], y[i]-yerror[i], y[i]-yerror[i])
                s = Shape(xl, yl)
                plot!(h, s; framestyle, label=(labelled ? "" : label), kwargs...)
                labelled = true
            end
        end
        return h
    end
    Thermochron.errorbox(c::MultipleDomain; kwargs...) = errorbox!(plot(), c; kwargs...)
    function Thermochron.errorbox!(h::Union{Plots.Plot, Plots.Subplot}, c::MultipleDomain; fillalpha=0.5, excludedalpha=0.15, color=:black, kwargs...)
        errorbox!(h, c.fraction_experimental, c.age, c.fit;
            yerror = 2*c.age_sigma,
            label = "Data (2σ analytical)",
            color, 
            fillalpha,
            kwargs...
        )
        errorbox!(h, c.fraction_experimental, c.age;
            yerror = 2*c.age_sigma,
            label = "Data (excluded)",
            color,
            fillalpha,
            alpha = excludedalpha,
            kwargs...
        )
    end

    # Diffusivity and annealing models
    function Plots.plot(d::ZRDAAM, dms::Vector{<:Thermochron.Model}; framestyle=:box, layout=(2,1), size=(600,800), alpha=0.75, kwargs...) 
        hd = plot(; framestyle, xlabel="Log10(D₀ [cm^2/s])", ylabel="Probability Density", kwargs...)
        D0 = dms .|> x->log10(x.DzD0)
        histogram!(hd, D0; normalized=true, lw=0, color=lines[1], label="", bins=(minimum(D0)-0.1):0.1:(maximum(D0)+0.1), alpha, kwargs...)
        D0₀ = Normal(log10(d.DzD0), d.DzD0_logsigma/log(10))
        x = range(mean(D0₀)-3std(D0₀), mean(D0₀)+3std(D0₀), length=100)
        plot!(hd, x, pdf.(D0₀,x), color=lines[1], label="Crystalline zircon", kwargs...)
        D0 = dms .|> x->log10(x.DN17D0)
        histogram!(hd, D0; normalized=true, lw=0, color=lines[2], label="", bins=(minimum(D0)-0.1):0.1:(maximum(D0)+0.1), alpha, kwargs...)
        D0₀ = Normal(log10(d.DN17D0), d.DN17D0_logsigma/log(10))
        x = range(mean(D0₀)-3std(D0₀), mean(D0₀)+3std(D0₀), length=100)
        plot!(hd, x, pdf.(D0₀,x), color=lines[2], label="Amorphous zircon", kwargs...)

        he = plot(;framestyle, xlabel="Log10(Eₐ [kj/mol])", ylabel="Probability Density", kwargs...)
        Ea = dms .|> x->log10(x.DzEa)
        histogram!(he, Ea; normalized=true, lw=0, color=lines[1], label="posterior", bins=(minimum(Ea)-0.01):0.01:(maximum(Ea)+0.01), alpha, kwargs...)
        Ea₀ = Normal(log10(d.DzEa), d.DzEa_logsigma/log(10))
        x = range(mean(Ea₀)-3std(Ea₀), mean(Ea₀)+3std(Ea₀), length=100)
        plot!(he, x, pdf.(Ea₀,x), color=lines[1], label="prior", kwargs...)
        Ea = dms .|> x->log10(x.DN17Ea)
        histogram!(he, Ea; normalized=true, lw=0, color=lines[2], label="", bins=(minimum(Ea)-0.01):0.01:(maximum(Ea)+0.01), alpha, kwargs...)
        Ea₀ = Normal(log10(d.DN17Ea), d.DN17Ea_logsigma/log(10))
        x = range(mean(Ea₀)-3std(Ea₀), mean(Ea₀)+3std(Ea₀), length=100)
        plot!(he, x, pdf.(Ea₀,x), color=lines[2], label="", kwargs...)

        return plot(hd, he; layout, size, kwargs...)
    end
    function Plots.plot(d::RDAAM, dms::Vector{<:Thermochron.Model}; framestyle=:box, layout=(2,1), size=(600,800), alpha=0.75, kwargs...) 
        hd = plot(; framestyle, xlabel="Log10(D₀ [cm^2/s])", ylabel="Probability Density", kwargs...)
        D0 = dms .|> x->log10(x.D0L)
        histogram!(hd, D0; normalized=true, lw=0, color=lines[1], label="D0L posterior", bins=(minimum(D0)-0.1):0.1:(maximum(D0)+0.1), alpha, kwargs...)
        D0₀ = Normal(log10(d.D0L), d.D0L_logsigma/log(10))
        x = range(mean(D0₀)-3std(D0₀), mean(D0₀)+3std(D0₀), length=100)
        plot!(hd, x, pdf.(D0₀,x), color=lines[1], label="D0L prior", kwargs...)

        he = plot(;framestyle, xlabel="Log10(Eₐ [kj/mol])", ylabel="Probability Density", kwargs...)
        Ea = dms .|> x->log10(x.EaL)
        histogram!(he, Ea; normalized=true, lw=0, color=lines[1], label="EaL posterior", bins=(minimum(Ea)-0.02):0.02:(maximum(Ea)+0.02), alpha, kwargs...)
        Ea₀ = Normal(log10(d.EaL), d.EaL_logsigma/log(10))
        x = range(mean(Ea₀)-3std(Ea₀), mean(Ea₀)+3std(Ea₀), length=100)
        plot!(he, x, pdf.(Ea₀,x), color=lines[1], label="EaL prior", kwargs...)
        Ea = dms .|> x->log10(x.EaTrap)
        histogram!(he, Ea; normalized=true, lw=0, color=lines[2], label="EaTrap posterior", bins=(minimum(Ea)-0.02):0.02:(maximum(Ea)+0.02), alpha, kwargs...)
        Ea₀ = Normal(log10(d.EaTrap), d.EaTrap_logsigma/log(10))
        x = range(mean(Ea₀)-3std(Ea₀), mean(Ea₀)+3std(Ea₀), length=100)
        plot!(he, x, pdf.(Ea₀,x), color=lines[2], label="EaTrap prior", kwargs...)

        return plot(hd, he; layout, size, kwargs...)
    end
    function Plots.plot(d::MDDiffusivity, dms::Vector{<:Thermochron.Model}, r=100.0; framestyle=:box, layout=(2,1), size=(600,800), alpha=0.75, kwargs...) 
        ndomains = length(d.D0)
        hd = plot(; framestyle, xlabel="Log10(D₀/a² [1/s])", ylabel="Probability Density", kwargs...)
        for j in 1:ndomains
            D0a2 = dms .|> x->log10(x.D0[j]/(r/10000)^2)
            histogram!(hd, D0a2; normalized=true, lw=0, color=lines[j], label="", bins=(minimum(D0a2)-0.05):0.1:(maximum(D0a2)+0.1), alpha, kwargs...)
            D0a2₀ = Normal(log10(d.D0[j]./(r/10000)^2), d.D0_logsigma[j]/log(10))
            x = range(mean(D0a2₀)-3std(D0a2₀), mean(D0a2₀)+3std(D0a2₀), length=100)
            plot!(hd, x, pdf.(D0a2₀,x), color=lines[j], label="domain $j", kwargs...)
        end

        he = plot(; framestyle, xlabel="Log10(Eₐ [kj/mol])", ylabel="Probability Density")
        for j in 1:ndomains
            Ea = dms .|> x->log10(x.Ea[j])
            histogram!(he, Ea; normalized=true, lw=0, color=lines[j], label=(j==1 ? "posterior" : ""), bins=(minimum(Ea)-0.0025):0.005:(maximum(Ea)+0.005), alpha, kwargs...)
            Ea₀ = Normal(log10(d.Ea[j]), d.Ea_logsigma[j]/log(10))
            x = range(mean(Ea₀)-3std(Ea₀), mean(Ea₀)+3std(Ea₀), length=100)
            plot!(he, x, pdf.(Ea₀,x), color=lines[j], label=(j==1 ? "prior" : ""), kwargs...)
        end
        
        return plot(hd, he; layout, size, kwargs...)
    end

    # Constraint boxes
    lowerbound(x::Uniform) = x.a
    lowerbound(x::Distribution) = quantile(x, 0.025)
    upperbound(x::Uniform) = x.b
    upperbound(x::Distribution) = quantile(x, 0.975)

    Plots.plot(c::Constraint; framestyle=:box, kwargs...) = plot!(plot(), c; framestyle, kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval function Plots.plot!(hdl::($P), c::Constraint; framestyle=:box, lw=2, fillalpha=0.1, color=:black, label="Constraints", kwargs...)
            for i in eachindex(c.agedist, c.Tdist)
                t, T = c.agedist[i], c.Tdist[i]
                x = [lowerbound(t), upperbound(t), upperbound(t), lowerbound(t)]
                y = [upperbound(T), upperbound(T), lowerbound(T), lowerbound(T)]
                plot!(hdl, Shape(x, y); framestyle, lw, fillalpha, color, label, kwargs...)
                label = ""
            end
            return hdl
        end
    end


end # module