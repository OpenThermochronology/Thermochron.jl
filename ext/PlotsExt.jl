module PlotsExt

    using Thermochron
    using Plots

    # Plot chronometers in general
    Plots.plot(y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), y, args...; framestyle, kwargs...)
    Plots.plot(x, y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)
    Plots.plot(x::Vector{<:Chronometer}, y::Vector{<:Chronometer}, args...; framestyle=:box, kwargs...) = plot!(plot(), x, y, args...; framestyle, kwargs...)
    for P in (Plots.Plot, Plots.Subplot)
        @eval Plots.plot!(hdl::($P), y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.val.(y); yerror=2*Thermochron.err.(y), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x, y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, x, Thermochron.val.(y); yerror=2*Thermochron.err.(y), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:Chronometer}, y; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.val.(x), y; xerror=2*Thermochron.err.(x), framestyle, kwargs...)
        @eval Plots.plot!(hdl::($P), x::Vector{<:Chronometer}, y::Vector{<:Chronometer}; framestyle=:box, kwargs...) = plot!(hdl, Thermochron.val.(x), Thermochron.val.(y); xerror=2*Thermochron.err.(x), yerror=2*Thermochron.err.(y), framestyle, kwargs...)
    end

    # Age-eU plots
    Thermochron.ageeuplot(x::Vector{<:Chronometer}, args...; framestyle=:box, xlabel="eU [ppm]", ylabel="Age [Ma]", kwargs...) = ageeuplot!(plot(), x, args...; framestyle, xlabel, ylabel, kwargs...)
    Thermochron.ageeuplot!(hdl::Plots.Plot, x::Vector{<:Chronometer}; seriestype=:scatter, mscolor=:auto, kwargs...) = plot!(hdl, Thermochron.eU.(x), Thermochron.val.(x); yerror=2*Thermochron.err.(x), seriestype, mscolor, kwargs...)
    Thermochron.ageeuplot!(hdl::Plots.Subplot, x::Vector{<:Chronometer}; seriestype=:scatter, mscolor=:auto, kwargs...) = plot!(hdl, Thermochron.eU.(x), Thermochron.val.(x); yerror=2*Thermochron.err.(x), seriestype, mscolor, kwargs...)

    # Error boxes for Ar-Ar age spectra
    Thermochron.errorbox(xc::AbstractVector, y::AbstractVector, t::BitVector=trues(length(y)); kwargs...) = errorbox!(plot(), xc, y, t; kwargs...)
    function Thermochron.errorbox!(h::Plots.Plot, xc::AbstractVector, y::AbstractVector, t::BitVector=trues(length(y)); yerror::AbstractVector=zeros(size(x)), startvalue=0, framestyle=:box, label="", kwargs...)
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