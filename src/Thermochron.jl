module Thermochron

    using LinearAlgebra
    using Statistics
    using LoopVectorization
    using ProgressMeter: @showprogress

    include("minerals.jl")
    include("zirconhelium.jl")
    include("inversion.jl")

end
