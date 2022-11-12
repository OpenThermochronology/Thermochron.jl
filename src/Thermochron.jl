module Thermochron

    using LinearAlgebra
    using Statistics
    using LoopVectorization
    using ProgressMeter: @showprogress

    include("minerals.jl")
    include("ZrnHe.jl")
    include("inversion.jl")

end
