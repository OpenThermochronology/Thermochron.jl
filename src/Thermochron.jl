module Thermochron

    using LinearAlgebra
    using VectorizedStatistics
    using LoopVectorization
    using ProgressMeter: Progress, update!
    using StatGeochemBase

    include("types.jl")
    include("utilities.jl")
    include("minerals.jl")
    include("helium.jl")
    include("inversion.jl")

end
