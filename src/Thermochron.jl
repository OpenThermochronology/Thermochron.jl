module Thermochron

    using LinearAlgebra
    using VectorizedStatistics
    using LoopVectorization
    using ProgressMeter: Progress, update!
    using StatGeochemBase

    include("types.jl")
    include("minerals.jl")
    include("zirconhelium.jl")
    include("inversion.jl")

end
