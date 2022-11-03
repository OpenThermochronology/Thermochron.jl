using Thermochron
using LinearAlgebra
using Test

@testset "ZrnHe.jl" begin include("testZrnHe.jl") end
@testset "inversion.jl" begin include("testinversion.jl") end
@testset "Examples" begin include("testexamples.jl") end
