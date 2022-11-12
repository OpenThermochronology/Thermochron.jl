using Thermochron
using LinearAlgebra
using Test

@testset "Zircon helium" begin include("testZrnHe.jl") end
@testset "Inversion" begin include("testinversion.jl") end
@testset "Integrated Examples" begin include("testexamples.jl") end
