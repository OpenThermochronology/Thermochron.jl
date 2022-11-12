using Thermochron
using LinearAlgebra
using Test

@testset "Zircon helium" begin include("testzirconhelium.jl") end
@testset "Inversion" begin include("testinversion.jl") end
@testset "Integrated Examples" begin include("testexamples.jl") end
