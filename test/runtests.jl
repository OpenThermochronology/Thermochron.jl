using Thermochron
using LinearAlgebra, Statistics
using StatGeochem
using Test

@testset "Apatite helium" begin include("testapatitehelium.jl") end
@testset "Zircon helium" begin include("testzirconhelium.jl") end
@testset "Inversion" begin include("testinversion.jl") end
@testset "Integrated Examples (zrn)" begin include("testcompletezrn.jl") end
@testset "Integrated Examples (all)" begin include("testcomplete.jl") end
