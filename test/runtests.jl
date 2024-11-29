using Thermochron
using LinearAlgebra, Statistics
using Test

@testset "Utilities" begin include("testutilities.jl") end
@testset "Fission Track" begin include("testfissiontrack.jl") end
@testset "Apatite Helium" begin include("testapatitehelium.jl") end
@testset "Zircon Helium" begin include("testzirconhelium.jl") end
@testset "Inversion" begin include("testinversion.jl") end
@testset "Integrated Examples (zrn)" begin include("testcompletezrn.jl") end
@testset "Integrated Examples (all)" begin include("testcomplete.jl") end
