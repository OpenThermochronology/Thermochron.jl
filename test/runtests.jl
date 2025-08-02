using Thermochron
using LinearAlgebra, Statistics
using Test

const GROUP = get(ENV, "GROUP", "All")
liveplot = GROUP == "All" || GROUP == "Plot"

@testset "Utilities" begin include("testutilities.jl") end
@testset "Fission Track" begin include("testfissiontrack.jl") end
@testset "Generic Helium" begin include("testhelium.jl") end
@testset "Apatite Helium" begin include("testheliumap.jl") end
@testset "Zircon Helium" begin include("testheliumzrn.jl") end
@testset "Argon" begin include("testargon.jl") end
@testset "MDD" begin include("testmdd.jl") end
@testset "Inversion" begin include("testinversion.jl") end
@testset "Integrated Examples (zrn)" begin include("testcompletezrn.jl") end
@testset "Integrated Examples (all)" begin include("testcomplete.jl") end
