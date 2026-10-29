using Test
using SMRPInversion

@testset "SMRPInversion" begin
    @testset "Deterministic (Occam) inversion" begin
        include("test_deterministic.jl")
    end
    @testset "Bootstrap inversion" begin
        include("test_bootstrap.jl")
    end
    @testset "Hierarchical noise (σᵢ² = s·σ̂ᵢ²)" begin
        include("test_hierarchical_noise.jl")
    end
end
