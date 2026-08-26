using Test
using SMRPInversion

@testset "SMRPInversion" begin
    @testset "Deterministic (Occam) inversion" begin
        include("test_deterministic.jl")
    end
end
