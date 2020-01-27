using MixedModelsSim
using Test

@testset "levels" begin
    levs = nlevels(100, 'k')
    @test length(levs) == 100
    @test length(unique(levs)) == 100
    @test first(levs) == "k001"
end
