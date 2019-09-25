using MixedModelsSim, PooledArrays, Test

@testset "levels" begin
    levs = makelevels(["N","Y"], inner=2)
    @test levs isa PooledVector
    @test length(levs) == 4
    @test length(unique(levs)) == 2
    @test isa(first(levs), String)
    @test first(levs) isa String
    @test_throws MethodError makelevels("foo")
    @test first(levs) == levs[2] == "N"
    @test levs[3] == last(levs) == "Y"
    nsubj = 20
    subj = makelevels(nsubj, outer=2)
    @test length(subj) == 2*nsubj
    @test first(subj) == "S01"
    @test subj[2] == "S02"
    @test last(subj) == "S$nsubj"
end
