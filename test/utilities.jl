using DataFrames
using LinearAlgebra
using MixedModelsSim
using PooledArrays
using Tables
using Test

@testset "factorproduct" begin
	items = (item = nlevels(50, 'I'), category = repeat(["ingroup", "outgroup"], inner=25))
	df = columntable(factorproduct(items, (subj = nlevels(100),)))
	@test length(df) == 3
	@test length(first(df)) == 5000
end

@testset "levels" begin
    levs = nlevels(100, 'k')
    @test length(levs) == 100
    @test length(unique(levs)) == 100
    @test first(levs) == "k001"
end

@testset "pooled!" begin
    subject = (subj = ["S1","S2","S3","S4","S5"], age=["O","O","O","Y","Y"]);

    @test_throws ArgumentError pooled!(subject)

    df = pooled!(DataFrame(subject))

    @test df.subj isa PooledArray
    @test df.age isa PooledArray
end

@testset "flatlowertri" begin
    l = LowerTriangular([1 0 0; 2 3 0; 4 5 6])
    @test all(flatlowertri(l) .== collect(1:6))
end
