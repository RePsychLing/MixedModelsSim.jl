using DataFrames
using LinearAlgebra
using MixedModels
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

@testset "nlevstbl" begin
    subj = DataFrame(nlevstbl(:subj, 20))
    @test nrow(subj) == 20
    @test isone(ncol(subj))
    @test propertynames(subj) == [:subj]
    item = DataFrame(nlevstbl(:item, 9, :lev => ["low", "medium", "high"]))
    @test nrow(item) == 9
    @test ncol(item) == 2
    @test propertynames(item) == [:item, :lev]
    @test first(first(subj)) == "S01"
    df = crossjoin(subj, item)
    @test nrow(df) == 180
    @test ncol(df) == 3
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
    @test all(flatlowertri(l) .== [1, 2, 4, 3, 5, 6])
end


@testset "update!" begin
    @testset "LMM" begin
        fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)),
                  MixedModels.dataset(:sleepstudy))
        update!(fm1; θ=[1.0, 0.0, 1.0])

        @test all(values(first(fm1.σs)) .== fm1.σ)
        @test only(VarCorr(fm1).σρ.subj.ρ) == 0.0
    end

    @testset "GLMM" begin
        cbpp = MixedModels.dataset(:cbpp)
        gm1 = fit(MixedModel, @formula((incid/hsz) ~ 1 + period + (1|herd)), cbpp, Binomial();
                wts=cbpp.hsz, fast=true)

        β = repeat([-1.], 4)
        θ = [0.5]

        β₀ = copy(fixef(gm1))
        update!(gm1; θ=θ)
        @test all(values(first(gm1.σs)) .== θ)
        @test all(gm1.β .== β₀)

        refit!(gm1, fast=false)
        β₀ = copy(fixef(gm1))
        update!(gm1; θ=θ)
        @test all(values(first(gm1.σs)) .== θ)
        @test all(gm1.β .== β₀)
    end

end

@testset "create_re" begin
    fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)),
              MixedModels.dataset(:sleepstudy))

    σs = values(first(fm1.σs)) ./ fm1.σ
    ρ = only(first(VarCorr(fm1).σρ).ρ)
    corr = [1 ρ; ρ 1]

    subjre = create_re(σs...; corrmat=corr)

    @test all(only(fm1.λ) .≈ subjre)
    @test create_re(σs...) == diagm([σs...])

    fm1unfitted = LinearMixedModel(@formula(reaction ~ 1 + days + (1 + days|subj)),
                                   MixedModels.dataset(:sleepstudy))

    update!(fm1unfitted, subjre)

    vcu, vc = VarCorr(fm1unfitted), VarCorr(fm1)

    @test all(values(vcu.σρ.subj.σ) .≈ values(vc.σρ.subj.σ))
    @test all(values(vcu.σρ.subj.ρ) .≈ values(vc.σρ.subj.ρ))
end
