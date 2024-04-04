using DataFrames
using LinearAlgebra
using MixedModels
using MixedModelsSim
using PooledArrays
using Tables
using Test

@testset "factorproduct" begin
    items = (item=nlevels(50, 'I'), category=repeat(["ingroup", "outgroup"]; inner=25))
    df = columntable(factorproduct(items, (subj=nlevels(100),)))
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
    @test_throws ArgumentError nlevstbl(:baditem, 9, :nlev => ["low", "high"])
    twofac = nlevstbl(:item, 18, :level => ["low", "medium", "high"], :blur => ["N", "Y"])
    @test length(twofac) == 3
    @test length(first(twofac)) == 18
    @test_throws ArgumentError nlevstbl(:item, 18, :level => ["low", "high"],
                                        :blur => ["N", "Y"])

    item = nlevstbl(:subj, 128, :trt => ["C", "T"])
    @test item.subj isa PooledVector{String,Int16,Vector{Int16}}
    @test item.trt isa PooledVector{String,Int8,Vector{Int8}}
    @test length(item.subj) == length(item.trt)
end

@testset "pooled!" begin
    subject = (subj=["S1", "S2", "S3", "S4", "S5"], age=["O", "O", "O", "Y", "Y"])

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
        fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days | subj)),
                  MixedModels.dataset(:sleepstudy); progress=false)
        update!(fm1; θ=[1.0, 0.0, 1.0])

        @test all(values(first(fm1.σs)) .== fm1.σ)
        @test only(VarCorr(fm1).σρ.subj.ρ) == 0.0
    end

    @testset "GLMM" begin
        cbpp = MixedModels.dataset(:cbpp)
        gm1 = fit(MixedModel, @formula((incid / hsz) ~ 1 + period + (1 | herd)), cbpp,
                  Binomial();
                  wts=cbpp.hsz, fast=true, progress=false)

        β = repeat([-1.0], 4)
        θ = [0.5]

        β₀ = copy(fixef(gm1))
        update!(gm1; θ=θ)
        @test all(values(first(gm1.σs)) .== θ)
        @test all(gm1.β .== β₀)

        refit!(gm1; fast=false, progress=false)
        β₀ = copy(fixef(gm1))
        update!(gm1; θ=θ)
        @test all(values(first(gm1.σs)) .== θ)
        @test all(gm1.β .== β₀)
    end
end

@testset "create_re|θ" begin
    fm1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days | subj)),
              MixedModels.dataset(:sleepstudy); progress=false)

    σs = values(first(fm1.σs)) ./ fm1.σ
    ρ = only(first(VarCorr(fm1).σρ).ρ)
    corr = [1 ρ; ρ 1]

    subjre = create_re(σs...; corrmat=corr)

    @test all(only(fm1.λ) .≈ subjre)
    @test create_re(σs...) == diagm([σs...])

    fm1unfitted = LinearMixedModel(@formula(reaction ~ 1 + days + (1 + days | subj)),
                                   MixedModels.dataset(:sleepstudy))

    @test_logs((:warn,
                "Specifying the random effects by position instead of name is deprecated"),
               update!(fm1unfitted, subjre))

    update!(fm1unfitted, subjre)

    vcu, vc = VarCorr(fm1unfitted), VarCorr(fm1)

    @test all(values(vcu.σρ.subj.σ) .≈ values(vc.σρ.subj.σ))
    @test all(values(vcu.σρ.subj.ρ) .≈ values(vc.σρ.subj.ρ))

    @test_throws ArgumentError update!(fm1unfitted; subjitem=subjre, θ=fm1.θ)
    @test_throws ArgumentError update!(fm1unfitted)
    @test_throws ArgumentError update!(fm1unfitted; item=subjre)

    fmkb = fit(MixedModel,
               @formula(rt_trunc ~ 1 + spkr + prec + load + (1 | subj) + (1 + prec | item)),
               MixedModels.dataset(:kb07); progress=false)

    fmkb2 = LinearMixedModel(@formula(rt_trunc ~ 1 + spkr + prec + load + (1 | subj) +
                                                 (1 + prec | item)),
                             MixedModels.dataset(:kb07))
    update!(fmkb2; subj=fmkb.reterms[2].λ, item=fmkb.reterms[1].λ)

    vcu, vc = VarCorr(fmkb2), VarCorr(fmkb2)

    @test all(values(vcu.σρ.subj.σ) .≈ values(vc.σρ.subj.σ))
    @test all(values(vcu.σρ.subj.ρ) .≈ values(vc.σρ.subj.ρ))
    @test all(values(vcu.σρ.item.σ) .≈ values(vc.σρ.item.σ))
    @test all(values(vcu.σρ.item.ρ) .≈ values(vc.σρ.item.ρ))
end
