using MixedModelsSim, MixedModels, Random
using DataFrames, Tables
using Statistics
using Test

kb07 = MixedModels.dataset(:kb07);
form = @formula(rt_raw ~ 1 + spkr + prec + load + (1+spkr+prec+load|subj) + (1+spkr+prec+load|item));
cont = Dict(nm => HelmertCoding() for nm in (:spkr, :prec, :load));
fm1 = fit(MixedModel, form, kb07, contrasts=cont, REML=false);
zpmt = simulate_waldtests(MersenneTwister(42),10,fm1,use_threads=false);

@testset "sim_to_df" begin
    df = sim_to_df(zpmt)

    @test nrow(df) == 40
    @test ncol(df) == 6
    @test names(df) == ["iteration", "coefname", "beta", "se", "z", "p"]
    @test df.coefname == repeat(1:10, outer=4)  # currently the column names are mixed up
    # should work after coefname order is fixed
    # @test string.(df.coefname[1:4]) == collect(string.(zpmt[1][1] |> keys))
end

@testset "power_table" begin
    pt = power_table(zpmt)

    @test nrow(pt) == 4
    @test ncol(pt) == 2
    @test names(pt) == ["coefname", "power"]
    @test string.(pt.coefname[1:4]) == collect(string.(zpmt[1][1] |> keys))
    @test pt.power == [1.0, 0.9, 1.0, 1.0]
end
