using MixedModelsSim, MixedModels, StableRNGs
using DataFrames, Tables
using Statistics
using Test

subj_n = 10
item_n = 10
subj_btwn = Dict(:age => ["O", "Y"],
				:pet => ["cat", "dog"])
item_btwn = Dict(:cond => ["A", "B"])
both_win = Dict(:time => ["morning", "evening"])
dat = simdat_crossed(StableRNG(42), subj_n, item_n,
					 subj_btwn = subj_btwn,
					 item_btwn = item_btwn,
                     both_win = both_win)
form = @formula(dv ~ 1 + age + pet + cond + time  + (1|subj) + (1|item));
cont = Dict(nm => HelmertCoding() for nm in (:age, :pet, :cond, :time));
fm1 = fit(MixedModel, form, dat, contrasts=cont, REML=false);
zpmt = simulate_waldtests(StableRNG(42),10,fm1,use_threads=false);
bsamp = parametricbootstrap(StableRNG(42),10,fm1,use_threads=false);


@testset "sim_to_df" begin
    df = sim_to_df(zpmt)

    @test nrow(df) == 50
    @test ncol(df) == 6
    @test names(df) == ["iteration", "coefname", "beta", "se", "z", "p"]
    # should work after coefname order is fixed
    # @test string.(df.coefname[1:4]) == collect(string.(zpmt[1][1] |> keys))
end

@testset "power_table" begin
    pt = power_table(zpmt)

    @test nrow(pt) == 5
    @test ncol(pt) == 2
    @test names(pt) == ["coefname", "power"]
    @test string.(pt.coefname) == collect(string.(zpmt[1][1] |> keys))
    @test pt.power == [0.3, 0.2, 0.1, 0.0, 0.3]
end
