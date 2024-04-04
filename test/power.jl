using MixedModelsSim, MixedModels, StableRNGs
using DataFrames, Tables
using Statistics
using Test

subj_n = 40
item_n = 20
subj_btwn = Dict(:age => ["O", "Y"],
                 :pet => ["cat", "dog"])
item_btwn = Dict(:cond => ["A", "B"])
both_win = Dict(:time => ["morning", "evening"])
dat = simdat_crossed(StableRNG(42), subj_n, item_n;
                     subj_btwn=subj_btwn,
                     item_btwn=item_btwn,
                     both_win=both_win)
form = @formula(dv ~ 1 + age + pet + cond + time + (1 | subj) + (1 | item));
cont = Dict(nm => HelmertCoding() for nm in (:age, :pet, :cond, :time));
fm1 = fit(MixedModel, form, dat; contrasts=cont, REML=false, progress=false);
sim = parametricbootstrap(StableRNG(42), 10, fm1; use_threads=false, hide_progress=true);

@testset "power_table" begin
    pt = power_table(sim)

    @test Tables.isrowtable(pt)

    ptdf = sort!(DataFrame(pt), :coefname)

    @test nrow(ptdf) == 5
    @test ncol(ptdf) == 2
    @test names(ptdf) == ["coefname", "power"]
    @test ptdf.power == [0.4, 0.2, 0.0, 0.3, 0.1]
end
