using MixedModelsSim
using MixedModels
using StableRNGs
using Statistics
using Tables
using Test


subj_n = 30
item_n = 30
subj_btwn = Dict(:age => ["O", "Y"],
				:pet => ["cat", "dog"])
item_btwn = Dict(:cond => ["A", "B"])
both_win = Dict(:time => ["morning", "evening"])
dat = simdat_crossed(subj_n, item_n, 
					 subj_btwn = subj_btwn, 
					 item_btwn = item_btwn, 
					 both_win = both_win)

@testset "powersimulation" begin
	form = @formula(dv ~ 1 + age + pet + cond + time  + (1|subj) + (1|item));
	cont = Dict(nm => HelmertCoding() for nm in (:age, :pet, :cond, :time));
	fm1 = fit(MixedModel, form, dat, contrasts=cont, REML=false);
	zpmt = simulate_waldtests(StableRNG(42),10,fm1);
	pvals = getindex.(columntable(zpmt).p, 1)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power ≈ 0.1
	@test getindex.(columntable(zpmt).p,1) == getindex.(columntable(zpmt).p,Symbol("(Intercept)"))
	pvals = getindex.(columntable(zpmt).p,2)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power ≈ 0.0
end
