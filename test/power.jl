using MixedModelsSim
using MixedModels
using Random
using Tables
using Test


@testset "powersimulation" begin
	kb07 = MixedModels.dataset(:kb07);
	form = @formula(rt_raw ~ 1 + spkr + prec + load + (1+spkr+prec+load|subj) + (1+spkr+prec+load|item));
	cont = Dict(nm => HelmertCoding() for nm in (:spkr, :prec, :load));
	fm1 = fit(MixedModel, form, kb07, contrasts=cont, REML=false);
	zpmt = simulate_waldtests(MersenneTwister(42),10,fm1,use_threads=false);
	pvals = getindex.(columntable(zpmt).p, 1)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power == 1.
	@test getindex.(columntable(zpmt).p,1) == getindex.(columntable(zpmt).p,Symbol("(Intercept)"))
	pvals = getindex.(columntable(zpmt).p,2)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power == 0.9
end
