using MixedModelsSim
using Tables
using Test


@testset "powersimulation" begin
	kb07 = MixedModels.dataset(:kb07);
	form = @formula(rt_raw ~ 1 + spkr + prec + load + (1+spkr+prec+load|subj) + (1+spkr+prec+load|item));
	cont = Dict(:spkr => HelmertCoding(),
	            :prec => HelmertCoding(),
	            :load => HelmertCoding())
	fm1 = fit(MixedModel, form, kb07, contrasts=cont, REML=false);
	zpmt = simulate_waldtests(MersenneTwister(42),10,fm1,use_threads=true);
	pvals = getindex.(columntable(zpmt).p,1)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power == 1.
	pvals = getindex.(columntable(zpmt).p,2)
	power = sum(pvals .< 0.05) / length(pvals)
	@test power == 0.6
end
