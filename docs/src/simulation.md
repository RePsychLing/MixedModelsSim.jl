Simulating an Experiment
=========================


```julia
using MixedModels, MixedModelsSim
using DataFrames, Gadfly, Random

kb07 = MixedModels.dataset(:kb07);
form = @formula(rt_raw ~ 1 + spkr + prec + load + (1+spkr+prec+load|subj) + (1+spkr+prec+load|item));
cont = Dict(:spkr => HelmertCoding(),
            :prec => HelmertCoding(),
            :load => HelmertCoding())
fm1 = fit(MixedModel, form, kb07, contrasts=cont, REML=false);
sim = parametricbootstrap(MersenneTwister(42),10,fm1,use_threads=true);


shortestcovint(samples[!,Symbol("prec: maintain")])

plot(samples, layer(x=Symbol("load: yes")),
    Geom.density,
    Guide.xlabel("Bootstrap estimates of 'load: yes'"))

plot(stack(samples),
   layer(x=:value, color=:variable, xgroup=:variable),
   Geom.subplot_grid(Geom.density, free_x_axis=true),
   Guide.xlabel("Bootstrap estimates of β"))
```
