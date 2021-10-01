Simulating an Experiment from Scratch
======================================

Here is a worked example of simulating a partially crossed design from scratch.

First, some setup:

```@example Main
using DataFrames
using MixedModels, MixedModelsSim
using Random
```

## Assemble the Design

We're going to do a 2 x 2 x 2 design.
For concreteness, let's think of this as a linguistic design:

- **age** `old` vs. `young`, _between subjects_
- **frequency** `high` vs. `low`, _between items_
- **context** `matched` vs. `unmatched`, _within both_.

Further, let's assume we want 40 subjects and 80 items.
We can specify this design as follows:

```@example Main
n_subj = 40
n_item = 80
subj_btwn = Dict(:age => ["old", "young"])
item_btwn = Dict(:frequency => ["high", "low"])
both_win = Dict(:context => ["matched", "unmatched"])
```

and then generate it:

```@example Main
rng = MersenneTwister(42)  # specify our random number generator for reproducibility
design = simdat_crossed(rng, n_subj, n_item;
                        subj_btwn = subj_btwn,
                        item_btwn = item_btwn,
                        both_win = both_win)
```

Note that `simdat_crossed` returns a row table, which `MixedModels.jl` can process directly.
For nicely displaying it, we can again use `pretty_table`:
```@example Main
pretty_table(first(design, 5))
```

We can also convert it to a DataFrame and change the factors to use pooled arrays to save a bit of memory.

```@example Main
df = pooled!(DataFrame(design))
first(df, 5)
```

Note that `simdat_crossed` generated a column `dv` for our dependent variable that has been pre-populated with noise from a standard normal distribution ($N(0,1)$).
Typically, we will want to scale that, but we can do that in the simulation step.
Also, this dependent variable is *pure noise*: we haven't yet added in effects.
Adding in effects also comes in the simulation step.

But before we get to simulating, let's fit the model to the noise, just to see how things look. We're going to use effects coding for our contrasts.

```@example Main
contrasts = Dict(:age => EffectsCoding(base="young"),
                 :frequency => EffectsCoding(base="high"),
                 :context => EffectsCoding(base="matched"))
form = @formula(dv ~ 1 + age * frequency * context +
                    (1 + frequency + context | subj) +
                    (1 + age + context | item))
m0 = fit(MixedModel, form, design; contrasts=contrasts)
```

## Assemble the Random Effects

The hard part in simulating right now is specifying the random effects.
We're working on making this bit easier, but you need to specify the variance-covariance matrix of the random effects. You can see what this
looks like:

```@example Main
vc = VarCorr(m0)
```

For each grouping variable (subjects and items), there are two major components: the standard deviations ahd the correlations.

For this example, we'll just assume all the correlations and thus the covariances are 0 in order to make things simple.
Then we only have to worry about the standard deviations.

Let's assume that the variability
- between items
  - in the intercept is 1.3 times the residual variability
  - in age is 0.35 times the residual variability
  - in context is 0.75 times the residual variability
- between subjects
  - in the intercept is 1.5 times the residual variability
  - in frequency is 0.5 times the residual variability
  - in context is 0.75 times the residual variability

Note these are always specified relative to the residual standard deviation.
In other words, we think about how big the between-subject and between-item differences are relative to the between-observation differences.

We can now create the associated covariance matrices.[^cholesky]

[^cholesky]: Technically, we're creating the lower Cholesky factor of these matrices, which is a bit like the matrix square root. In other words, we're creating the matrix form of standard deviations instead of the matrix form of the variances.

```@example Main
re_item = create_re(1.3, 0.35, 0.75)
```

```@example Main
re_subj = create_re(1.5, 0.5, 0.75)
```

We can check that we got these right by installing these parameter values into the model.
```@example Main
update!(m0; subj=re_subj, item=re_item)
VarCorr(m0)
```

Looks good. The values don't exactly match the values in our parameter vector because the
residual standard deviation isn't exactly 1.0.

For the actual simulation, we'll need the compact form of these covariance matrices that MixedModels.jl uses internally.
This compact form is the parameter vector θ and we can get it back out of the model where we just installed it:
```@example Main
show(m0.θ)
```

Alternatively, we can also just generate θ directly from the random-effects matrices:
```@example Main
θ = createθ(m0; subj=re_subj, item=re_item)
show(θ)
```

We could also create it directly from the covariance matrices we created, but in this case we need to make sure they're in the same order as in the `VarCorr` output:
```@example Main
θhand = vcat( flatlowertri(re_item), flatlowertri(re_subj) )
show(θhand)
```

!!! warning
    In storing the parameter vector θ, MixedModels.jl uses an ordering that yields maximum sparseness, which enables better computational efficiency.
    The ordering is thus dependent on the entirety of the design -- both the choice of the random effects and the relative number of subjects and items.
    For this reason, we strongly recommend using the helper methods that allow specifying the grouping variable by name.

## Assemble the Fixed Effects

The last two components we need are the residual variance and the effect sizes for the fixed effects.

```@example Main
σ = 5;
β = [1.0, -1.0, 2.0, -1.5, 0.3, -1.3, 1.4, 0];
nothing # hide
```

The entries in the β correspond to the coefficients in the model given by

```@example Main
coefnames(m0)
```

## Simulate

Now we're ready to actually simulate our data.
We can use `parametricbootstrap` to do this: the parametric bootstrap actually works by simulating new data from an existing model and then looking at how the estimates fit to that new data look.
In MixedModels.jl, you can specify different parameter values, such as the ones
 we made up for our fake data.

```@example Main
# typically we would use a lot more simulations
# but we want to be quick in this example
sim = parametricbootstrap(MersenneTwister(12321), 20, m0; β=β, σ=σ, θ=θ)
```

As mentioned above, the ordering within θ is dependent on the entire design, so if you embed the simulation code in a loop iterating over different numbers of subjects and items, it's probably better to write it as:
```julia
sim = parametricbootstrap(MersenneTwister(12321), 20, m0;
                          β=β, σ=σ, θ=createθ(m0; subj=re_subj, item=re_item))
```

## See your power and profit!

Finally, we can turn this into a power table:

```@example Main
ptbl = power_table(sim)
```

For nicely displaying it, we can again use `pretty_table`:
```@example Main
pretty_table(ptbl)
```

We can of course convert the row table into a DataFrame:
```@example Main
DataFrame(ptbl)
```
