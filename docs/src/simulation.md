Simulating an Experiment from Scratch
======================================

Here is a worked example of simulating a partially crossed design from scratch.

First, some setup:

```@example Main
using DataFrames
using LinearAlgebra
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
For nicely displaying it, we can also convert it to a dataframe and change the factors to use pooled arrays to save a bit of memory.

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

```@example Main
vc.σρ
```

In matrix form:

```@example Main
m0.λ
```

```@example Main
m0.λ[1]
```

```@example Main
m0.λ[2]
```

Note that the entries in `m0.λ` are stored in the same order as the random effects are listed in `VarCorr(m0)`.  The off-diagonal elements are covariances and correspond to the correlations (the ρ's).
The on-diagonal elements are just the standard deviations (the σ's).
For this example, we'll just assume all the correlations and thus the covariances are 0 in order to make things simple.
Then we only have to worry about the diagonal elements.

Let's assume that the variability
- between items
    - in the intercept is 1.3 times the residual variability
    - in age is 0.35 times the residual variability
    - in context is 0.75 times the residual variability
- between subjects
    - in the intercept is 1.5 times the residual variability
    - in frequency is 0.5 times the residual variability
    - in context is 0.75 times the residual variability

Then we'll have the following λ matrices:

```@example Main
λitem = LowerTriangular(diagm([1.3, 0.35, 0.75]))
```

```@example Main
λsubj = LowerTriangular(diagm([1.5, 0.5, 0.75]))
```

For the actual optimization process, MixedModels.jl actually uses a flattened version of these stored in the θ vector. (More precisely, these λ matrices are derived from the θ vector.)

```@example Main
isapprox(m0.θ,  [flatlowertri(m0.λ[1]); flatlowertri(m0.λ[2])]; atol=0.1)
```

With that in mind, we can assemble our θ vector for simulation:

```@example Main
# make sure these are in the same order as in the model summary!
θ = [flatlowertri(λitem); flatlowertri(λsubj)]
```

## Assemble the Fixed Effects

The last two components we need are the residual variance and the effect sizes for the fixed effects.

```@example Main
σ = 5
β = [1.0, -1.0, 2.0, -1.5, 0.3, -1.3, 1.4, 0]
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
## See your power and profit!

Finally, we can turn this into a power table:

```@example Main
ptbl = power_table(sim)
```

We can of course convert the row table into a DataFrame:
```@example Main
DataFrame(ptbl)
```
