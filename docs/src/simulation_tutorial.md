Power Analysis and Simulation Tutorial
======================================

This tutorial demonstrates how to conduct power analyses and data simulation using Julia and the MixedModelsSim package.

Power analysis is an important tool for planning an experimental design. Here we show how to
1. Take existing data and calculate power by simulate new data with bootstrapping.
2. Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set. 
3. Create a (simple) fully crossed dataset from scratch and analyze power.
4. Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.




### Load the packages we'll be using in Julia
First, here are the packages needed in this example.

```@example Main
using MixedModels        # run mixed models
using MixedModelsSim     # simulation functions for mixed models
using RCall              # call R functions from inside Julia
using DataFrames, Tables # work with data tables
using StableRNGs         # random number generator
using CSV                # write CSV files
using Markdown
using Statistics         # basic math funcions
using DataFramesMeta     # dplyr-like operations
using Gadfly             # plotting package

using LinearAlgebra      # not used yet, for specifing θ

```

### Define number of iterations
Here we define how many model simulations we want to do. A large number will give more reliable results, but will take longer to compute. It is useful to set it to a low number for testing, and increase it for your final analysis.

```@example Main
nsims = 1000
```

# 1. Take existing data and calculate power by simulate new data with bootstrapping.
## **1.1. Build a Linear Mixed Model from existing data.**

For the first example we are going to simulate bootstrapped data from an existing data set:

*Experiment 2 from Kronmüller, E., & Barr, D. J. (2007). Perspective-free pragmatics: Broken precedents and the recovery-from-preemption hypothesis. Journal of Memory and Language, 56(3), 436-455.*

The data we will be using through out this tutorial is a study about how in a conversation the change of a speaker or the change of precedents (which are patterns of word usage to discribe an object, e.g. one can refer to the same object "white shoes", "runners", "sneakers") affects the understanding.

Objects are presented on a screen while participants listen to instructions to move the objects around. Participants eye movements are tracked.
The dependent variable is response time, defined as the latency between the onset of the test description and the moment at which the target was selected.
The independent variables are speaker (old vs. new), precedents (maintain vs. break)  and cognitive load (a secondary memory task).

We have to load the data and define some characteristics like the contrasts and the underlying model.

Load existing data:
```@example Main
kb07 = MixedModels.dataset("kb07");
```

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());
```

The chosen LMM for this dataset is defined by the following model formula:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Fit the model
```@example Main
kb07_m = fit(MixedModel, kb07_f, kb07, contrasts=contrasts);
print(kb07_m)
```

## **1.2 Simulate from existing data with same model parameters**

We will first look at the power of the dataset with the same parameters as in the original data set. This means that each dataset will have the exact number of observations as the original data. Here, we use the model `kb07_m` we fitted above to our dataset `kb07`.

You can use the `parameparametricbootstrap()` function to run `nsims` iterations of data sampled using the parameters from `kb07_m`.
Set up a random seed to make the simulation reproducible. You can use your favourite number.

To use multithreading, you need to set the number of cores you want to use.
E.g. in Visual Studio Code, open the settings (gear icon in the lower left corner or cmd-,) and search for "thread".
Set `julia.NumThreads` to the number of cores you want to use (at least 1 less than your total number).


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Run nsims iterations:
```@example Main
kb07_sim = parametricbootstrap(rng, nsims, kb07_m, use_threads = false);
```
**Try**: Run the code above with or without `use_threads = true`.

The output DataFrame `kb07_sim` contains the results of the bootstrapping procedure. 

```@example Main
df = DataFrame(kb07_sim.allpars);
first(df, 9)
nrow(df)
```

The dataframe df has 9000 rows: 9 parameters, in 1000 iterations.

Plot some bootstrapped parameter:
```@example Main
σres = @where(df, :type .== "σ", :group .== "residual").value
plot(x = σres, Geom.density, Guide.xlabel("Parametric bootstrap estimates of σ"), Guide.ylabel("Density"))

βInt = @where(df, :type .== "β", :names .== "(Intercept)").value
plot(x = βInt, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β (Intercept)"), Guide.ylabel("Density"))

βSpeaker = @where(df, :type .== "β", :names .== "spkr: old").value
plot(x = βSpeaker, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Speaker"), Guide.ylabel("Density"))

βPrecedents = @where(df, :type .== "β", :names .== "prec: maintain").value
plot(x = βPrecedents, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Precedents"), Guide.ylabel("Density"))

βLoad = @where(df, :type .== "β", :names .== "load: yes").value
plot(x = βLoad, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Load"), Guide.ylabel("Density"))
```

Convert p-values to dataframe and save it as CSV
```@example Main
kb07_sim_df = DataFrame(kb07_sim.coefpvalues);
CSV.write("kb07_sim.csv", kb07_sim_df);
```

Have a look at your simluated data:
```@example Main
print(first(kb07_sim_df, 8))
```

Now that we have a bootstrapped data, we can start our power calculation.

### Power calculation

The function `power_table()` from `MixedModelsSim` takes the output of `parametricbootstrap()` and calculates the proportion of simulations where the p-value is less than alpha for each coefficient.
You can set the `alpha` argument to change the default value of 0.05 (justify your alpha ;).

```@example Main
ptbl = power_table(kb07_sim, 0.05)
```

A powervalue of 1 maens that in every iteration the spefific parameter we are looking at was below our alpha.
A powervalue of 0.5 means that in half of our iterations the spefific parameter we are looking at was below our alpha.

You can also do it manually:
```@example Main
kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:]

mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05)
````

For nicely displaying, you can use `pretty_table`:
```@example Main
pretty_table(ptbl)
```

# 2. Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set.

Let's say we want to check our power to detect effects of spkr, prec, and load
that are only half the size as in our pilot data. We can set a new vector of beta values
with the `β` argument to `parametricbootstrap()`.


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Specify β:
```@example Main
new_beta = kb07_m.β
new_beta[2:4] = kb07_m.β[2:4]/2
```

Run nsims iterations:
```@example Main
kb07_sim_half = parametricbootstrap(rng, nsims, kb07_m, β = new_beta, use_threads = false);
```

### Power calculation

```@example Main
power_table(kb07_sim_half)
```

# 3. Create a (simple) fully crossed dataset from scratch and analyze power.

In some situations, instead of using an existing dataset it may be useful to simulate the data from scratch. This could be the case when the original data is not available but the effect sizes are known. That means that we have to:

a) specify the effect sizes manually
b) manually create an experimental design, according to which data can be simulated

If we simulate data from scratch, aside from subject and item number, we can manipulate the arguments `β`, `σ` and `θ`.
Lets have a closer look at them, define their meaning and we will see where the corresponding values in the model output are.

### **Beta**
`β` are our effect sizes. If we look again on our LMM summary from the `kb07`-dataset `kb07_m`
we see our four `β` under fixed-effects parameters in the `Coef.`-column.

```@example Main
kb07_m
kb07_m.β
```

### **Sigma**
`σ` is the residual-standard deviation listed under the variance components.
```@example Main
kb07_m
kb07_m.σ
```

### **Theta**
The meaning of `θ` is a bit less intuitive. In a less complex model (one that only has intercepts for the random effects) or if we supress the correlations in the formula with `zerocorr()` then `θ` describes the relationship between the random effects standard deviation and the standard deviation of the residual term.
In our `kb07_m` example:
The residual standard deviation is `680.032`.
The standard deviation of our first variance component *`item - (Intercept)`* is `364.713`.
Thus our first `θ` is the relationship: variance component devided by residual standard deviation
364.713 /  680.032 =  `0.53631`

```@example Main
kb07_m.θ
```

We also can calculate the `θ` for variance component *`subj - (Intercept)`*.
The residual standard deviation is `680.032`.
The standard deviation of our variance component *`subj - (Intercept)`* is `298.026`.
Thus, the related θ is the relationship: variance component devided by residual standard deviation
298.026 /  680.032 =  `0.438252`

```@example Main
kb07_m.θ
```

We can not calculate the `θ` for variance component *`item - prec: maintain`* yet, because it includes the correlation of
*`item - prec: maintain`* and *`item - (Intercept)`*.
The `θ` vector is the flattened version of the variance-covariance matrix - a lowertrinangular matrix.
The on-diagonal elements are just the standard deviations (the `σ`'s), If all off-diagonal elements are zero, we can use our
calculation above. The off-diagonal elements are covariances and correspond to the correlations (the `ρ`'s).
If they are unequal to zero, as it is in our `kb07`-dataset, we cannot recreate the variance-covariance matrix having the model output.
We just take it from the model we already fitted

See the two inner values:
```@example Main
kb07_m.θ
```

Play with theta
```@example Main
kb07_m2 = kb07_m
update!(kb07_m2, re_item, re_subj)
```

make a lower triangular matrix
```@example Main
re_item = create_re(0.5363168233715857,0.37133693708531357)
re_item[2]=-0.70
re_item

re_subj = create_re(0.4382528181348316)
```

make the compact form out of it = is equal to θ
```@example Main
newθ= vcat( flatlowertri(re_item), flatlowertri(re_subj) )
```


## *A simple example*

Having this knowledge about the parameters we can now **simulate data from scratch**

The `simdat_crossed()` function from `MixedModelsSim` lets you set up a data frame with a specified experimental design.
For now, it only makes fully balanced crossed designs!, but you can generate an unbalanced design by simulating data for the largest cell and deleting extra rows.

TODO: NEED HELP: is that still right? I think it is possible now to make partially crossed designs? 
===========

Firstly we will set an easy design where `subj_n` subjects per `age` group (O or Y) respond to `item_n` items in each of two `condition`s (A or B).

Your factors need to be specified separately for between-subject, between-item, and within-subject/item factors using `Dict` with the name of each factor as the keys and vectors with the names of the levels as values.

First we have to define factors in a dict.

We start with the between subject factors:
```@example Main
subj_btwn = Dict("age" => ["O", "Y"])
```

There are no between-item factors in this design so you can omit it or set it to nothing.
Note that if you have factors which are between subject and between item, you need to put them in both dicts.
```@example Main
item_btwn = nothing
```

Next, we put within-subject/item factors in a dict:
```@example Main
both_win = Dict("condition" => ["A", "B"])
```

Define subject and item number:
 ```@example Main
subj_n = 10
item_n = 30
```

Simulate data:
```@example Main
dat = simdat_crossed(subj_n,
                     item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);
```

Have a look:
```@example Main
first(DataFrame(dat),8)
```
The values we see in the column `dv` is just random noise.

Set contrasts:
```@example Main
contrasts = Dict(:age => HelmertCoding(),
                 :condition => HelmertCoding());
```

Define formula:
```@example Main
f1 = @formula dv ~ 1 + age * condition + (1|item) + (1|subj);
```

Fit the model:
```@example Main
m1 = fit(MixedModel, f1, dat, contrasts=contrasts)
print(m1)
```

Because the `dv` is just random noise from N(0,1), there will be basically no subject or item random variance, residual variance will be near 1.0, and the estimates for all effects should be small.
Don't worry, we'll specify fixed and random effects directly in `parametricbootstrap()`.


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Specify `β`, `σ`, and `θ`, we just made up this parameter:
```@example Main
new_beta = [0., 0.25, 0.25, 0.]
new_sigma = 2.0
new_theta = [1.0, 1.0]
```

Run nsims iterations:
```@example Main
sim1 = parametricbootstrap(rng, nsims, m1,
                        β = new_beta,
                        σ = new_sigma,
                        θ = new_theta,
                        use_threads = false);
```

### Power calculation

```@example Main
ptbl= power_table(sim1)
```

For nicely displaying it, you can use pretty_table:
```@example Main
pretty_table(ptbl)
```

# 4. Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.

# *4.1 Recreate the `kb07`-dataset from scratch*
For full control over all parameters in our `kb07` data set we will recreate the design using the method shown above.


Define subject and item number:
```@example Main
subj_n = 56
item_n = 32
```

Define factors in a dict:
```@example Main
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);
```

Simulate data:
```@example Main
fake_kb07 = simdat_crossed(subj_n, item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);
```

Make a dataframe:
```@example Main
fake_kb07_df = DataFrame(fake_kb07)
```

Have a look:
```@example Main
first(fake_kb07_df,8)
```

The function `simdat_crossed` generates a fully crossed design.
Unfortunately, our original design is not fully crossed. Every subject saw an image only once, thus in one of eight possible conditions. To simulate that we only keep one of every eight lines.

TODO: NEED HELP, is that correct? or is it possible to do it in simdat_crossed?
========

We sort the dataframe to enable easier selection
```
fake_kb07_df = sort(fake_kb07_df, [:subj, :item, :load, :prec, :spkr])
```

In order to select only the relevant rows of the data set we define an index which represets a random choice of one of every eight rows. First we generate a vector `idx` which represents which row to keep in each set of 8.
```@example Main
len = convert(Int64,(length(fake_kb07)/8))
idx = rand(rng, 1:8 , len)
```

Then we create an array `A`, of the same length that is populated multiples of the number 8. Added together `A` and `idx` give the indexes of one row from each set of 8s.
```@example Main
A = repeat([8], inner=len-1)
A = append!( [0], A )
A = cumsum(A)
idx = idx+A
```

Reduce the fully crossed design to the original experimental design:
```@example Main
fake_kb07_df= fake_kb07_df[idx, :]
rename!(fake_kb07_df, :dv => :rt_trunc)
```

Write a CSV:
```@example Main
CSV.write("fake_kb07_df.csv", fake_kb07_df);
```


Now we can use the simulated data in the same way as above.

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());
```

Define formula, same as above:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Fit the model:
```@example Main
fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df, contrasts=contrasts);
print(fake_kb07_m)
```

Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Then, again, we specify `β`, `σ`, and `θ`.
Here we use the values that we found in the model of the existing dataset:

```@example Main
new_beta = [2181.85, 67.879, -333.791, 78.5904]
new_beta = kb07_m.β

new_sigma = 680.032
new_sigma = kb07_m.σ

new_theta = [0.5363168233715857,
           -0.25981993107379,
            0.2653016002105174,
            0.4382528181348316]
new_theta = kb07_m.θ
```


Run nsims iterations:
```@example Main
fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m,
                        β = new_beta,
                        σ = new_sigma,
                        θ = new_theta,
                        use_threads = false);
```

### Power calculation

```@example Main
power_table(fake_kb07_sim)
```

Compare to the powertable from the existing data:
```@example Main
power_table(kb07_sim)
```

We have successfully recreated the power simulation of an existing dataset from scratch. This has the advantage, that we now can iterate over different numbers of subjects and items.

# *4.2 Loop over subject and item sizes*

When designing a study, you may be interested in trying various numbers of subjects and items to see how that affects the power of your study.
To do this you can use a loop to run simulations over a range of values for any parameter.

### We first define every fixed things outside the loop (same as above):

Define factors in a dict:
```@example Main
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);
```

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());
```

Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Define formula:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Specify `β`, `σ`, and `θ`:
```@example Main
new_beta = [2181.85, 67.879, -333.791, 78.5904]
new_beta = kb07_m.β

new_sigma = 680.032
new_sigma = kb07_m.σ

new_theta = [0.5363168233715857,
           -0.25981993107379,
            0.2653016002105174,
            0.4382528181348316]
new_theta = kb07_m.θ
```


### Then we define the variables that out loop will iterate over

Define subject and item numbers as arrays:
```@example Main
sub_ns = [20, 30, 40];
item_ns = [10, 20, 30];
```

Make an emty dataframe:
```@example Main
d = DataFrame();
```

### Run the loop:
```@example Main
for subj_n in sub_ns
    for item_n in item_ns

    # Make fully crossed data:
    fake_kb07 = simdat_crossed(subj_n, item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);
    fake_kb07_df = DataFrame(fake_kb07)

    # Reduce the fully crossed design to the original experimental design:
    fake_kb07_df = sort(fake_kb07_df, [:subj, :item, :load, :prec, :spkr])
    len = convert(Int64,(length(fake_kb07)/8))
    idx = rand(rng, 1:8 , len)
    A = repeat([8], inner=len-1)
    A = append!( [0], A )
    A = cumsum(A)
    idx = idx+A
    fake_kb07_df= fake_kb07_df[idx, :]
    rename!(fake_kb07_df, :dv => :rt_trunc)

    # Fit the model:
    fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df, contrasts=contrasts);

    # Run nsims iterations:
    fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m,
                        β = new_beta,
                        σ = new_sigma,
                        θ = new_theta,
                        use_threads = false);

    # Power calculation
    ptbl = power_table(fake_kb07_sim)
    ptdf = DataFrame(ptbl)
    ptdf[!, :item_n] .= item_n
    ptdf[!, :subj_n] .= subj_n
    append!(d, ptdf)

    end
end
```

Our dataframe `d` now contains the power information for each combination of subjects and items.

Save the powertable as CSV
```@example Main
CSV.write("power.csv", d)
```

TODO: NEED Help: it would be nice to make a plot of this!
=========

Here we show how to make a quick plot of the output
```
using StatsPlots

@df d scatter(
    :subj_n,
    :power,
    group = :item_n,
    title = "Power of different designs",
    xlabel = "# of Subjects",
    ylabel = "power",
    m = (0.5, [:cross :hex :star7], 12),
    bg = RGB(0.2, 0.2, 0.2)
)

```

# Credit
This tutorial was conceived for ZiF research and tutorial workshop by Lisa DeBruine (Feb. 2020) presented again by Phillip Alday during the SMLP Summer School (Sep. 2020).

Updated and extended by Lisa Schwetlick & Daniel Backhaus, with the kind help of Phillip Alday,  after changes to the package.