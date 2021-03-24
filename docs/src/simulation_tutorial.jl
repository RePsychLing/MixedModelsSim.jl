### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ b8016ebc-8cab-11eb-0bd3-5bea764c4251
begin
using MixedModels        # run mixed models
using MixedModelsSim     # simulation functions for mixed models
using DataFrames, Tables # work with data tables
using StableRNGs         # random number generator
using CSV                # write CSV files
#using Markdown
using Statistics         # basic math funcions
using DataFramesMeta     # dplyr-like operations
using Gadfly             # plotting package
end

# ╔═╡ 3edca7b8-8cab-11eb-3982-0bfba7ea51b3
md"""

Power Analysis and Simulation Tutorial
======================================

This tutorial demonstrates how to conduct power analyses and data simulation using Julia and the MixedModelsSim package.

Power analysis is an important tool for planning an experimental design. Here we show how to
1. Take existing data and calculate power by simulate new data with bootstrapping.
2. Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set.
3. Create a (simple) balanced fully crossed dataset from scratch and analyze power.
4. Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.




### Load the packages we'll be using in Julia
First, here are the packages needed in this example.
"""

# ╔═╡ 7cf5bc5a-8cac-11eb-1b2d-25987ca50e72
md"""
### Define number of iterations
Here we define how many model simulations we want to do. A large number will give more reliable results, but will take longer to compute. It is useful to set it to a low number for testing, and increase it for your final analysis.
"""

# ╔═╡ 77079142-8cac-11eb-2894-3b1af57b13cd
nsims = 500

# ╔═╡ a6254e9c-8cac-11eb-2ecc-8f141b7fb022
md"""
# 1. Take existing data and calculate power by simulate new data with bootstrapping.
## **1.1. Build a Linear Mixed Model from existing data.**

For the first example we are going to simulate bootstrapped data from an existing data set:

*Experiment 2 from Kronmüller, E., & Barr, D. J. (2007). Perspective-free pragmatics: Broken precedents and the recovery-from-preemption hypothesis. Journal of Memory and Language, 56(3), 436-455.*

The data we will be using through out this tutorial is a study about how in a conversation the change of a speaker or the change of precedents (which are patterns of word usage to discribe an object, e.g. one can refer to the same object "white shoes", "runners", "sneakers") affects the understanding.

Objects are presented on a screen while participants listen to instructions to move the objects around. Participants eye movements are tracked.
The dependent variable is response time, defined as the latency between the onset of the test description and the moment at which the target was selected.
The independent variables are speaker (old vs. new), precedents (maintain vs. break) and cognitive load (a secondary memory task).

We have to load the data and define some characteristics like the contrasts and the underlying model.

Load existing data:
"""

# ╔═╡ c2a36428-8cac-11eb-1eba-85a5bba4f980
kb07 = MixedModels.dataset(:kb07);

# ╔═╡ dd95cf64-8cac-11eb-14b4-a91cf59e9094
md"""
Set contrasts:
"""

# ╔═╡ ef62641e-8cac-11eb-1b85-7be4b0eea3ea
md"""
The chosen LMM for this dataset is defined by the following model formula:
"""

# ╔═╡ 0d827ad8-8cad-11eb-2c4c-9348874e0343
md"""
Fit the model
"""

# ╔═╡ 2702a6d6-8cad-11eb-3922-5fa14ce5c7b8
md"""
## **1.2 Simulate from existing data with same model parameters**

We will first look at the power of the dataset with the same parameters as in the original data set. This means that each dataset will have the exact number of observations as the original data. Here, we use the model `kb07_m` we fitted above to our dataset `kb07`.

You can use the `parameparametricbootstrap()` function to run `nsims` iterations of data sampled using the parameters from `kb07_m`.
Set up a random seed to make the simulation reproducible. You can use your favourite number.

To use multithreading, you need to set the number of cores you want to use.
E.g. in Visual Studio Code, open the settings (gear icon in the lower left corner or cmd-,) and search for "thread".
Set `julia.NumThreads` to the number of cores you want to use (at least 1 less than your total number).


Set random seed for reproducibility:
"""

# ╔═╡ 3eceee46-8cad-11eb-2229-fdd1dc488c31
md"""
Run nsims iterations:
"""

# ╔═╡ 7d3b3aae-8cad-11eb-209a-8d5a4b726ec6
md"""
**Try**: Run the code above with or without `use_threads = true`.

The output DataFrame `kb07_sim` contains the results of the bootstrapping procedure.
"""

# ╔═╡ 9b7d0cf4-8cad-11eb-3ecc-d3e4e6287837
md"""
The dataframe df has 4500 rows: 9 parameters, each from 500 iterations.

Plot some bootstrapped parameter:
"""

# ╔═╡ 165c3f30-8cae-11eb-037a-f337833479cf
md"""
Convert p-values of your fixed-effects parameters to dataframe 
"""

# ╔═╡ 4be8e2e8-8cae-11eb-2843-c79f6cbad261
md"""
Have a look at your simulated data:
"""

# ╔═╡ 4b7d4632-8cae-11eb-236e-3da3b6e9d9fb
md"""

Now that we have a bootstrapped data, we can start our power calculation.

### Power calculation

The function `power_table()` from `MixedModelsSim` takes the output of `parametricbootstrap()` and calculates the proportion of simulations where the p-value is less than alpha for each coefficient.
You can set the `alpha` argument to change the default value of 0.05 (justify your alpha).

"""

# ╔═╡ 4ad40e64-8cae-11eb-0304-fbc12a867ef1
md"""
An estimated power of 1 means that in every iteration the specific parameter we are looking at was below our alpha.
An estimated power of 0.5 means that in half of our iterations the specific parameter we are looking at was below our alpha.
An estimated power of 0 means that for none of our iterations the specific parameter we are looking at was below our alpha.

You can also do it manually:
"""

# ╔═╡ 6cd06d82-8cae-11eb-0d5b-bdf66b8ddd86
md"""
For nicely displaying, you can use `pretty_table`:
"""

# ╔═╡ 4683428a-8cae-11eb-2cdb-3516864105ce
md"""

# 2. Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set.

Let's say we want to check our power to detect effects of spkr, prec, and load
that are only half the size as in our pilot data. We can set a new vector of beta values
with the `β` argument to `parametricbootstrap()`.


Set random seed for reproducibility:

"""

# ╔═╡ f60e9542-8cae-11eb-3708-7d07915f3c73
#rng = StableRNG(42);

# ╔═╡ f5f3bd44-8cae-11eb-06b5-b330398ec687
md"""
Specify β:
"""

# ╔═╡ f5c0d794-8cae-11eb-13a1-af079fabc06c
md"""
Run nsims iterations:
"""

# ╔═╡ f590368e-8cae-11eb-00cc-2142ff14576c
md"""
### Power calculation

"""

# ╔═╡ f55dd39c-8cae-11eb-17e2-5b4e85500ded
md"""

# 3. Create a (simple) balanced fully crossed dataset from scratch and analyze power.

In some situations, instead of using an existing dataset it may be useful to simulate the data from scratch. This could be the case when the original data is not available but the effect sizes are known. That means that we have to:

a) specify the effect sizes manually

b) manually create an experimental design, according to which data can be simulated

If we simulate data from scratch, aside from subject and item number, we can manipulate the arguments `β`, `σ` and `θ`.
Lets have a closer look at them, define their meaning and we will see where the corresponding values in the model output are.

### **Beta**
`β` are our effect sizes. If we look again on our LMM summary from the `kb07`-dataset `kb07_m`
we see our four `β` under fixed-effects parameters in the `Coef.`-column.

"""

# ╔═╡ f52d4376-8cae-11eb-22db-3daaf156824a
md"""
### **Sigma**
`σ` is the `residual`-standard deviation listed under the variance components.

"""

# ╔═╡ f4ff8cb0-8cae-11eb-2e53-919c341e56a4
md"""

### **Theta**
The meaning of `θ` is a bit less intuitive. In a less complex model (one that only has intercepts for the random effects) or if we suppress the correlations in the formula with `zerocorr()` then `θ` describes the relationship between the random effects standard deviation and the standard deviation of the residual term.
In our `kb07_m` example:
The `residual` standard deviation is `680.032`.
The standard deviation of our first variance component *`item - (Intercept)`* is `364.713`.
Thus our first `θ` is the relationship: variance component devided by `residual` standard deviation
364.713 /  680.032 =  `0.53631`
"""

# ╔═╡ f4cb5df0-8cae-11eb-1f0f-c3dbe99ea6e1
md"""

We also can calculate the `θ` for variance component *`subj - (Intercept)`*.
The `residual` standard deviation is `680.032`.
The standard deviation of our variance component *`subj - (Intercept)`* is `298.026`.
Thus, the related θ is the relationship: variance component devided by `residual` standard deviation
298.026 /  680.032 =  `0.438252`
"""

# ╔═╡ f496b1d6-8cae-11eb-1b3d-818efd566e11
md"""

We can not calculate the `θ`s for variance component *`item - prec: maintain`* this way, because it includes the correlation of
*`item - prec: maintain`* and *`item - (Intercept)`*. But keep in mind that the relation of  *`item - prec: maintain`*-variability (`252.521`)
and the `residual`-variability (`680.032`) is 252.521  /  680.032 =  `0.3713369`.

The `θ` vector is the flattened version of the variance-covariance matrix - a lowertrinangular matrix.
The on-diagonal elements are just the standard deviations (the `σ`'s). If all off-diagonal elements are zero, we can use our
calculation above. The off-diagonal elements are covariances and correspond to the correlations (the `ρ`'s).
If they are unequal to zero, as it is in our `kb07`-dataset, one way to get the two missing `θ`-values is to take the values directly from the model we have already fitted.

See the two inner values:

"""

# ╔═╡ be62068c-8caf-11eb-0700-d9517d0ab339
md"""

Another way is to make use of the `create_re()` function.
Here you have to define the relation of all random effects variabilities to the variability of the residuals, as shown above,
and the correlation-matrices.

Let's start by defining the correlation matrix for the `item`-part.

The diagonal is always `1.0` because everything is perfectly correlated with itself.
The elements below the diagonal follow the same form as the `Corr.` entries in the output of `VarCorr()`. 
In our example the correlation of
*`item - prec: maintain`* and *`item - (Intercept)`* is `-0.7`.
The elements above the diagonal are just a mirror image.

"""

# ╔═╡ be2fd36a-8caf-11eb-32ca-89721b1a4af9
md"""

Now we put together all relations of standard deviations and the correlation-matrix for the `item`-group.
This calculates the covariance factorization which is the theta matrix.
"""

# ╔═╡ bdfcee32-8caf-11eb-3976-89e317980665
md"""

![](ThetaExplanation.png)


Note: Don't be too specific with your values in create_re(). If there are rounding errors, you will get the error-message: 
`PosDefException: matrix is not Hermitian; Cholesky factorization failed.`

But you can extract the exact values like shown below:

"""

# ╔═╡ bdcc1098-8caf-11eb-0a59-913d4469624a
md"""
Let's continue with the `subj`-part.

We haven't had any correlation, because we only have one intercept term, but you can set it to `1.0`.

"""

# ╔═╡ bdb53588-8caf-11eb-3559-f93ab518b372
re_subj_corr = [1.0]

# ╔═╡ bd9cb5d0-8caf-11eb-1b71-33a671d45eda
md"""

Now we put together all relations of standard deviations and the correlation-matrix for the `subj`-group:

This calculates the covariance factorization which is the theta matrix.

"""

# ╔═╡ bd82ac6c-8caf-11eb-0483-ad95e5567d7c
re_subj = create_re(0.438; corrmat = re_subj_corr)

# ╔═╡ bd6b7cb6-8caf-11eb-0394-713880cf2005
md"""
If you want the exact value you can use

"""

# ╔═╡ bd50de9e-8caf-11eb-33e2-f1f851657ab6
σ_residuals_exact = VarCorr(kb07_m).s
σ_3_exact = VarCorr(kb07_m).σρ[2][1][1] / σ_residuals_exact
re_subj = create_re(σ_3_exact; corrmat = re_subj_corr)

# ╔═╡ bd350412-8caf-11eb-3c83-3b6085365d56
md"""
As mentioned above `θ` is the compact form of these covariance matrices:

"""

# ╔═╡ cc13dd14-8caf-11eb-31fd-6b157205a14f
md"""
We can install these parameter in the `parametricbootstrap()`-function or in the model like this:
"""

# ╔═╡ cbdf51f2-8caf-11eb-0f7b-233b77883cab
md"""

## *A simple example*

Having this knowledge about the parameters we can now **simulate data from scratch**

The `simdat_crossed()` function from `MixedModelsSim` lets you set up a data frame with a specified experimental design.
For now, it only makes fully balanced crossed designs!, but you can generate an unbalanced design by simulating data for the largest cell and deleting extra rows.


Firstly we will set an easy design where `subj_n` subjects per `age` group (O or Y) respond to `item_n` items in each of two `condition`s (A or B).

Your factors need to be specified separately for between-subject, between-item, and within-subject/item factors using `Dict` with the name of each factor as the keys and vectors with the names of the levels as values.

We start with the between subject factors:
"""

# ╔═╡ cbc406cc-8caf-11eb-2e3e-05bf2933c51c
subj_btwn = Dict("age" => ["O", "Y"])

# ╔═╡ cbab8ec8-8caf-11eb-272c-73ccc1889895
md"""
There are no between-item factors in this design so you can omit it or set it to nothing.
Note that if you have factors which are between subject and between item, you need to put them in both dicts.

"""

# ╔═╡ bd192666-8caf-11eb-00fc-9f1b94efad9e
item_btwn = nothing

# ╔═╡ d2db85a6-8cb0-11eb-1cc6-a319ced77e50
md"""
Next, we put within-subject/item factors in a dict:
"""

# ╔═╡ d2bf03cc-8cb0-11eb-3033-997c3653c8ea
both_win = Dict("condition" => ["A", "B"])

# ╔═╡ d2a3e7fe-8cb0-11eb-34a6-cf056d3b9245
md"""
Define subject and item number:
"""

# ╔═╡ d289eef6-8cb0-11eb-35aa-cb4a31363282
subj_n = 10
item_n = 30

# ╔═╡ d26df392-8cb0-11eb-2d29-47ea195a090e
md"""
Simulate data:
"""

# ╔═╡ d254f4b4-8cb0-11eb-05a5-73adb8d66fe9
dat = simdat_crossed(subj_n,
                     item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);

# ╔═╡ d239a768-8cb0-11eb-2ad8-ebae77b9560e
md"""
Have a look:
"""

# ╔═╡ d221998e-8cb0-11eb-384a-e58e153ee139
first(DataFrame(dat),8)

# ╔═╡ d20e0784-8cb0-11eb-3ba7-2d2d2a3d919f
md"""
The values we see in the column `dv` is just random noise.

Set contrasts:

"""

# ╔═╡ d1dacc5c-8cb0-11eb-29dc-3d1978835633
md"""
Define formula:
"""

# ╔═╡ d1c4a04e-8cb0-11eb-04f8-dbfaede3de19
f1 = @formula dv ~ 1 + age * condition + (1|item) + (1|subj);

# ╔═╡ d1ac1916-8cb0-11eb-2f16-774ad7b6eab1
md"""
Note that we did not include condition as random slopes for item and subject.
This is mainly to keep the example simple and to keep the parameter `θ` easier to understand (see Section 3 for the explanation of theta).


Fit the model:
"""

# ╔═╡ d1933d42-8cb0-11eb-0c3a-05d7fc631088
m1 = fit(MixedModel, f1, dat, contrasts=contrasts);
print(m1)

# ╔═╡ d17cbd6a-8cb0-11eb-3dd3-db7a60f2a2ed
md"""

Because the `dv` is just random noise from N(0,1), there will be basically no subject or item random variance, residual variance will be near 1.0, and the estimates for all effects should be small.
Don't worry, we'll specify fixed and random effects directly in `parametricbootstrap()`.


Set random seed for reproducibility:
"""

# ╔═╡ d14960be-8cb0-11eb-0840-f9d41ecdb349
md"""
Specify `β`, `σ`, and `θ`, we just made up this parameter:

"""

# ╔═╡ eca35f9a-8cb0-11eb-0fc8-bfc8d863a61f
new_beta = [0., 0.25, 0.25, 0.]
new_sigma = 2.0
new_theta = [1.0, 1.0]


# ╔═╡ ec88a11e-8cb0-11eb-370f-95b2eb3d6323
md"""
Run nsims iterations:

"""

# ╔═╡ ec421b18-8cb0-11eb-0364-83385e09a16b
md"""
### Power calculation
"""

# ╔═╡ ec13bb76-8cb0-11eb-09b4-252f6e72bf0e
md"""
For nicely displaying it, you can use pretty_table:
"""

# ╔═╡ ebda2058-8cb0-11eb-0b69-2b0602d80b5e
md"""
# 4. Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.

# *4.1 Recreate the `kb07`-dataset from scratch*
For full control over all parameters in our `kb07` data set we will recreate the design using the method shown above.


Define subject and item number:
"""

# ╔═╡ ebc60148-8cb0-11eb-08fa-19c0146dc551
subj_n = 56
item_n = 32

# ╔═╡ eba7356c-8cb0-11eb-316b-218eac0d8ef4
md"""
Define factors in a dict:
"""

# ╔═╡ eb8b3204-8cb0-11eb-3022-3b661160c116
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);

# ╔═╡ eb71936c-8cb0-11eb-24a1-d1d7c35a2ede
md"""
### Play with simdat_crossed

"""

# ╔═╡ 643101b4-8cb1-11eb-195f-db166a2c91f3
subj_btwn = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);
item_btwn = nothing
both_win = nothing;

subj_n = 56
item_n = 32

fake_kb07 = simdat_crossed(subj_n, item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);


fake_kb07_df = DataFrame(fake_kb07)

# ╔═╡ 63fc9674-8cb1-11eb-39e8-499f019f122a
md"""

Simulate data:
"""

# ╔═╡ 63e1bd18-8cb1-11eb-139d-3bb39de43442
fake_kb07 = simdat_crossed(subj_n, item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);

# ╔═╡ 63bf3612-8cb1-11eb-378b-75cf2584ae40
md"""
Make a dataframe:

"""

# ╔═╡ 63668f08-8cb1-11eb-13da-e160a5481fc1
md"""
Have a look:
"""

# ╔═╡ 6324ad2c-8cb1-11eb-08c6-d785cf336cb0
md"""

The function `simdat_crossed` generates a balanced fully crossed design.
Unfortunately, our original design is not fully crossed. Every subject saw an image only once, thus in one of eight possible conditions. To simulate that we only keep one of every eight lines.

We sort the dataframe to enable easier selection
"""

# ╔═╡ 634422d8-8cb1-11eb-2c90-75e0bcc32073
first(fake_kb07_df,8)

# ╔═╡ 62eb8aae-8cb1-11eb-09d9-fffd19330418
md"""

In order to select only the relevant rows of the data set we define an index which represets a random choice of one of every eight rows. First we generate a vector `idx` which represents which row to keep in each set of 8.
"""

# ╔═╡ 62cd000e-8cb1-11eb-2781-1fcbc178c761
len = convert(Int64,(length(fake_kb07)/8))
idx = rand(rng, 1:8 , len)

# ╔═╡ 62acd414-8cb1-11eb-33e9-77e17de7c261
md"""
Then we create an array `A`, of the same length that is populated multiples of the number 8. Added together `A` and `idx` give the indexes of one row from each set of 8s.
"""

# ╔═╡ 628e7604-8cb1-11eb-3df0-1d6c2934cd2d
A = repeat([8], inner=len-1)
A = append!( [0], A )
A = cumsum(A)
idx = idx+A

# ╔═╡ 6253cec8-8cb1-11eb-0465-eb899c3f9fcd
md"""
Reduce the fully crossed design to the original experimental design:
"""

# ╔═╡ 62358b22-8cb1-11eb-2842-d5edfb9c930a
fake_kb07_df= fake_kb07_df[idx, :]
rename!(fake_kb07_df, :dv => :rt_trunc)

# ╔═╡ 6218f3a2-8cb1-11eb-0da5-7b2da2de25a9
md"""
Now we can use the simulated data in the same way as above.

Set contrasts:
"""

# ╔═╡ 61cd0866-8cb1-11eb-2e9d-1363dd663697
md"""
Define formula, same as above:
"""

# ╔═╡ 61947500-8cb1-11eb-3200-b3fcebe16bc7
md"""
Fit the model:
"""

# ╔═╡ 877879d8-8cb1-11eb-3b15-6bc4912abd43
fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df, contrasts=contrasts);
print(fake_kb07_m)

# ╔═╡ 875b0b58-8cb1-11eb-3b7b-75b80a443da9
md"""
Set random seed for reproducibility:
"""

# ╔═╡ 86efe7ec-8cb1-11eb-005f-d94414fa275f
md"""
Then, again, we specify `β`, `σ`, and `θ`.
Here we use the values that we found in the model of the existing dataset:
"""

# ╔═╡ 86b64694-8cb1-11eb-3336-e72f6c8af905
#beta
new_beta = [2181.85, 67.879, -333.791, 78.5904] #manual
new_beta = kb07_m.β #grab from existing model  

#sigma
new_sigma = 680.032 #manual
new_sigma = kb07_m.σ #grab from existing model

#theta
re_item_corr = [1.0 -0.7; -0.7 1.0]
re_item = create_re(0.536, 0.371; corrmat = re_item_corr)
re_subj_corr = [1.0]
re_subj = create_re(0.438; corrmat = re_subj_corr)
new_theta = vcat( flatlowertri(re_item), flatlowertri(re_subj) )  #manual
new_theta = kb07_m.θ #grab from existing model

# ╔═╡ 86a02c18-8cb1-11eb-3589-31916dc01f6e
md"""
Run nsims iterations:
"""

# ╔═╡ 864b574c-8cb1-11eb-06ca-01f99d845a60
md"""
### Power calculation
"""

# ╔═╡ 8c0a0214-8cb1-11eb-2c21-437c819bc6bd
md"""
Compare to the powertable from the existing data:
"""

# ╔═╡ 8bb7f0c8-8cb1-11eb-351f-8f44d0f30dcb
md"""
We have successfully recreated the power simulation of an existing dataset from scratch. This has the advantage, that we now can iterate over different numbers of subjects and items.
"""

# ╔═╡ 8b854e70-8cb1-11eb-0af5-1d4c6e133d4b
md"""
# *4.2 Loop over subject and item sizes*

When designing a study, you may be interested in trying various numbers of subjects and items to see how that affects the power of your study.
To do this you can use a loop to run simulations over a range of values for any parameter.

### We first define every fixed things outside the loop (same as above):

Define factors in a dict:
"""

# ╔═╡ 8b67b7ac-8cb1-11eb-30c2-518e707ccfbc
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);

# ╔═╡ 862bc706-8cb1-11eb-3f99-b7721ef04165
md"""
Set contrasts:
"""

# ╔═╡ 598a9758-8cb2-11eb-2b4e-f9f2e86f4ba1
md"""
Set random seed for reproducibility:
"""

# ╔═╡ 59159c28-8cb2-11eb-3665-1164efea81ce
md"""
Define formula:
"""

# ╔═╡ 1031aca4-8cad-11eb-0a14-e1595ba754b6
kb07_m = fit(MixedModel, kb07_f, kb07; contrasts=contrasts)

# ╔═╡ 45870a2a-8cad-11eb-0fd0-fd70946df45c
kb07_sim = parametricbootstrap(rng, nsims, kb07_m; use_threads = false);

# ╔═╡ 890fa504-8cad-11eb-247c-955a9b083616
begin
	df = DataFrame(kb07_sim.allpars);
	first(df, 9)
	nrow(df)
end

# ╔═╡ ad317520-8cad-11eb-3eb5-45e72a39623f
begin
	σres = @where(df, :type .== "σ", :group .== "residual").value
	plot(x = σres, Geom.density, Guide.xlabel("Parametric bootstrap estimates of σ"), Guide.ylabel("Density"))
end

# ╔═╡ e0584840-8cad-11eb-27e4-e56c7356ef1d
begin
	βInt = @where(df, :type .== "β", :names .== "(Intercept)").value
	plot(x = βInt, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β 	(Intercept)"), Guide.ylabel("Density"))
end

# ╔═╡ fa13cb4a-8cad-11eb-3812-21ca536d38db
begin
	βSpeaker = @where(df, :type .== "β", :names .== "spkr: old").value
	plot(x = βSpeaker, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Speaker"), Guide.ylabel("Density"))
end 

# ╔═╡ 04b24bd0-8cae-11eb-1281-31727cd028fb
begin
	βPrecedents = @where(df, :type .== "β", :names .== "prec: maintain").value
	plot(x = βPrecedents, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Precedents"), Guide.ylabel("Density"))
end

# ╔═╡ 0e635c8c-8cae-11eb-3d39-4fe9d2c34813
begin
	βLoad = @where(df, :type .== "β", :names .== "load: yes").value
	plot(x = βLoad, Geom.density, Guide.xlabel("Parametric bootstrap estimates of β Load"), Guide.ylabel("Density"))
end

# ╔═╡ 30978c4c-8cae-11eb-0506-299579d1686a
kb07_sim_df = DataFrame(kb07_sim.coefpvalues);

# ╔═╡ 4b9f5de4-8cae-11eb-3bef-6f25e1e2d0d0
first(kb07_sim_df, 8)

# ╔═╡ 6cfdc3cc-8cae-11eb-111a-37898255fdaf
	kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:]

# ╔═╡ df391c0c-8cae-11eb-2565-135c8743b69f
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05)

# ╔═╡ 2c1beefc-8cb2-11eb-019c-c9cecfb8fcee
power_table(kb07_sim)

# ╔═╡ f5da0f5c-8cae-11eb-3a4c-d90fec6aa543
new_beta = kb07_m.β

# ╔═╡ ec5abd62-8cb0-11eb-3913-7f9a344f7da6
sim1 = parametricbootstrap(rng, nsims, m1,
                        β = new_beta,
                        σ = new_sigma,
                        θ = new_theta,
                        use_threads = false);

# ╔═╡ 6f10bce6-8cae-11eb-33d5-abf752392c44
pretty_table(ptbl)

# ╔═╡ ebf7ef20-8cb0-11eb-1962-f5fbb7c18a10
pretty_table(ptbl)

# ╔═╡ 86671d56-8cb1-11eb-0a44-59c8797fc9f1
fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m,
                        β = new_beta,
                        σ = new_sigma,
                        θ = new_theta,
                        use_threads = false);

# ╔═╡ 8c20fc4e-8cb1-11eb-0626-61fc6029241c
power_table(fake_kb07_sim)

# ╔═╡ 8bd460a0-8cb1-11eb-1e3a-0b666f276101
power_table(fake_kb07_sim)

# ╔═╡ 6c81882e-8caf-11eb-236a-634ffba3187d
new_beta[2:4] = kb07_m.β[2:4]/2

# ╔═╡ f5a53fc0-8cae-11eb-227b-adeae8b329f4
kb07_sim_half = parametricbootstrap(rng, nsims, kb07_m; β = new_beta, use_threads = false);

# ╔═╡ f5755698-8cae-11eb-1b35-2f5808c8a3ef
power_table(kb07_sim_half)

# ╔═╡ f5463ae8-8cae-11eb-20a0-19b54e79984f
kb07_m

# ╔═╡ ae24d704-8caf-11eb-0c12-b9237a676ff4
kb07_m.β

# ╔═╡ f517b1d2-8cae-11eb-1b78-e73c3f529834
kb07_m

# ╔═╡ b1209682-8caf-11eb-3f3d-e735c6c6b379
kb07_m.σ

# ╔═╡ f4e4a6c0-8cae-11eb-1e9a-67b07d7a8f97
kb07_m.θ

# ╔═╡ f4b1e21c-8cae-11eb-0a84-474375b8992c
kb07_m.θ

# ╔═╡ f47c05f2-8cae-11eb-1f57-c54a459031fc
kb07_m.θ

# ╔═╡ cc25963a-8caf-11eb-209d-032effb04a4c
kb07_m.θ = vcat( flatlowertri(re_item), flatlowertri(re_subj) )

# ╔═╡ cbf7f830-8caf-11eb-35ce-c960d58c6550
update!(kb07_m, re_item, re_subj)

# ╔═╡ 58afbabe-8cb2-11eb-08eb-33aeeb1b99f5
md"""
Specify `β`, `σ`, and `θ`:
"""

# ╔═╡ 58814f0a-8cb2-11eb-2874-3d0cf06e52fb
#beta
new_beta = [2181.85, 67.879, -333.791, 78.5904]
new_beta = kb07_m.β

#sigma
new_sigma = 680.032
new_sigma = kb07_m.σ

#theta
re_item_corr = [1.0 -0.7; -0.7 1.0]
re_item = create_re(0.536, 0.371; corrmat = re_item_corr)
re_subj_corr = [1.0]
re_subj = create_re(0.438; corrmat = re_subj_corr)
new_theta = vcat( flatlowertri(re_item), flatlowertri(re_subj) )
new_theta = kb07_m.θ

# ╔═╡ 581c415a-8cb2-11eb-281e-cb673a54e27e
md"""
### Then we define the variables that out loop will iterate over

Define subject and item numbers as arrays:
"""

# ╔═╡ 57ed1fc4-8cb2-11eb-37d4-05d44823b83d
sub_ns = [20, 30, 40];
item_ns = [16, 24, 32];

# ╔═╡ 57b64e10-8cb2-11eb-0b72-f3b0ddc7f93a
md"""
Make an emty dataframe:
"""

# ╔═╡ 57896b5a-8cb2-11eb-167f-ad146697a564
d = DataFrame();

# ╔═╡ 576d6ac2-8cb2-11eb-2708-c74c10a0bd88
md"""
### Run the loop:
"""

# ╔═╡ 57399788-8cb2-11eb-3954-87e32d726689
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

# ╔═╡ 570c01d8-8cb2-11eb-3752-b3740a41cc41
md"""
Our dataframe `d` now contains the power information for each combination of subjects and items.
"""

# ╔═╡ 56dc62b6-8cb2-11eb-190a-3129097f1e6c
print(d)

# ╔═╡ 56acfabe-8cb2-11eb-1f33-3387291ff9a0
md"""
Lastly we plot our results:
"""

# ╔═╡ 568e26dc-8cb2-11eb-0e2d-c37960cf8ebb
categorical!(d, :item_n)

plot(d, x="subj_n", y="power",xgroup= "coefname",color="item_n", Geom.subplot_grid(Geom.point, Geom.line), Guide.xlabel("Number of subjects by parameter"), Guide.ylabel("Power"))

# ╔═╡ 565d819e-8cb2-11eb-1597-3533294d81dc
md"""
# Credit
This tutorial was conceived for ZiF research and tutorial workshop by Lisa DeBruine (Feb. 2020) presented again by Phillip Alday during the SMLP Summer School (Sep. 2020).

Updated and extended by Lisa Schwetlick & Daniel Backhaus, with the kind help of Phillip Alday, after changes to the package.
"""

# ╔═╡ ec2a8f2a-8cb0-11eb-0f20-5dd450803442
ptbl= power_table(sim1)

# ╔═╡ be178ef4-8caf-11eb-1c7d-99ce46120231
re_item = create_re(0.536, 0.371; corrmat = re_item_corr)

# ╔═╡ 873ab954-8cb1-11eb-1db9-2f95dd1c9086
rng = StableRNG(42);

# ╔═╡ bde81fa2-8caf-11eb-1e73-75ca31dbca9b
begin
	corr_exact = VarCorr(kb07_m).σρ[1][2][1][1]
	σ_residuals_exact = VarCorr(kb07_m).s
	σ_1_exact = VarCorr(kb07_m).σρ[1][1][1] / residuals_exact
	σ_2_exact = VarCorr(kb07_m).σρ[1][1][2] / residuals_exact
	
	re_item_corr = [1.0 corr_exact; corr_exact 1.0]
	re_item = create_re(σ_1_exact, σ_2_exact; corrmat = re_item_corr)
end

# ╔═╡ d1f4be3c-8cb0-11eb-2dfc-9d1af5fcb02e
contrasts = Dict(:age => HelmertCoding(),
                 :condition => HelmertCoding());


# ╔═╡ 6387cb32-8cb1-11eb-2a55-f5f73b024f3f
fake_kb07_df = DataFrame(fake_kb07);

# ╔═╡ 58e30664-8cb2-11eb-336f-bfa2d55e1cb2
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );

# ╔═╡ 595439e2-8cb2-11eb-3808-4f209b1bd3d7
rng = StableRNG(42);

# ╔═╡ e7c3c40a-8cac-11eb-3a5c-a99081421b16
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());

# ╔═╡ 6307f5b0-8cb1-11eb-3dea-a9f1397a826f
fake_kb07_df = sort(fake_kb07_df, [:subj, :item, :load, :prec, :spkr])

# ╔═╡ 37d145da-8cad-11eb-0cbe-a108993570f2
rng = StableRNG(42);

# ╔═╡ 054e7b6e-8cad-11eb-0eb3-73c98860c609
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );

# ╔═╡ be49a9a2-8caf-11eb-00c6-c50e74fc7667
re_item_corr = [1.0 -0.7; -0.7 1.0]

# ╔═╡ 4b5cd0d2-8cae-11eb-313b-d7599e45e40a
ptbl = power_table(kb07_sim, 0.05)

# ╔═╡ d162e674-8cb0-11eb-3a48-e912d191d4ef
rng = StableRNG(42);

# ╔═╡ 61b17fce-8cb1-11eb-070d-9d34832e7190
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );

# ╔═╡ 61eba730-8cb1-11eb-3336-817bf39c41d5
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());

# ╔═╡ 59a940c2-8cb2-11eb-1708-9fab465f5527
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());

# ╔═╡ Cell order:
# ╟─3edca7b8-8cab-11eb-3982-0bfba7ea51b3
# ╠═b8016ebc-8cab-11eb-0bd3-5bea764c4251
# ╟─7cf5bc5a-8cac-11eb-1b2d-25987ca50e72
# ╠═77079142-8cac-11eb-2894-3b1af57b13cd
# ╟─a6254e9c-8cac-11eb-2ecc-8f141b7fb022
# ╠═c2a36428-8cac-11eb-1eba-85a5bba4f980
# ╟─dd95cf64-8cac-11eb-14b4-a91cf59e9094
# ╠═e7c3c40a-8cac-11eb-3a5c-a99081421b16
# ╟─ef62641e-8cac-11eb-1b85-7be4b0eea3ea
# ╠═054e7b6e-8cad-11eb-0eb3-73c98860c609
# ╟─0d827ad8-8cad-11eb-2c4c-9348874e0343
# ╠═1031aca4-8cad-11eb-0a14-e1595ba754b6
# ╟─2702a6d6-8cad-11eb-3922-5fa14ce5c7b8
# ╠═37d145da-8cad-11eb-0cbe-a108993570f2
# ╟─3eceee46-8cad-11eb-2229-fdd1dc488c31
# ╠═45870a2a-8cad-11eb-0fd0-fd70946df45c
# ╟─7d3b3aae-8cad-11eb-209a-8d5a4b726ec6
# ╠═890fa504-8cad-11eb-247c-955a9b083616
# ╟─9b7d0cf4-8cad-11eb-3ecc-d3e4e6287837
# ╠═ad317520-8cad-11eb-3eb5-45e72a39623f
# ╠═e0584840-8cad-11eb-27e4-e56c7356ef1d
# ╠═fa13cb4a-8cad-11eb-3812-21ca536d38db
# ╠═04b24bd0-8cae-11eb-1281-31727cd028fb
# ╠═0e635c8c-8cae-11eb-3d39-4fe9d2c34813
# ╟─165c3f30-8cae-11eb-037a-f337833479cf
# ╠═30978c4c-8cae-11eb-0506-299579d1686a
# ╟─4be8e2e8-8cae-11eb-2843-c79f6cbad261
# ╠═4b9f5de4-8cae-11eb-3bef-6f25e1e2d0d0
# ╟─4b7d4632-8cae-11eb-236e-3da3b6e9d9fb
# ╠═4b5cd0d2-8cae-11eb-313b-d7599e45e40a
# ╟─4ad40e64-8cae-11eb-0304-fbc12a867ef1
# ╠═6cfdc3cc-8cae-11eb-111a-37898255fdaf
# ╠═df391c0c-8cae-11eb-2565-135c8743b69f
# ╟─6cd06d82-8cae-11eb-0d5b-bdf66b8ddd86
# ╠═6f10bce6-8cae-11eb-33d5-abf752392c44
# ╟─4683428a-8cae-11eb-2cdb-3516864105ce
# ╠═f60e9542-8cae-11eb-3708-7d07915f3c73
# ╟─f5f3bd44-8cae-11eb-06b5-b330398ec687
# ╠═f5da0f5c-8cae-11eb-3a4c-d90fec6aa543
# ╠═6c81882e-8caf-11eb-236a-634ffba3187d
# ╟─f5c0d794-8cae-11eb-13a1-af079fabc06c
# ╠═f5a53fc0-8cae-11eb-227b-adeae8b329f4
# ╟─f590368e-8cae-11eb-00cc-2142ff14576c
# ╠═f5755698-8cae-11eb-1b35-2f5808c8a3ef
# ╟─f55dd39c-8cae-11eb-17e2-5b4e85500ded
# ╠═f5463ae8-8cae-11eb-20a0-19b54e79984f
# ╠═ae24d704-8caf-11eb-0c12-b9237a676ff4
# ╟─f52d4376-8cae-11eb-22db-3daaf156824a
# ╠═f517b1d2-8cae-11eb-1b78-e73c3f529834
# ╠═b1209682-8caf-11eb-3f3d-e735c6c6b379
# ╟─f4ff8cb0-8cae-11eb-2e53-919c341e56a4
# ╠═f4e4a6c0-8cae-11eb-1e9a-67b07d7a8f97
# ╟─f4cb5df0-8cae-11eb-1f0f-c3dbe99ea6e1
# ╠═f4b1e21c-8cae-11eb-0a84-474375b8992c
# ╟─f496b1d6-8cae-11eb-1b3d-818efd566e11
# ╠═f47c05f2-8cae-11eb-1f57-c54a459031fc
# ╟─be62068c-8caf-11eb-0700-d9517d0ab339
# ╠═be49a9a2-8caf-11eb-00c6-c50e74fc7667
# ╟─be2fd36a-8caf-11eb-32ca-89721b1a4af9
# ╠═be178ef4-8caf-11eb-1c7d-99ce46120231
# ╟─bdfcee32-8caf-11eb-3976-89e317980665
# ╠═bde81fa2-8caf-11eb-1e73-75ca31dbca9b
# ╟─bdcc1098-8caf-11eb-0a59-913d4469624a
# ╠═bdb53588-8caf-11eb-3559-f93ab518b372
# ╟─bd9cb5d0-8caf-11eb-1b71-33a671d45eda
# ╠═bd82ac6c-8caf-11eb-0483-ad95e5567d7c
# ╟─bd6b7cb6-8caf-11eb-0394-713880cf2005
# ╠═bd50de9e-8caf-11eb-33e2-f1f851657ab6
# ╟─bd350412-8caf-11eb-3c83-3b6085365d56
# ╠═cc25963a-8caf-11eb-209d-032effb04a4c
# ╟─cc13dd14-8caf-11eb-31fd-6b157205a14f
# ╠═cbf7f830-8caf-11eb-35ce-c960d58c6550
# ╟─cbdf51f2-8caf-11eb-0f7b-233b77883cab
# ╠═cbc406cc-8caf-11eb-2e3e-05bf2933c51c
# ╟─cbab8ec8-8caf-11eb-272c-73ccc1889895
# ╠═bd192666-8caf-11eb-00fc-9f1b94efad9e
# ╠═d2db85a6-8cb0-11eb-1cc6-a319ced77e50
# ╠═d2bf03cc-8cb0-11eb-3033-997c3653c8ea
# ╠═d2a3e7fe-8cb0-11eb-34a6-cf056d3b9245
# ╠═d289eef6-8cb0-11eb-35aa-cb4a31363282
# ╠═d26df392-8cb0-11eb-2d29-47ea195a090e
# ╠═d254f4b4-8cb0-11eb-05a5-73adb8d66fe9
# ╠═d239a768-8cb0-11eb-2ad8-ebae77b9560e
# ╠═d221998e-8cb0-11eb-384a-e58e153ee139
# ╠═d20e0784-8cb0-11eb-3ba7-2d2d2a3d919f
# ╠═d1f4be3c-8cb0-11eb-2dfc-9d1af5fcb02e
# ╠═d1dacc5c-8cb0-11eb-29dc-3d1978835633
# ╠═d1c4a04e-8cb0-11eb-04f8-dbfaede3de19
# ╠═d1ac1916-8cb0-11eb-2f16-774ad7b6eab1
# ╠═d1933d42-8cb0-11eb-0c3a-05d7fc631088
# ╠═d17cbd6a-8cb0-11eb-3dd3-db7a60f2a2ed
# ╠═d162e674-8cb0-11eb-3a48-e912d191d4ef
# ╠═d14960be-8cb0-11eb-0840-f9d41ecdb349
# ╠═eca35f9a-8cb0-11eb-0fc8-bfc8d863a61f
# ╠═ec88a11e-8cb0-11eb-370f-95b2eb3d6323
# ╠═ec5abd62-8cb0-11eb-3913-7f9a344f7da6
# ╠═ec421b18-8cb0-11eb-0364-83385e09a16b
# ╠═ec2a8f2a-8cb0-11eb-0f20-5dd450803442
# ╠═ec13bb76-8cb0-11eb-09b4-252f6e72bf0e
# ╠═ebf7ef20-8cb0-11eb-1962-f5fbb7c18a10
# ╠═ebda2058-8cb0-11eb-0b69-2b0602d80b5e
# ╠═ebc60148-8cb0-11eb-08fa-19c0146dc551
# ╠═eba7356c-8cb0-11eb-316b-218eac0d8ef4
# ╠═eb8b3204-8cb0-11eb-3022-3b661160c116
# ╠═eb71936c-8cb0-11eb-24a1-d1d7c35a2ede
# ╠═643101b4-8cb1-11eb-195f-db166a2c91f3
# ╠═63fc9674-8cb1-11eb-39e8-499f019f122a
# ╠═63e1bd18-8cb1-11eb-139d-3bb39de43442
# ╠═63bf3612-8cb1-11eb-378b-75cf2584ae40
# ╠═6387cb32-8cb1-11eb-2a55-f5f73b024f3f
# ╠═63668f08-8cb1-11eb-13da-e160a5481fc1
# ╠═634422d8-8cb1-11eb-2c90-75e0bcc32073
# ╠═6324ad2c-8cb1-11eb-08c6-d785cf336cb0
# ╠═6307f5b0-8cb1-11eb-3dea-a9f1397a826f
# ╠═62eb8aae-8cb1-11eb-09d9-fffd19330418
# ╠═62cd000e-8cb1-11eb-2781-1fcbc178c761
# ╠═62acd414-8cb1-11eb-33e9-77e17de7c261
# ╠═628e7604-8cb1-11eb-3df0-1d6c2934cd2d
# ╠═6253cec8-8cb1-11eb-0465-eb899c3f9fcd
# ╠═62358b22-8cb1-11eb-2842-d5edfb9c930a
# ╠═6218f3a2-8cb1-11eb-0da5-7b2da2de25a9
# ╠═61eba730-8cb1-11eb-3336-817bf39c41d5
# ╠═61cd0866-8cb1-11eb-2e9d-1363dd663697
# ╠═61b17fce-8cb1-11eb-070d-9d34832e7190
# ╠═61947500-8cb1-11eb-3200-b3fcebe16bc7
# ╠═877879d8-8cb1-11eb-3b15-6bc4912abd43
# ╠═875b0b58-8cb1-11eb-3b7b-75b80a443da9
# ╠═873ab954-8cb1-11eb-1db9-2f95dd1c9086
# ╠═86efe7ec-8cb1-11eb-005f-d94414fa275f
# ╠═86b64694-8cb1-11eb-3336-e72f6c8af905
# ╠═86a02c18-8cb1-11eb-3589-31916dc01f6e
# ╠═86671d56-8cb1-11eb-0a44-59c8797fc9f1
# ╠═864b574c-8cb1-11eb-06ca-01f99d845a60
# ╠═8c20fc4e-8cb1-11eb-0626-61fc6029241c
# ╠═8c0a0214-8cb1-11eb-2c21-437c819bc6bd
# ╠═8bd460a0-8cb1-11eb-1e3a-0b666f276101
# ╠═2c1beefc-8cb2-11eb-019c-c9cecfb8fcee
# ╠═8bb7f0c8-8cb1-11eb-351f-8f44d0f30dcb
# ╠═8b854e70-8cb1-11eb-0af5-1d4c6e133d4b
# ╠═8b67b7ac-8cb1-11eb-30c2-518e707ccfbc
# ╠═862bc706-8cb1-11eb-3f99-b7721ef04165
# ╠═59a940c2-8cb2-11eb-1708-9fab465f5527
# ╠═598a9758-8cb2-11eb-2b4e-f9f2e86f4ba1
# ╠═595439e2-8cb2-11eb-3808-4f209b1bd3d7
# ╠═59159c28-8cb2-11eb-3665-1164efea81ce
# ╠═58e30664-8cb2-11eb-336f-bfa2d55e1cb2
# ╠═58afbabe-8cb2-11eb-08eb-33aeeb1b99f5
# ╠═58814f0a-8cb2-11eb-2874-3d0cf06e52fb
# ╠═581c415a-8cb2-11eb-281e-cb673a54e27e
# ╠═57ed1fc4-8cb2-11eb-37d4-05d44823b83d
# ╠═57b64e10-8cb2-11eb-0b72-f3b0ddc7f93a
# ╠═57896b5a-8cb2-11eb-167f-ad146697a564
# ╠═576d6ac2-8cb2-11eb-2708-c74c10a0bd88
# ╠═57399788-8cb2-11eb-3954-87e32d726689
# ╠═570c01d8-8cb2-11eb-3752-b3740a41cc41
# ╠═56dc62b6-8cb2-11eb-190a-3129097f1e6c
# ╠═56acfabe-8cb2-11eb-1f33-3387291ff9a0
# ╠═568e26dc-8cb2-11eb-0e2d-c37960cf8ebb
# ╠═565d819e-8cb2-11eb-1597-3533294d81dc
