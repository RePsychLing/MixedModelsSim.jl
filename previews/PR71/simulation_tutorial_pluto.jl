### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 89495294-10ed-11ec-04a9-db58dba9f3d1
begin
	using MixedModels        # run mixed models
	using MixedModelsSim     # simulation functions for mixed models
	using DataFrames, Tables # work with data tables
	using StableRNGs         # random number generator
	using CSV                # write CSV files
	using Statistics         # basic math functions
	using DataFramesMeta     # dplyr-like operations
	using Gadfly             # plotting package
	using PlutoUI
end

# ╔═╡ c071c03f-6ddb-415e-8445-b4fd8d45abc1
md"""
Power Analysis and Simulation Tutorial
======================================

This tutorial demonstrates how to conduct power analyses and data simulation using Julia and the MixedModelsSim package.

Power analysis is an important tool for planning an experimental design. Here we show how to

1. Take existing data and calculate power by simulating new data.
2. Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set.
3. Create a (simple) balanced fully crossed dataset from scratch and analyze power.
4. Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.

### Load the packages we'll be using in Julia

First, here are the packages needed in this example.
"""

# ╔═╡ 14ae0269-ea23-4c73-a04f-40eece7a5bc5
TableOfContents()

# ╔═╡ 76df9a69-7ee6-4876-9dcc-e6178f59d3ad
md"
### Define number of iterations

Here we define how many model simulations we want to do. A large number will give more reliable results, but will take longer to compute. It is useful to set it to a low number for testing, and increase it for your final analysis.
"

# ╔═╡ 16822577-6f44-4330-9181-678bb2c7cce2
nsims = 500

# ╔═╡ d4e9bb7e-3f70-48fb-a73c-9606cfb39f39
md"""
##  Take existing data and calculate power by simulate new data with bootstrapping.

### Build a Linear Mixed Model from existing data.

For the first example we are going to simulate bootstrapped data from an existing data set:

*Experiment 2 from Kronmüller, E., & Barr, D. J. (2007). Perspective-free pragmatics: Broken precedents and the recovery-from-preemption hypothesis. Journal of Memory and Language, 56(3), 436-455.*

The data we will be using through out this tutorial is a study about how in a conversation the change of a speaker or the change of precedents (which are patterns of word usage to describe an object, e.g. one can refer to the same object "white shoes", "runners", "sneakers") affects the understanding.

Objects are presented on a screen while participants listen to instructions to move the objects around. Participants eye movements are tracked.
The dependent variable is response time, defined as the latency between the onset of the test description and the moment at which the target was selected.
The independent variables are speaker (old vs. new), precedents (maintain vs. break) and cognitive load (a secondary memory task).

We have to load the data and define some characteristics like the contrasts and the underlying model.

Load existing data:
"""

# ╔═╡ 4c00bf63-1d33-4f16-9776-e2f3915850f9
kb07 = MixedModels.dataset(:kb07);

# ╔═╡ 87af988a-2d55-4c24-a132-3c76f1ae562b
md"""
Set contrasts:
"""

# ╔═╡ 8caad72c-a8d7-4e23-8eb6-b9096141b471
contrasts = Dict(:spkr => HelmertCoding(),
                 :prec => HelmertCoding(),
                 :load => HelmertCoding());

# ╔═╡ 29b32ea9-0018-4ec3-9095-1e527e11fc4f
md"""
The chosen LMM for this dataset is defined by the following model formula:
"""

# ╔═╡ 3763460c-046d-46ce-b9ef-4faf1fbd0428
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );

# ╔═╡ f05a1f73-503a-4ed7-8879-1da242b2b224
md"""
Fit the model:
"""

# ╔═╡ a75ad071-0789-41ab-9253-27c8bca615ad
kb07_m = fit(MixedModel, kb07_f, kb07; contrasts=contrasts)

# ╔═╡ 8e4629aa-58b7-4ab8-8954-74ee579a3f9b
md"""
### Simulate from existing data with same model parameters

We will first look at the power of the dataset with the same parameters as in the original data set. This means that each dataset will have the exact number of observations as the original data. Here, we use the model `kb07_m` we fitted above to our dataset `kb07`.

You can use the `parameparametricbootstrap()` function to run `nsims` iterations of data sampled using the parameters from `kb07_m`.
Set up a random seed to make the simulation reproducible. You can use your favourite number.

To use multithreading, you need to set the number of cores you want to use.
E.g. in Visual Studio Code, open the settings (gear icon in the lower left corner or cmd-,) and search for "thread".
Set `julia.NumThreads` to the number of cores you want to use (at least 1 less than your total number).


Set random seed for reproducibility:
"""

# ╔═╡ 3e0e9908-002c-4def-a2eb-c8f7d775f205
rng = StableRNG(42);

# ╔═╡ 33d9ecae-ac5c-4fdf-bd1a-4d68f785084a
md"""
Run nsims iterations:
"""

# ╔═╡ b3e2dea0-bfb0-4d5d-bd9c-8ed0e0272cc0
kb07_sim = parametricbootstrap(rng, nsims, kb07_m; use_threads = false);

# ╔═╡ 1b814e14-136e-45ad-83e9-c783e30dc4b1
md"""
**Try**: Run the code above with or without `use_threads = true`.

The output DataFrame `kb07_sim` contains the results of the bootstrapping procedure.
"""

# ╔═╡ fc141d66-3252-4989-b6e4-98bbfe82385a
df = DataFrame(kb07_sim.allpars);

# ╔═╡ 5d7c1822-453c-4483-a1a1-32e149145afb
first(df, 9)

# ╔═╡ ef3c9d81-38cf-45d6-9258-926ad209e36f
nrow(df)

# ╔═╡ 146c6a52-a512-4601-9b58-bc56575bfb2e
md"""
The dataframe df has 4500 rows: 9 parameters, each from 500 iterations.

Plot some bootstrapped parameter:"""

# ╔═╡ 4598ae9a-8d39-4280-a74a-4d7145d64920
begin
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
end

# ╔═╡ fa91a0f8-675e-4155-9ff4-58090689c397
md"""
Convert p-values of your fixed-effects parameters to dataframe
"""

# ╔═╡ 811a81ac-50d9-4d07-98d1-3fc8e1586685
kb07_sim_df = DataFrame(kb07_sim.coefpvalues);

# ╔═╡ d8257d21-8dad-4533-ba93-78ba7c3202b3
md"""
Have a look at your simulated data: 
"""

# ╔═╡ ac2f3b62-e86c-43b2-8109-558db79e5773
first(kb07_sim_df, 8)

# ╔═╡ 67493854-920e-477e-805b-21d58cd19ade
md"""
Now that we have a bootstrapped data, we can start our power calculation.

### Power calculation

The function `power_table()` from `MixedModelsSim` takes the output of `parametricbootstrap()` and calculates the proportion of simulations where the p-value is less than alpha for each coefficient.
You can set the `alpha` argument to change the default value of 0.05 (justify your alpha).
"""

# ╔═╡ 8c06724f-a947-44c2-8541-7aa50c915356
ptbl = power_table(kb07_sim,0.05)

# ╔═╡ fc239c6e-ca35-4df6-81b2-115be160437b
md"""
An estimated power of 1 means that in every iteration the specific parameter we are looking at was below our alpha.
An estimated power of 0.5 means that in half of our iterations the specific parameter we are looking at was below our alpha.
An estimated power of 0 means that for none of our iterations the specific parameter we are looking at was below our alpha.

You can also do it manually:
"""

# ╔═╡ 97907ccc-e4b8-43a3-9e5a-ae1da59a203a
kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:]

# ╔═╡ a383994a-e4f4-4c9e-9f17-32ec11e32e87
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("(Intercept)"),:p] .< 0.05)

# ╔═╡ 523c1c21-43d2-469e-b138-3dd5cab1e6b9
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("spkr: old"),:p] .< 0.05)

# ╔═╡ 3419a020-7333-4cb5-981b-84e7fa788280
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05)

# ╔═╡ 231ce20d-ea0c-4314-a7b0-0e296bc2160b
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("load: yes"),:p] .< 0.05)

# ╔═╡ b66ab7be-0547-4535-9d33-cfae62441c61
md"""
For nicely displaying, you can use `pretty_table`:
"""

# ╔═╡ 356751b1-73cc-4dae-ad5a-b1cd0b45ede0
pretty_table(ptbl)

# ╔═╡ 910ca4be-d4bf-4231-95d6-de8b3e97d4d5
md"""
## Adapt parameters in a given Linear Mixed Model to analyze power without changing the existing data set.

Let's say we want to check our power to detect effects of spkr, prec, and load
that are only half the size as in our pilot data. We can set a new vector of beta values
with the `β` argument to `parametricbootstrap()`.


Set random seed for reproducibility:
"""

# ╔═╡ 11eacaa9-f309-4c72-a659-3fc262f3111f
md"""
Specify β:
"""

# ╔═╡ bc6ace0a-2776-409c-892c-6892df79f0e7
new_beta = kb07_m.β

# ╔═╡ bbe93bad-4d41-4300-910a-ec16a660def3
new_beta[2:4] = kb07_m.β[2:4]/2;

# ╔═╡ 9510ce0a-6ba7-4903-8e2a-84f6fc69492c
new_beta

# ╔═╡ f5436bb7-9087-437f-bc43-1b113c907695
md"""Run nsims iterations:"""

# ╔═╡ 9fbfff0f-289c-4c92-9ffa-615ad1b69ff2
kb07_sim_half = parametricbootstrap(rng, nsims, kb07_m; β = new_beta, use_threads = false);

# ╔═╡ 7aa2cea2-9136-42ad-bc62-9966aa93b84c
md"""### Power calculation
"""

# ╔═╡ f7a5cffe-2819-4a9f-a56e-45fefb2c4201
power_table(kb07_sim_half)

# ╔═╡ c3c9ab9c-c81d-44d2-bd45-fb72a7c5da8a
md"""
Convert p-values of your fixed-effects parameters to dataframe
"""

# ╔═╡ 9e9c337e-556f-44a8-ad5f-0476be1458b0
kb07_sim_half_df = DataFrame(kb07_sim_half.coefpvalues);

# ╔═╡ 6b7fdacd-1343-42a2-8e7e-d06bd559765d
mean(kb07_sim_half_df[kb07_sim_half_df.coefname .== Symbol("(Intercept)"),:p] .< 0.05)

# ╔═╡ f232988a-c3cf-404f-9f0d-fe0ef591e8bb
mean(kb07_sim_half_df[kb07_sim_half_df.coefname .== Symbol("spkr: old"),:p] .< 0.05)

# ╔═╡ ed309fad-0058-4069-8e8e-e526b2b36370
mean(kb07_sim_half_df[kb07_sim_half_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05)

# ╔═╡ 2fc97a13-abc5-42ff-a1eb-28e24271fae0
mean(kb07_sim_half_df[kb07_sim_half_df.coefname .== Symbol("load: yes"),:p] .< 0.05)

# ╔═╡ 42f4e755-c7b1-46e2-91fd-31a19ccd5685
md"""## Create a (simple) balanced fully crossed dataset from scratch and analyze power.
"""

# ╔═╡ fcd25a28-cb0d-4077-8b20-f5837c4e99f8
md"""
In some situations, instead of using an existing dataset it may be useful to simulate the data from scratch. This could be the case when the original data is not available but the effect sizes are known. That means that we have to:

a) specify the effect sizes manually

b) manually create an experimental design, according to which data can be simulated

If we simulate data from scratch, aside from subject and item number, we can manipulate the arguments `β`, `σ` and `θ`.
Lets have a closer look at them, define their meaning and we will see where the corresponding values in the model output are.
"""

# ╔═╡ ae7efbb0-aa35-4128-aebf-10a1deab5fee
md"""
### Fixed Effects (βs)
`β` are our effect sizes. If we look again on our LMM summary from the `kb07`-dataset `kb07_m`
we see our four `β` under fixed-effects parameters in the `Est.`-column.
"""

# ╔═╡ 21567edd-0593-4823-b538-15dcb0de48f0
kb07_m

# ╔═╡ 5eb6074a-a695-4404-841e-21b7b82f1150
kb07_m.β

# ╔═╡ a02e6191-239f-4b0b-996d-f73bfc54ea28
md"""
### Residual Variance (σ)
`σ` is the `residual`-standard deviation also listed in the `Est.`-column.
"""

# ╔═╡ ce275097-63e5-410e-a23b-9fdd4250e35c
kb07_m

# ╔═╡ c5f9f6ad-0b18-434b-a368-42afd8cb398e
kb07_m.σ

# ╔═╡ 9adf46f1-d84c-4ae3-af1b-2087adbf0315
md"""
### Random Effects (θ)
The meaning of `θ` is a bit less intuitive. In a less complex model (one that only has intercepts for the random effects) or if we suppress the correlations in the formula with `zerocorr()` then `θ` describes the relationship between the random effects standard deviation and the standard deviation of the residual term.
In our `kb07_m` example:
The `residual` standard deviation is `680.032`.
The standard deviation of our first variance component *`item - (Intercept)`* is `364.713`.
Thus our first `θ` is the relationship: variance component divided by `residual` standard deviation
364.713 /  680.032 =  `0.53631`
"""

# ╔═╡ 077b22ae-814f-4b18-8297-2124317cbff6
kb07_m.θ

# ╔═╡ f33d0b96-3b3c-4651-8b39-46461445cd09
md"""
We also can calculate the `θ` for variance component *`subj - (Intercept)`*.
The `residual` standard deviation is `680.032`.
The standard deviation of our variance component *`subj - (Intercept)`* is `298.026`.
Thus, the related θ is the relationship: variance component divided by `residual` standard deviation
298.026 /  680.032 =  `0.438252`
"""

# ╔═╡ 8103cf6e-18a2-40a1-8043-c87ed92e1c11
kb07_m.θ

# ╔═╡ cdee3c0c-6c85-4515-8ffd-1b1c71018e09
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

# ╔═╡ f5535f46-e919-4dad-9f0e-337d68ab669b
kb07_m.θ

# ╔═╡ 6f2c0659-3f5d-45d8-8a40-b3124ce49fa3
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

# ╔═╡ 9b5601c0-613e-426f-9318-c87ba9108f42
re_item_corr = [1.0 -0.7; -0.7 1.0]

# ╔═╡ d0493428-916e-4bfe-8f5a-586e21af0f5c
md"""
Now we put together all relations of standard deviations and the correlation-matrix for the `item`-group.
This calculates the covariance factorization which is the theta matrix.
"""

# ╔═╡ ca0b404c-50da-43eb-8bb6-d6037e09ba64
re_item = create_re(0.536, 0.371; corrmat = re_item_corr)

# ╔═╡ 2a04994e-bb88-455b-a732-bae9ef1c5422
md" ![](https://github.com/lschwetlick/MixedModelsSim.jl/blob/new_tutorial/docs/src/ThetaExplanation.png?raw=true)"

# ╔═╡ 5e364a80-e5d9-458b-be1e-d00dafa04033
md"""
Note: Don't be too specific with your values in create_re(). If there are rounding errors, you will get the error-message:
`PosDefException: matrix is not Hermitian; Cholesky factorization failed.`

But you can extract the exact values like shown below:
"""

# ╔═╡ 026c051d-bcfa-4b1d-8648-bdd0d8cdf894
corr_exact = VarCorr(kb07_m).σρ.item.ρ[1]

# ╔═╡ 2951d7fe-0373-4f52-86e3-c11bd7a57f1a
σ_residuals_exact = VarCorr(kb07_m).s

# ╔═╡ 1ddbde42-6dcc-425d-b906-92880bb15e89
σ_1_exact = VarCorr(kb07_m).σρ.item.σ[1] / σ_residuals_exact

# ╔═╡ 29ef3cab-e774-4cf8-a392-92eec8acc98b
σ_2_exact = VarCorr(kb07_m).σρ.item.σ[2] / σ_residuals_exact

# ╔═╡ 37f07a78-3616-41c4-90ad-6fc52d9223b3
re_item_corr_exact = [1.0 corr_exact; corr_exact 1.0]

# ╔═╡ f4929a81-71db-4a78-86f9-22c7ed89f04d
re_item_exact = create_re(σ_1_exact, σ_2_exact; corrmat = re_item_corr_exact)

# ╔═╡ dd8aa406-ef50-4036-9d6f-f499934a84c3
md"""
Let's continue with the `subj`-part.

Since there the by-subject random effects have only one entry (the intercept), there are no correlations to specify and we can omit the `corrmat` argument.

Now we put together all relations of standard deviations and the correlation-matrix for the `subj`-group:

This calculates the covariance factorization which is the theta matrix.
"""

# ╔═╡ f8432a72-da18-4540-9a26-2c4091bd43ac
re_subj = create_re(0.438)

# ╔═╡ 1fb0fa62-9757-4a0c-ac67-89a6f4b68caf
md"""
If you want the exact value you can use
"""

# ╔═╡ 47ddae23-9cfd-4db1-bbe7-446c456a2fb5
md"""
We already extracted
"""

# ╔═╡ ce10d865-0edf-4f34-9142-48297fb79810
σ_residuals_exact

# ╔═╡ c1055a38-1f1a-4ef2-9fda-fe14670821b9
md"""
Following the above logic, but of course you can only do the first step:"""

# ╔═╡ 651a7d74-715a-4562-85b1-97a8bdaf27a2
σ_3_exact = VarCorr(kb07_m).σρ.subj.σ[1] / σ_residuals_exact

# ╔═╡ 3e182015-b555-4b9a-91b9-5e6a58e020b3
re_subj_exact = create_re(σ_3_exact)

# ╔═╡ 4ca9d39a-c061-4f99-b240-bd8fdd773530
md"""
As mentioned above `θ` is the compact form of these covariance matrices:
"""

# ╔═╡ 4d5abbf8-5c04-41e6-a49d-8798c985f953
kb07_m.θ = vcat( flatlowertri(re_item), flatlowertri(re_subj) )

# ╔═╡ b6d4353c-c4bf-421e-abd6-e1c4e427f1fe
md"""
We can install these parameter in the `parametricbootstrap()`-function or in the model like this:
"""

# ╔═╡ 223f7d34-e433-49f4-bae7-818ec91985df
update!(kb07_m, re_item, re_subj)

# ╔═╡ 82f87d91-963a-4fb9-96c1-162db2cd9c58
md"""
## A simple example from scratch
"""

# ╔═╡ d9eb46c8-c37a-4037-aca6-addf66690a10
md"""
Having this knowledge about the parameters we can now **simulate data from scratch**

The `simdat_crossed()` function from `MixedModelsSim` lets you set up a data frame with a specified experimental design.
For now, it only makes fully balanced crossed designs!, but you can generate an unbalanced design by simulating data for the largest cell and deleting extra rows.


Firstly we will set an easy design where `subj_n` subjects per `age` group (O or Y) respond to `item_n` items in each of two `condition`s (A or B).

Your factors need to be specified separately for between-subject, between-item, and within-subject/item factors using `Dict` with the name of each factor as the keys and vectors with the names of the levels as values.

We start with the between subject factors:
"""

# ╔═╡ 2a6a4bf8-93f6-4185-8510-bb0c0cb979c8
subj_btwn_smpl = Dict("age" => ["O", "Y"])

# ╔═╡ cd479f01-7b04-46c7-9f23-8615dee536de
md"""
There are no between-item factors in this design so you can omit it or set it to nothing.
Note that if you have factors which are between subject and between item, you need to put them in both dicts.
"""

# ╔═╡ 356b2ee9-aecc-4819-b535-eb11e46d2e52
item_btwn_smpl = nothing

# ╔═╡ ec41b5cf-461f-4655-a3c0-dc7db29e3ac6
md"""
Next, we put within-subject/item factors in a dict:
"""

# ╔═╡ 811669c1-ea76-4db2-a09a-de3480cb428b
both_win_smpl = Dict("condition" => ["A", "B"])

# ╔═╡ 20b38246-966b-4cc3-82ce-4cdb76be058e
md"""
Define subject and item number:
"""

# ╔═╡ 2a3411a4-42df-4f6c-bf38-cc7114c5ec87
begin
	subj_n_smpl = 10
	item_n_smpl = 30
end

# ╔═╡ 48f404de-f654-4021-ba98-edeaabec6e5b
md"""
Simulate data:
"""

# ╔═╡ 2164a837-b0e9-4207-90fd-4db090021963
dat_smpl = simdat_crossed(subj_n_smpl,
                     item_n_smpl,
                     subj_btwn = subj_btwn_smpl,
                     item_btwn = item_btwn_smpl,
                     both_win = both_win_smpl);

# ╔═╡ 873d708d-3f3a-4510-877d-25f73710a35b
md"Have a look:
"

# ╔═╡ d5b7f38f-2e72-44b6-9d37-131a14e47658
first(DataFrame(dat_smpl),8)

# ╔═╡ f3465ca6-8838-4f71-9f8a-df727c0f3780
first(sort(DataFrame(dat_smpl),[:subj,:item]),18)

# ╔═╡ 6155d5d8-cf14-4c73-8f41-b0ac5f40244a
md"""
The values we see in the column `dv` is just random noise.

Set contrasts:
"""

# ╔═╡ 91120f25-068e-402c-9b13-2eb8d3833b0c
contrasts_smpl = Dict(:age => HelmertCoding(),
                 :condition => HelmertCoding());

# ╔═╡ 0f627b49-72c6-4235-a7f2-bd3bd758704a
md"""Define formula:
"""

# ╔═╡ 7a9c7cd3-851b-4d2b-be12-cea014e834e5
f1 = @formula dv ~ 1 + age * condition + (1|item) + (1|subj);

# ╔═╡ b4f9eff5-3626-47ed-ab8a-8d70f926c831
md"""Note that we did not include condition as random slopes for item and subject.
This is mainly to keep the example simple and to keep the parameter `θ` easier to understand (see Section 3 for the explanation of theta).


Fit the model:"""

# ╔═╡ 5caa14d8-73b8-46d7-b0e1-bc58522f24eb
m1 = fit(MixedModel, f1, dat_smpl, contrasts=contrasts_smpl)

# ╔═╡ 97af5c8a-ab21-45c1-8fbc-39f749c4d0c7
md"""
Because the `dv` is just random noise from N(0,1), there will be basically no subject or item random variance, residual variance will be near 1.0, and the estimates for all effects should be small.
Don't worry, we'll specify fixed and random effects directly in `parametricbootstrap()`.


Set random seed for reproducibility:
"""

# ╔═╡ 7dfc2f19-bcb8-4669-a119-5238e248ac5c
rng_smpl = StableRNG(42);

# ╔═╡ c88f5177-1a24-4ffc-9228-3f8bd873eb07
md"""
Specify `β`, `σ`, and `θ`, we just made up this parameter:
"""

# ╔═╡ 72304a3d-6a3a-40a5-933f-4cc88852db11
new_beta_smpl = [0.0, 0.25, 0.25, 0.]

# ╔═╡ 99f81384-06c3-45a4-8ccf-305ab7622fe3
new_sigma_smpl = 2.0

# ╔═╡ a171e0d3-a406-42cf-9567-4d33d7be3159
new_theta_smpl = [1.0, 1.0]

# ╔═╡ e7740a25-15c2-4997-b838-22f4e8b9bf78
md"""
Run nsims iterations:
"""

# ╔═╡ 08d1781d-dad4-47c8-861f-d745d0f406bf
sim_smpl = parametricbootstrap(rng_smpl, nsims, m1,
                        β = new_beta_smpl,
                        σ = new_sigma_smpl,
                        θ = new_theta_smpl,
                        use_threads = false);

# ╔═╡ 7b7423be-dbbf-4e13-8c98-2eea4df0d02f
md"""
### Power calculation
"""

# ╔═╡ 97130e8f-f155-4c4f-a061-696bdbb3c588
ptbl_smpl= power_table(sim_smpl)

# ╔═╡ 55096ce7-b7f2-4884-a853-a87316cd171a
DataFrame(shortestcovint(sim_smpl))


# ╔═╡ 7dd1cd44-ee63-4c14-9e82-d2b7c7cef7f8
md"""
For nicely displaying it, you can use pretty_table:
"""

# ╔═╡ 63f69580-c5bc-463e-a42b-757b3d897be1
pretty_table(ptbl_smpl)

# ╔═╡ d57de97c-c21f-4e74-b093-5d3a552609e7
md"""
Convert p-values of your fixed-effects parameters to dataframe
"""

# ╔═╡ ac8ba024-a466-419d-9c80-ef9212228361
sim_smpl_df = DataFrame(sim_smpl.coefpvalues);

# ╔═╡ 0cff3191-cdd7-43dc-8761-7ad8bfc64b8c
[mean(sim_smpl_df[sim_smpl_df.coefname .== Symbol("(Intercept)"),:p] .< 0.05),
 mean(sim_smpl_df[sim_smpl_df.coefname .== Symbol("age: Y"),:p] .< 0.05),
 mean(sim_smpl_df[sim_smpl_df.coefname .== Symbol("condition: B"),:p] .< 0.05),
 mean(sim_smpl_df[sim_smpl_df.coefname .== Symbol("age: Y & condition: B"),:p] .< 0.05)]	

# ╔═╡ fc75435a-3cce-42ab-9e7d-3b9126d634b6
md"""
## Recreate a more complex dataset from scratch and analyze power for specific model parameter but various sample sizes.

### Recreate the `kb07`-dataset from scratch

For full control over all parameters in our `kb07` data set we will recreate the design using the method shown above.

Define subject and item number:
"""

# ╔═╡ a48b29ef-94da-4917-b1fe-e9101e56f319
begin
	subj_n_cmpl = 56
	item_n_cmpl = 32
end

# ╔═╡ 9dcbb3f7-c2d2-4b05-aa25-5adb4ca72d7f
md"""
Define factors in a dict:
"""

# ╔═╡ 8d13df4c-8468-482b-92f7-62f28a1bf995
begin
	subj_btwn_cmpl = nothing
	item_btwn_cmpl = nothing
	both_win_cmpl = Dict("spkr" => ["old", "new"],
	                "prec" => ["maintain", "break"],
	                "load" => ["yes", "no"]);
end

# ╔═╡ d63c51d1-1ef8-4aa9-b0f1-288c387f6d72
md"### **Try**: Play with `simdat_crossed`.
"

# ╔═╡ 7bbe248a-6ff5-44f6-9cd4-4fdbd48b81ec
begin
	subj_btwn_play = nothing
	item_btwn_play = nothing
	both_win_play  = Dict("spkr" => ["old", "new"],
					 "prec" => ["maintain", "break"],
	                 "load" => ["yes", "no"]);
	
	subj_n_play = 56
	item_n_play = 32
	
	fake_kb07_play = simdat_crossed(subj_n_play, item_n_play,
	                     subj_btwn = subj_btwn_play,
	                     item_btwn = item_btwn_play,
	                     both_win = both_win_play);
	
	
	fake_kb07_df_play = DataFrame(fake_kb07_play);
end

# ╔═╡ 462c5d9d-e579-4cc4-8f99-af7486ee9714
md"""
Simulate data:
"""

# ╔═╡ f5d79038-e98a-472b-b82c-699ee90112bf
fake_kb07 = simdat_crossed(subj_n_cmpl, item_n_cmpl,
                     subj_btwn = subj_btwn_cmpl,
                     item_btwn = item_btwn_cmpl,
                     both_win = both_win_cmpl);

# ╔═╡ 4ce55c9a-f91b-4724-a3af-349289599d45
md"""
Make a dataframe:
"""

# ╔═╡ 1346fd99-6cec-43b7-bca3-c8dfcf3bb207
fake_kb07_df = DataFrame(fake_kb07);

# ╔═╡ 18dc8d6d-7200-424a-b715-5ae14f4f178f
md"""
Have a look:
"""

# ╔═╡ 410e4de9-a672-4774-b281-562e444299a3
first(fake_kb07_df,8)

# ╔═╡ b93bcbd0-3bd1-45fc-b7f8-27319b05b290
md"""
The function `simdat_crossed` generates a balanced fully crossed design.
Unfortunately, our original design is not balanced fully crossed. Every subject saw an image only once, thus in one of eight possible conditions. To simulate that we only keep one of every eight lines.

We sort the dataframe to enable easier selection
"""

# ╔═╡ 26b5c5b4-6ae2-4bc6-83c9-1f1b0a485d5d
fake_kb07_df_sort = sort(fake_kb07_df, [:subj, :item, :load, :prec, :spkr])

# ╔═╡ 8d487d48-e80c-4ad0-8977-eaea7a6ed1e9
md"""
In order to select only the relevant rows of the data set we define an index which represents a random choice of one of every eight rows. First we generate a vector `idx` which represents which row to keep in each set of 8.
"""

# ╔═╡ d223b588-e380-458c-a546-f552f51fdda5
len = convert(Int64,(length(fake_kb07)/8));

# ╔═╡ 7d6ef8df-3230-4b90-8756-1f41b73670d4
idx_0 = rand(rng, collect(1:8) , len);

# ╔═╡ 8a04a59a-5ddd-4783-83ae-d73008673aa1
md"""
Then we create an array `A`, of the same length that is populated multiples of the number 8. Added together `A` and `idx` give the indexes of one row from each set of 8s.
"""

# ╔═╡ 0d24be79-1d59-4d44-bf08-6fa3826847ac
begin
	A_0 = repeat([8], inner=len-1)
	A_1 = append!( [0], A_0 )
	A_2 = cumsum(A_1)
	idx_1 = idx_0 + A_2
end

# ╔═╡ 7e993ae9-7a84-42f8-98e8-9195e987d942
md"""
Reduce the balanced fully crossed design to the original experimental design:
"""

# ╔═╡ f6ae17ef-3ab3-469f-970b-8aee0e3adfcd
fake_kb07_df_final= fake_kb07_df_sort[idx_1, :];

# ╔═╡ 5916d763-658b-44f6-baea-6cd11898954f
rename!(fake_kb07_df_final, :dv => :rt_trunc)

# ╔═╡ fa219b9c-52e9-45b1-9bf4-6daf1f7983c7
md"""
Now we can use the simulated data in the same way as above.

Set contrasts:
"""

# ╔═╡ cbe91c03-352c-42fe-a410-3371e023a2c8
contrasts_cmpl = Dict(:spkr => HelmertCoding(),
                 	  :prec => HelmertCoding(),
                      :load => HelmertCoding());

# ╔═╡ 35938e63-7b16-40b9-a321-fce1f3f933f7
md"""
Use formula, same as above `kb07_f`
"""

# ╔═╡ 35fe77d9-81bd-4a8e-9b5a-fe5baf28b2f6
fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df_final, contrasts=contrasts_cmpl)

# ╔═╡ a1640dce-dd08-427e-bb9d-7d3ef8495617
md"""
Use random seed for reproducibility `rng`
"""

# ╔═╡ 63220c73-f12a-4ec5-aa73-5efe9b6e2fcf
md"""
Then, again, we specify `β`, `σ`, and `θ`.
Here we use the values that we found in the model of the existing dataset:
"""

# ╔═╡ edb5e147-0fed-41d0-b48e-24537e2463f9
kb07_m

# ╔═╡ 69fc8f6b-8862-44d0-aca2-20def7028487
md"`β`"


# ╔═╡ 9ba9c512-ba39-4260-a4f3-1eebc3678afa
new_beta_cmpl = [2181.85, 67.879, -333.791, 78.5904]; #manual

# ╔═╡ c4a252b5-58df-428c-81e1-0e22353e41fc
new_beta_cmpl_exact = kb07_m.β; #grab from existing model

# ╔═╡ b350df8d-e351-43c3-870e-bae8ee3916b1
md"`σ`"

# ╔═╡ 721140a5-4986-46df-bc18-9834506e1e2e
new_sigma_cmpl = 680.032; #manual

# ╔═╡ 62dc7c9c-8fdf-42e2-a6a3-ec954e538bdc
new_sigma_cmpl_exact = kb07_m.σ; #grab from existing model

# ╔═╡ df24c706-2455-4324-8941-e66ffc084890
md"""
Have a look on original model again:
"""

# ╔═╡ adc5f8fa-cdf6-4b19-843d-a0cddd26b9ae
VarCorr(kb07_m)

# ╔═╡ 526dab6b-3409-4ffd-b1af-56eb73f5b7cb
md"`θ`"

# ╔═╡ 57efe6d8-16dc-4f60-9b78-f840f1d44fed
re_item_corr_cmpl = [1.0 -0.7; -0.7 1.0];

# ╔═╡ 71ab4138-bb3a-47a3-9d12-9a7004624dea
re_item_cmpl = create_re(0.536, 0.371; corrmat = re_item_corr_cmpl);

# ╔═╡ fd90e01d-9e14-4bac-9e1f-3f2d36cf1507
re_subj_cmpl = create_re(0.438);

# ╔═╡ fea3d964-4190-4413-aff1-5aa9ae4df0c0
new_theta_cmpl = vcat( flatlowertri(re_item_cmpl), flatlowertri(re_subj_cmpl) ) #manual

# ╔═╡ 845a4983-74ce-4cfa-84d9-ed1f043ce9f1
new_theta_cmpl_exact = kb07_m.θ; #grab from existing model

# ╔═╡ 9120bd5f-f6b3-417e-8b81-286c1c9a7a4e
md"""
Run nsims iterations:
"""

# ╔═╡ 5ab9112e-328c-44c2-bf1d-527576aedb38
fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m,
                        β = new_beta_cmpl,
                        σ = new_sigma_cmpl,
                        θ = new_theta_cmpl,
                        use_threads = false)

# ╔═╡ e9e851c1-e74f-4b8a-9a8f-c817a6c751fc
md"""
### Power calculation
"""

# ╔═╡ 5070c540-4d53-4107-8748-5e3f56fea37f
power_table(fake_kb07_sim)

# ╔═╡ 0ae98b2b-0e3f-4974-aa67-39a0882469fa
fake_kb07_sim_df = DataFrame(fake_kb07_sim.coefpvalues);

# ╔═╡ 6cb151cf-2f03-44f2-8a7a-5c8ef0236079
[mean(fake_kb07_sim_df[fake_kb07_sim_df.coefname .== Symbol("(Intercept)"),:p] .< 0.05),
mean(fake_kb07_sim_df[fake_kb07_sim_df.coefname .== Symbol("spkr: old"),:p] .< 0.05),
mean(fake_kb07_sim_df[fake_kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05),
mean(fake_kb07_sim_df[fake_kb07_sim_df.coefname .== Symbol("load: yes"),:p] .< 0.05)]

# ╔═╡ 2bc411ac-e4b8-4e85-9243-4a222935290d
md"""
Compare to the powertable from the existing data:
"""

# ╔═╡ 4ccc0b79-7f26-4994-8126-88955460e43a
[mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("(Intercept)"),:p] .< 0.05),
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("spkr: old"),:p] .< 0.05),
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05),
mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("load: yes"),:p] .< 0.05)]

# ╔═╡ ea53f15b-db89-444c-bdde-25747751d505
md"""
We have successfully recreated the power simulation of an existing dataset from scratch. This has the advantage, that we now can iterate over different numbers of subjects and items.
"""

# ╔═╡ de802935-40fe-4010-b361-5ce1d1a6e19e
md"""
Compare to the powertable from the existing data:
"""

# ╔═╡ f273f990-f5b3-4f7e-a066-4da7bf71c165
power_table(kb07_sim)

# ╔═╡ f4b1e8af-7a3a-43ae-a189-080fa9886c44
md"""
We have successfully recreated the power simulation of an existing dataset from scratch. This has the advantage, that we now can iterate over different numbers of subjects and items.

## Loop over subject and item sizes

When designing a study, you may be interested in trying various numbers of subjects and items to see how that affects the power of your study.
To do this you can use a loop to run simulations over a range of values for any parameter.

### We first define every fixed things outside the loop (same as above):

Define factors in a dict:

"""

# ╔═╡ eca8b9d0-3c9b-4be4-bb2b-bbbaed248bfc
md"""
`subj_btwn_cmpl`, `item_btwn_cmpl`, `both_win_cmpl`
"""

# ╔═╡ b53158c1-80f6-4286-ae89-6273e47d4b95
md"""
`contrasts_cmpl`
"""

# ╔═╡ a0df1382-68d2-4b29-be29-3b3ebb94cb1e
md"`rng`"

# ╔═╡ 98244829-7c17-42f0-a10b-1d86408a4fb7
md"""
`kb07_f`
"""

# ╔═╡ 02d09aac-70ac-48ed-829a-3be9b2fa605c
md"""
`new_beta_cmpl`, `new_sigma_cmpl`, `new_theta_cmpl`
"""

# ╔═╡ 806491f8-e739-4537-8768-125d9124d157
md"""
### Then we define the variables that out loop will iterate over

Define subject and item numbers as arrays:
"""

# ╔═╡ 1f026d1a-9d5c-4732-9a42-bb85be9cff0b
sub_ns = [20, 30, 40];

# ╔═╡ 1ec45e86-fce4-43a6-b6dd-74972b316113
item_ns = [16, 24, 32];

# ╔═╡ 34b70baa-46f2-42cc-98bb-06b6d2c38ee3
md"""
Make an empty dataframe:
"""

# ╔═╡ 77aa8869-176d-4bf4-a021-b588b769e44f
d = DataFrame();

# ╔═╡ 2f8f5fa3-892a-4641-84e9-b0dfa59d2b27
md"""
### Run the loop:
"""

# ╔═╡ 789c6ecd-9933-49a5-b3fa-5c3ff86488b3
for subj_n in sub_ns
    for item_n in item_ns

    # Make balanced fully crossed data:
    fake_kb07 = simdat_crossed(subj_n, item_n,
                     subj_btwn = subj_btwn,
                     item_btwn = item_btwn,
                     both_win = both_win);
    fake_kb07_df = DataFrame(fake_kb07)

    # Reduce the balanced fully crossed design to the original experimental design:
    fake_kb07_df = sort(fake_kb07_df, [:subj, :item, :load, :prec, :spkr])
    len = convert(Int64,(length(fake_kb07)/8))
    idx = rand(rng, collect(1:8) , len)
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

# ╔═╡ f117902c-7a84-4a6c-ab6a-0dd11cddee16
md"""
Our dataframe `d` now contains the power information for each combination of subjects and items.
"""

# ╔═╡ dbc5569b-cb73-44c0-970c-16f713349163
print(d)

# ╔═╡ a562900c-93af-4562-8f92-445e59bc9481
md"""
Lastly we plot our results:
"""

# ╔═╡ 9165b8ac-4962-432f-9643-525f441d9615
begin
	categorical!(d, :item_n)
	
	plot(d, x="subj_n", y="power",xgroup= "coefname",color="item_n", Geom.subplot_grid(Geom.point, Geom.line), Guide.xlabel("Number of subjects by parameter"), Guide.ylabel("Power"))
end

# ╔═╡ e54b4939-f389-49f0-90c0-1a0b2ff87774
md"""
# Credit
This tutorial was conceived for ZiF research and tutorial workshop by Lisa DeBruine (Feb. 2020) presented again by Phillip Alday during the SMLP Summer School (Sep. 2020).

Updated and extended by Lisa Schwetlick & Daniel Backhaus, with the kind help of Phillip Alday, after changes to the package.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
MixedModels = "ff71e718-51f3-5ec2-a782-8ffcbfa3c316"
MixedModelsSim = "d5ae56c5-23ca-4a1f-b505-9fc4796fc1fe"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StableRNGs = "860ef19b-820b-49d6-a774-d7a799459cd3"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[compat]
CSV = "~0.9.1"
DataFrames = "~1.2.2"
DataFramesMeta = "~0.9.1"
Gadfly = "~1.3.3"
MixedModels = "~4.1.1"
MixedModelsSim = "~0.2.3"
PlutoUI = "~0.7.9"
StableRNGs = "~1.0.0"
Tables = "~1.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs"]
git-tree-sha1 = "b00e6eaba895683867728e73af78a00218f0db10"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "1.6.2"

[[ArrowTypes]]
deps = ["UUIDs"]
git-tree-sha1 = "a0633b6d6efabf3f76dacd6eb1b3ec6c42ab0552"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "1.2.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "42ac5e523869a84eac9669eaceed9e4aa0e1587b"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.4"

[[BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "f50b5a99aa6ff9db7bf51255b5c21c8bc871ad54"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.2.5"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "c907e91e253751f5840135f4c9deb1308273338d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.9.1"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "RecipesBase", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "1562002780515d2573a4fb0c3715e4e57481075e"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.0"

[[Chain]]
git-tree-sha1 = "cac464e71767e8a04ceee82a889ca56502795705"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.8"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "30ee06de5ff870b45c78f529a6b093b3323256a3"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.1"

[[CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "59fe0cb37784288d6b9f1baebddbf75457395d40"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.0"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "d19cd9ae79ef31774151637492291d75194fc5fa"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.7.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "6071cb87be6a444ac75fdbf51b8e7273808ce62f"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.35.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoupledFields]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "6c9671364c68c1158ac2524ac881536195b7e7bc"
uuid = "7ad07ef1-bdf2-5661-9d2b-286fd4296dac"
version = "0.2.0"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "bec2532f8adb82005476c141ec23e921fc20971b"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.8.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "Reexport"]
git-tree-sha1 = "29e71b438935977f8905c0cb3a8a84475fc70101"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.9.1"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "9f46deb4d4ee4494ffb5a40a27a2aced67bdd838"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.4"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a837fdf80f333415b69684ba8e8ae6ba76de6aaa"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.24.18"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FilePathsBase]]
deps = ["Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "0f5e8d0cb91a6386ba47bd1527b240bd5725fbae"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.10"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "693210145367e7685d8604aee33d9bfb85db8b31"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.11.9"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "f564ce4af5e79bb88ff1f4488e64363487674278"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.5.1"

[[Gadfly]]
deps = ["Base64", "CategoricalArrays", "Colors", "Compose", "Contour", "CoupledFields", "DataAPI", "DataStructures", "Dates", "Distributions", "DocStringExtensions", "Hexagons", "IndirectArrays", "IterTools", "JSON", "Juno", "KernelDensity", "LinearAlgebra", "Loess", "Measures", "Printf", "REPL", "Random", "Requires", "Showoff", "Statistics"]
git-tree-sha1 = "96da4818e4d481a29aa7d66aac1eb778432fb89a"
uuid = "c91e804a-d5a3-530f-b6f0-dfbca275c004"
version = "1.3.3"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[Hexagons]]
deps = ["Test"]
git-tree-sha1 = "de4a6f9e7c4710ced6838ca906f81905f7385fd6"
uuid = "a1b4810d-1bce-5fbd-ac56-80944d57a21f"
version = "0.2.0"

[[IndirectArrays]]
git-tree-sha1 = "c2a145a145dc03a7620af1444e0264ef907bd44f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JSON3]]
deps = ["Dates", "Mmap", "Parsers", "StructTypes", "UUIDs"]
git-tree-sha1 = "b3e5984da3c6c95bcf6931760387ff2e64f508f3"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.9.1"

[[Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "b5254a86cf65944c68ed938e575f5c81d5dfe4cb"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.3"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "86197a8ecb06e222d66797b0c2d2f0cc7b69e42b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "c253236b0ed414624b083e6b72bfe891fbd2c7af"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+1"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "a1f9933fa00624d8c97301253d14710b80fa08ee"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.10.1"

[[MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

[[MixedModels]]
deps = ["Arrow", "DataAPI", "Distributions", "GLM", "JSON3", "LazyArtifacts", "LinearAlgebra", "Markdown", "NLopt", "PooledArrays", "ProgressMeter", "Random", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "StatsModels", "StructTypes", "Tables"]
git-tree-sha1 = "f318e42a48ec0a856292bafeec6b07aed3f6d600"
uuid = "ff71e718-51f3-5ec2-a782-8ffcbfa3c316"
version = "4.1.1"

[[MixedModelsSim]]
deps = ["LinearAlgebra", "MixedModels", "PooledArrays", "PrettyTables", "Random", "Statistics", "Tables"]
git-tree-sha1 = "ad4eaa164a5ab5fd22effcb86e7c991192ed3488"
uuid = "d5ae56c5-23ca-4a1f-b505-9fc4796fc1fe"
version = "0.2.3"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NLopt]]
deps = ["MathOptInterface", "MathProgBase", "NLopt_jll"]
git-tree-sha1 = "f115030b9325ca09ef1619ba0617b2a64101ce84"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "0.6.4"

[[NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2b597c46900f5f811bec31f0dcc88b45744a2a09"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.7.0+0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c870a0d713b51e4b49be6432eff0e26a4325afee"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.6"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "7dff99fbc740e2f8228c6878e2aad6d7c2678098"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.1"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "54f37736d8934a12a200edea2f9206b03bdf3159"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.7"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[StableRNGs]]
deps = ["Random", "Test"]
git-tree-sha1 = "3be7d49667040add7ee151fefaf1f8c04c8c8276"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "46d7ccc7104860c38b11966dd1f72ff042f382e4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.10"

[[StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "3fa15c1f8be168e76d59097f66970adc86bfeb95"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.25"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "8445bf99a36d703a09c601f9a57e2f83000ef2ae"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "368d04a820fe069f9080ff1b432147a6203c3c89"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[WeakRefStrings]]
deps = ["DataAPI", "Parsers"]
git-tree-sha1 = "4a4cfb1ae5f26202db4f0320ac9344b3372136b0"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.3.0"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─c071c03f-6ddb-415e-8445-b4fd8d45abc1
# ╠═89495294-10ed-11ec-04a9-db58dba9f3d1
# ╠═14ae0269-ea23-4c73-a04f-40eece7a5bc5
# ╟─76df9a69-7ee6-4876-9dcc-e6178f59d3ad
# ╠═16822577-6f44-4330-9181-678bb2c7cce2
# ╟─d4e9bb7e-3f70-48fb-a73c-9606cfb39f39
# ╠═4c00bf63-1d33-4f16-9776-e2f3915850f9
# ╟─87af988a-2d55-4c24-a132-3c76f1ae562b
# ╠═8caad72c-a8d7-4e23-8eb6-b9096141b471
# ╟─29b32ea9-0018-4ec3-9095-1e527e11fc4f
# ╠═3763460c-046d-46ce-b9ef-4faf1fbd0428
# ╟─f05a1f73-503a-4ed7-8879-1da242b2b224
# ╠═a75ad071-0789-41ab-9253-27c8bca615ad
# ╟─8e4629aa-58b7-4ab8-8954-74ee579a3f9b
# ╠═3e0e9908-002c-4def-a2eb-c8f7d775f205
# ╟─33d9ecae-ac5c-4fdf-bd1a-4d68f785084a
# ╠═b3e2dea0-bfb0-4d5d-bd9c-8ed0e0272cc0
# ╟─1b814e14-136e-45ad-83e9-c783e30dc4b1
# ╠═fc141d66-3252-4989-b6e4-98bbfe82385a
# ╠═5d7c1822-453c-4483-a1a1-32e149145afb
# ╠═ef3c9d81-38cf-45d6-9258-926ad209e36f
# ╟─146c6a52-a512-4601-9b58-bc56575bfb2e
# ╠═4598ae9a-8d39-4280-a74a-4d7145d64920
# ╟─fa91a0f8-675e-4155-9ff4-58090689c397
# ╠═811a81ac-50d9-4d07-98d1-3fc8e1586685
# ╟─d8257d21-8dad-4533-ba93-78ba7c3202b3
# ╠═ac2f3b62-e86c-43b2-8109-558db79e5773
# ╟─67493854-920e-477e-805b-21d58cd19ade
# ╠═8c06724f-a947-44c2-8541-7aa50c915356
# ╟─fc239c6e-ca35-4df6-81b2-115be160437b
# ╠═97907ccc-e4b8-43a3-9e5a-ae1da59a203a
# ╠═a383994a-e4f4-4c9e-9f17-32ec11e32e87
# ╠═523c1c21-43d2-469e-b138-3dd5cab1e6b9
# ╠═3419a020-7333-4cb5-981b-84e7fa788280
# ╠═231ce20d-ea0c-4314-a7b0-0e296bc2160b
# ╟─b66ab7be-0547-4535-9d33-cfae62441c61
# ╠═356751b1-73cc-4dae-ad5a-b1cd0b45ede0
# ╟─910ca4be-d4bf-4231-95d6-de8b3e97d4d5
# ╟─11eacaa9-f309-4c72-a659-3fc262f3111f
# ╠═bc6ace0a-2776-409c-892c-6892df79f0e7
# ╠═bbe93bad-4d41-4300-910a-ec16a660def3
# ╠═9510ce0a-6ba7-4903-8e2a-84f6fc69492c
# ╟─f5436bb7-9087-437f-bc43-1b113c907695
# ╠═9fbfff0f-289c-4c92-9ffa-615ad1b69ff2
# ╟─7aa2cea2-9136-42ad-bc62-9966aa93b84c
# ╠═f7a5cffe-2819-4a9f-a56e-45fefb2c4201
# ╠═c3c9ab9c-c81d-44d2-bd45-fb72a7c5da8a
# ╠═9e9c337e-556f-44a8-ad5f-0476be1458b0
# ╠═6b7fdacd-1343-42a2-8e7e-d06bd559765d
# ╠═f232988a-c3cf-404f-9f0d-fe0ef591e8bb
# ╠═ed309fad-0058-4069-8e8e-e526b2b36370
# ╠═2fc97a13-abc5-42ff-a1eb-28e24271fae0
# ╟─42f4e755-c7b1-46e2-91fd-31a19ccd5685
# ╟─fcd25a28-cb0d-4077-8b20-f5837c4e99f8
# ╟─ae7efbb0-aa35-4128-aebf-10a1deab5fee
# ╠═21567edd-0593-4823-b538-15dcb0de48f0
# ╠═5eb6074a-a695-4404-841e-21b7b82f1150
# ╟─a02e6191-239f-4b0b-996d-f73bfc54ea28
# ╠═ce275097-63e5-410e-a23b-9fdd4250e35c
# ╠═c5f9f6ad-0b18-434b-a368-42afd8cb398e
# ╟─9adf46f1-d84c-4ae3-af1b-2087adbf0315
# ╠═077b22ae-814f-4b18-8297-2124317cbff6
# ╟─f33d0b96-3b3c-4651-8b39-46461445cd09
# ╠═8103cf6e-18a2-40a1-8043-c87ed92e1c11
# ╟─cdee3c0c-6c85-4515-8ffd-1b1c71018e09
# ╠═f5535f46-e919-4dad-9f0e-337d68ab669b
# ╟─6f2c0659-3f5d-45d8-8a40-b3124ce49fa3
# ╠═9b5601c0-613e-426f-9318-c87ba9108f42
# ╟─d0493428-916e-4bfe-8f5a-586e21af0f5c
# ╠═ca0b404c-50da-43eb-8bb6-d6037e09ba64
# ╟─2a04994e-bb88-455b-a732-bae9ef1c5422
# ╟─5e364a80-e5d9-458b-be1e-d00dafa04033
# ╠═026c051d-bcfa-4b1d-8648-bdd0d8cdf894
# ╠═2951d7fe-0373-4f52-86e3-c11bd7a57f1a
# ╠═1ddbde42-6dcc-425d-b906-92880bb15e89
# ╠═29ef3cab-e774-4cf8-a392-92eec8acc98b
# ╠═37f07a78-3616-41c4-90ad-6fc52d9223b3
# ╠═f4929a81-71db-4a78-86f9-22c7ed89f04d
# ╟─dd8aa406-ef50-4036-9d6f-f499934a84c3
# ╠═f8432a72-da18-4540-9a26-2c4091bd43ac
# ╟─1fb0fa62-9757-4a0c-ac67-89a6f4b68caf
# ╟─47ddae23-9cfd-4db1-bbe7-446c456a2fb5
# ╠═ce10d865-0edf-4f34-9142-48297fb79810
# ╟─c1055a38-1f1a-4ef2-9fda-fe14670821b9
# ╠═651a7d74-715a-4562-85b1-97a8bdaf27a2
# ╠═3e182015-b555-4b9a-91b9-5e6a58e020b3
# ╟─4ca9d39a-c061-4f99-b240-bd8fdd773530
# ╠═4d5abbf8-5c04-41e6-a49d-8798c985f953
# ╟─b6d4353c-c4bf-421e-abd6-e1c4e427f1fe
# ╠═223f7d34-e433-49f4-bae7-818ec91985df
# ╟─82f87d91-963a-4fb9-96c1-162db2cd9c58
# ╟─d9eb46c8-c37a-4037-aca6-addf66690a10
# ╠═2a6a4bf8-93f6-4185-8510-bb0c0cb979c8
# ╟─cd479f01-7b04-46c7-9f23-8615dee536de
# ╠═356b2ee9-aecc-4819-b535-eb11e46d2e52
# ╟─ec41b5cf-461f-4655-a3c0-dc7db29e3ac6
# ╠═811669c1-ea76-4db2-a09a-de3480cb428b
# ╟─20b38246-966b-4cc3-82ce-4cdb76be058e
# ╠═2a3411a4-42df-4f6c-bf38-cc7114c5ec87
# ╟─48f404de-f654-4021-ba98-edeaabec6e5b
# ╠═2164a837-b0e9-4207-90fd-4db090021963
# ╟─873d708d-3f3a-4510-877d-25f73710a35b
# ╠═d5b7f38f-2e72-44b6-9d37-131a14e47658
# ╠═f3465ca6-8838-4f71-9f8a-df727c0f3780
# ╟─6155d5d8-cf14-4c73-8f41-b0ac5f40244a
# ╠═91120f25-068e-402c-9b13-2eb8d3833b0c
# ╟─0f627b49-72c6-4235-a7f2-bd3bd758704a
# ╠═7a9c7cd3-851b-4d2b-be12-cea014e834e5
# ╟─b4f9eff5-3626-47ed-ab8a-8d70f926c831
# ╠═5caa14d8-73b8-46d7-b0e1-bc58522f24eb
# ╟─97af5c8a-ab21-45c1-8fbc-39f749c4d0c7
# ╠═7dfc2f19-bcb8-4669-a119-5238e248ac5c
# ╟─c88f5177-1a24-4ffc-9228-3f8bd873eb07
# ╠═72304a3d-6a3a-40a5-933f-4cc88852db11
# ╠═99f81384-06c3-45a4-8ccf-305ab7622fe3
# ╠═a171e0d3-a406-42cf-9567-4d33d7be3159
# ╟─e7740a25-15c2-4997-b838-22f4e8b9bf78
# ╠═08d1781d-dad4-47c8-861f-d745d0f406bf
# ╟─7b7423be-dbbf-4e13-8c98-2eea4df0d02f
# ╠═97130e8f-f155-4c4f-a061-696bdbb3c588
# ╠═55096ce7-b7f2-4884-a853-a87316cd171a
# ╟─7dd1cd44-ee63-4c14-9e82-d2b7c7cef7f8
# ╠═63f69580-c5bc-463e-a42b-757b3d897be1
# ╟─d57de97c-c21f-4e74-b093-5d3a552609e7
# ╠═ac8ba024-a466-419d-9c80-ef9212228361
# ╠═0cff3191-cdd7-43dc-8761-7ad8bfc64b8c
# ╟─fc75435a-3cce-42ab-9e7d-3b9126d634b6
# ╠═a48b29ef-94da-4917-b1fe-e9101e56f319
# ╟─9dcbb3f7-c2d2-4b05-aa25-5adb4ca72d7f
# ╠═8d13df4c-8468-482b-92f7-62f28a1bf995
# ╟─d63c51d1-1ef8-4aa9-b0f1-288c387f6d72
# ╠═7bbe248a-6ff5-44f6-9cd4-4fdbd48b81ec
# ╟─462c5d9d-e579-4cc4-8f99-af7486ee9714
# ╠═f5d79038-e98a-472b-b82c-699ee90112bf
# ╟─4ce55c9a-f91b-4724-a3af-349289599d45
# ╠═1346fd99-6cec-43b7-bca3-c8dfcf3bb207
# ╟─18dc8d6d-7200-424a-b715-5ae14f4f178f
# ╠═410e4de9-a672-4774-b281-562e444299a3
# ╟─b93bcbd0-3bd1-45fc-b7f8-27319b05b290
# ╠═26b5c5b4-6ae2-4bc6-83c9-1f1b0a485d5d
# ╟─8d487d48-e80c-4ad0-8977-eaea7a6ed1e9
# ╠═d223b588-e380-458c-a546-f552f51fdda5
# ╠═7d6ef8df-3230-4b90-8756-1f41b73670d4
# ╟─8a04a59a-5ddd-4783-83ae-d73008673aa1
# ╠═0d24be79-1d59-4d44-bf08-6fa3826847ac
# ╟─7e993ae9-7a84-42f8-98e8-9195e987d942
# ╠═f6ae17ef-3ab3-469f-970b-8aee0e3adfcd
# ╠═5916d763-658b-44f6-baea-6cd11898954f
# ╟─fa219b9c-52e9-45b1-9bf4-6daf1f7983c7
# ╠═cbe91c03-352c-42fe-a410-3371e023a2c8
# ╟─35938e63-7b16-40b9-a321-fce1f3f933f7
# ╠═35fe77d9-81bd-4a8e-9b5a-fe5baf28b2f6
# ╟─a1640dce-dd08-427e-bb9d-7d3ef8495617
# ╟─63220c73-f12a-4ec5-aa73-5efe9b6e2fcf
# ╠═edb5e147-0fed-41d0-b48e-24537e2463f9
# ╟─69fc8f6b-8862-44d0-aca2-20def7028487
# ╠═9ba9c512-ba39-4260-a4f3-1eebc3678afa
# ╠═c4a252b5-58df-428c-81e1-0e22353e41fc
# ╟─b350df8d-e351-43c3-870e-bae8ee3916b1
# ╠═721140a5-4986-46df-bc18-9834506e1e2e
# ╠═62dc7c9c-8fdf-42e2-a6a3-ec954e538bdc
# ╟─df24c706-2455-4324-8941-e66ffc084890
# ╠═adc5f8fa-cdf6-4b19-843d-a0cddd26b9ae
# ╟─526dab6b-3409-4ffd-b1af-56eb73f5b7cb
# ╠═57efe6d8-16dc-4f60-9b78-f840f1d44fed
# ╠═71ab4138-bb3a-47a3-9d12-9a7004624dea
# ╠═fd90e01d-9e14-4bac-9e1f-3f2d36cf1507
# ╠═fea3d964-4190-4413-aff1-5aa9ae4df0c0
# ╠═845a4983-74ce-4cfa-84d9-ed1f043ce9f1
# ╟─9120bd5f-f6b3-417e-8b81-286c1c9a7a4e
# ╠═5ab9112e-328c-44c2-bf1d-527576aedb38
# ╟─e9e851c1-e74f-4b8a-9a8f-c817a6c751fc
# ╠═5070c540-4d53-4107-8748-5e3f56fea37f
# ╠═0ae98b2b-0e3f-4974-aa67-39a0882469fa
# ╠═6cb151cf-2f03-44f2-8a7a-5c8ef0236079
# ╟─2bc411ac-e4b8-4e85-9243-4a222935290d
# ╠═4ccc0b79-7f26-4994-8126-88955460e43a
# ╠═ea53f15b-db89-444c-bdde-25747751d505
# ╟─de802935-40fe-4010-b361-5ce1d1a6e19e
# ╠═f273f990-f5b3-4f7e-a066-4da7bf71c165
# ╟─f4b1e8af-7a3a-43ae-a189-080fa9886c44
# ╟─eca8b9d0-3c9b-4be4-bb2b-bbbaed248bfc
# ╟─b53158c1-80f6-4286-ae89-6273e47d4b95
# ╟─a0df1382-68d2-4b29-be29-3b3ebb94cb1e
# ╟─98244829-7c17-42f0-a10b-1d86408a4fb7
# ╟─02d09aac-70ac-48ed-829a-3be9b2fa605c
# ╟─806491f8-e739-4537-8768-125d9124d157
# ╠═1f026d1a-9d5c-4732-9a42-bb85be9cff0b
# ╠═1ec45e86-fce4-43a6-b6dd-74972b316113
# ╟─34b70baa-46f2-42cc-98bb-06b6d2c38ee3
# ╠═77aa8869-176d-4bf4-a021-b588b769e44f
# ╟─2f8f5fa3-892a-4641-84e9-b0dfa59d2b27
# ╠═789c6ecd-9933-49a5-b3fa-5c3ff86488b3
# ╟─f117902c-7a84-4a6c-ab6a-0dd11cddee16
# ╠═dbc5569b-cb73-44c0-970c-16f713349163
# ╠═a562900c-93af-4562-8f92-445e59bc9481
# ╠═9165b8ac-4962-432f-9643-525f441d9615
# ╟─e54b4939-f389-49f0-90c0-1a0b2ff87774
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
