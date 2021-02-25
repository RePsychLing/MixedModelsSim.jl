### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 40bf1d04-76d0-11eb-0d98-cf2945aa3fd4
using DataFrames, MixedModels, MixedModelsSim, Random

# ╔═╡ 09276d18-76d2-11eb-3b91-b77f2d532e82
md"""
# Simulation of linear mixed models

This notebook shows the process of simulating responses from a mixed-effects model with crossed random effects for subject (`subj`) and item (`item`).
The general approach is to create a data frame with the appropriate experimental and grouping factors and with an _inert_ dependent variable (`dv`), meaning that, at this stage `dv` is just random noise.

This model is then used in a parametric bootstrap simulation with non-inert coefficient values (β) and random-effects variances and covariances.

Begin by loading the packages to be used.
"""

# ╔═╡ 28a121e2-76d3-11eb-0896-43da695be0a4
md"""
And initialize a random number generator.
"""

# ╔═╡ b64a889c-76d0-11eb-1936-5545db219def
rng = MersenneTwister(8675309);

# ╔═╡ ce3be2ec-76e0-11eb-3fb5-67bd31816bca
md"""
Next we set the number of subjects, `sub_n`, and the number of items, `item_n`
"""

# ╔═╡ 21958b24-76ed-11eb-11a2-41223c4864f8
begin
	sub_n = 200
	item_n = 50
end;

# ╔═╡ b8fe4dd2-76e2-11eb-3a7b-21e2b17043f7
md"A between-subject factor, `cond`, with two levels, `easy` and `hard` is specified as"

# ╔═╡ feb89de4-76d0-11eb-115e-89e7ed2fc483
subj_bt = Dict(:cond => ["easy", "hard"]);

# ╔═╡ eb911720-76e2-11eb-002d-21a35fd7d4b2
md"These settings are used to create a Table with an _inert_ dependent variable, `dv`, as"

# ╔═╡ 286561a4-76d1-11eb-2211-61d84bc1b9a9
frm = simdat_crossed(rng, sub_n, item_n, subj_btwn = subj_bt);

# ╔═╡ aa915fd2-76ec-11eb-2cec-bff137f10b18
md"This data table is in the form of a _row table_, which is a vector of `NamedTuple`s"

# ╔═╡ cb5f962a-76ec-11eb-3b1e-6f7f1f68bca5
typeof(frm)

# ╔═╡ 6914997e-76e3-11eb-2f2b-0501a82bee66
md"It is generally best to convert this table to a `DataFrame` for display or summary"

# ╔═╡ 9f327062-76e3-11eb-16dc-7be3a2b1c09f
df = DataFrame(frm)

# ╔═╡ b0202450-76e3-11eb-3218-8bddeed1c50d
describe(df)

# ╔═╡ 0fce98d0-76e4-11eb-0ed4-430ae7972fce
md"Note that if you go back to the entry field of `sub_n` or `item_n` and change the value, the other results are automatically updated"

# ╔═╡ 34b8677c-76e4-11eb-1807-9306323a7c9f
md"""
## Fitting a model to the inert data

In the `MixedModels` package the formula language is similar to that in the `lme4` package for `R` with the exception that the formula must be wrapped in a call to the macro `@formula`
"""

# ╔═╡ 48e92618-76d1-11eb-053c-e77449325b9c
m0form = @formula dv ~ 1 + cond + (1|subj) + (1|item);

# ╔═╡ 7292d5fa-76e4-11eb-243f-8f7ad404d683
md"Also, contrasts are specified in the call to fit the model instead of being part of the metadata for the factor"

# ╔═╡ 761fb192-76d1-11eb-2f8e-890af7ae08d8
m0 = fit(MixedModel, m0form, frm, contrasts=Dict(:cond => EffectsCoding()))

# ╔═╡ 91169880-76d1-11eb-28da-0f4135f0e995
md"""
As this dependent variable is inert (i.e. simulated from a standard Normal distribution without any fixed- or random-effects contributions) it is not surprising the the parameter estimates for the variance components for subject and item are close to zero as are the estimates for the fixed-effects.

## Simulated model fits for non-inert data

To simulate responses from non-intert data we must specify values for the parameters in the form of a scalar, `σ`, the standard deviation of the per-observation noise, and two vectors: `β`, the fixed-effects parameter vector and `θ`, the vector of relative standard deviations.

The fixed-effects parameter is relative to the model matrix, `X`, which, in this case, is the overall mean DV for a typical condition and half the difference between the `hard` and `easy` conditions.  It is half the difference because we use a +/- 1 coding for the `cond` factor, producing a model matrix, `X`, of
"""

# ╔═╡ d217cc16-76e8-11eb-2f64-6ff4931c87e9
Int.(m0.X')   # convert to integers to save space in the display

# ╔═╡ 9d46c892-76e9-11eb-008e-67d5100eb3fa
md"Suppose we set the population mean response time to 400 ms. and the difference between `hard` and `easy` to be 50 ms.
The β vector is then
"

# ╔═╡ dbd15424-76e9-11eb-2306-436929c4e756
β = [400.0, 25.0]   # Include the decimal point to force Float64 values

# ╔═╡ fbdfaf02-76e9-11eb-2404-d9a288cf1d31
md"We also set"

# ╔═╡ 167234ae-76ea-11eb-2589-93a3e6940453
begin
	σ  =  25.0  # s.d. of per-observation noise
	σ₁ = 100.0  # s.d. of random effects for subject
	σ₂ =  50.0  # s.d. of random effects for item
	θ = [σ₁, σ₂] ./ σ  # relative standard deviations
end;

# ╔═╡ 5d1b6440-7783-11eb-3bae-45d008439b01
md"""
To simulate a single response vector from `m0` with these parameters and re-estimate the parameters we use
"""

# ╔═╡ ac43c242-7783-11eb-17d4-19939f282c52
simulate!(MersenneTwister(12321), m0; β=β, σ=σ, θ=θ) |> refit!

# ╔═╡ fabb4eae-7783-11eb-13de-157f033b8c12
md"""
We see that the estimates of the standard deviations of the random effects and the residual standard deviation are very close to the values used in the simulation: 100.0, 50.0, and 25.0.

The fixed-effects coefficients are also close to those used in the simulation, [400.0, 25.0].
The standard errors of these estimates indicate that this design and data will have considerable power in testing whether the main effect for `cond` could be zero.

In general we are not interested in the power of a test of whether the `(Intercept)` could be zero.
We do not expect response times of zero or, even more peculiar, negative response times.
"""

# ╔═╡ b75c8580-76eb-11eb-2cca-37a0702e44a5
md"""
## Simulation of many cases

We could repeat this process of `simulate!` and `refit!` to create a collection of parameter estimates and derived statistics but it is easier to use the `parametricbootstrap` function to do this simulation and package the results.

We start the random number generator at the same seed as before so we can check that the individual parameter estimates are as expected.

We will simulate 2000 cases.
"""

# ╔═╡ cd924f92-76eb-11eb-3275-eb4f55282d5d
sim = parametricbootstrap(MersenneTwister(12321), 2000, m0; β=β, σ=σ, θ=θ);

# ╔═╡ 91c8aa1a-777c-11eb-0989-45fd3982ce8c
md"""
`sim` contains results of the simulation in a compact form with several "properties" defined that allow for extraction of quantities of interest.
"""

# ╔═╡ 0a8b0560-777d-11eb-2a90-25e0581a9407
typeof(sim)

# ╔═╡ 10ce3212-777d-11eb-13dd-2b02eaa97ef7
propertynames(sim)

# ╔═╡ 22151cfa-777d-11eb-2769-75cc409a9902
md"""
If, for example, we wanted to plot the p-values for the fixed-effects coefficients, we extract them as
"""

# ╔═╡ 598ccad6-777d-11eb-3cac-3fcd0f0f3969
pvals = sim.coefpvalues;

# ╔═╡ 6a6a3e0e-777d-11eb-09f7-99de6d4a20be
typeof(pvals)

# ╔═╡ 6fe7f0c6-777d-11eb-0d78-dd28ca754de3
length(pvals)

# ╔═╡ 802d53a4-777d-11eb-25cc-efc4d9dcf292
pvalscond = last(groupby(select(DataFrame(pvals), :iter, :coefname, :p), :coefname))

# ╔═╡ 49df45ee-7782-11eb-2363-d562aa91ca14
extrema(pvalscond.p)

# ╔═╡ 3f7a0d3e-76ec-11eb-368e-17428098bc72
ptbl = power_table(sim);

# ╔═╡ 52e23b14-76ec-11eb-1787-01b2cf9888eb
DataFrame(ptbl)

# ╔═╡ Cell order:
# ╟─09276d18-76d2-11eb-3b91-b77f2d532e82
# ╠═40bf1d04-76d0-11eb-0d98-cf2945aa3fd4
# ╟─28a121e2-76d3-11eb-0896-43da695be0a4
# ╠═b64a889c-76d0-11eb-1936-5545db219def
# ╟─ce3be2ec-76e0-11eb-3fb5-67bd31816bca
# ╠═21958b24-76ed-11eb-11a2-41223c4864f8
# ╟─b8fe4dd2-76e2-11eb-3a7b-21e2b17043f7
# ╠═feb89de4-76d0-11eb-115e-89e7ed2fc483
# ╟─eb911720-76e2-11eb-002d-21a35fd7d4b2
# ╠═286561a4-76d1-11eb-2211-61d84bc1b9a9
# ╟─aa915fd2-76ec-11eb-2cec-bff137f10b18
# ╠═cb5f962a-76ec-11eb-3b1e-6f7f1f68bca5
# ╟─6914997e-76e3-11eb-2f2b-0501a82bee66
# ╠═9f327062-76e3-11eb-16dc-7be3a2b1c09f
# ╠═b0202450-76e3-11eb-3218-8bddeed1c50d
# ╟─0fce98d0-76e4-11eb-0ed4-430ae7972fce
# ╟─34b8677c-76e4-11eb-1807-9306323a7c9f
# ╠═48e92618-76d1-11eb-053c-e77449325b9c
# ╟─7292d5fa-76e4-11eb-243f-8f7ad404d683
# ╠═761fb192-76d1-11eb-2f8e-890af7ae08d8
# ╟─91169880-76d1-11eb-28da-0f4135f0e995
# ╠═d217cc16-76e8-11eb-2f64-6ff4931c87e9
# ╟─9d46c892-76e9-11eb-008e-67d5100eb3fa
# ╠═dbd15424-76e9-11eb-2306-436929c4e756
# ╟─fbdfaf02-76e9-11eb-2404-d9a288cf1d31
# ╠═167234ae-76ea-11eb-2589-93a3e6940453
# ╟─5d1b6440-7783-11eb-3bae-45d008439b01
# ╠═ac43c242-7783-11eb-17d4-19939f282c52
# ╟─fabb4eae-7783-11eb-13de-157f033b8c12
# ╟─b75c8580-76eb-11eb-2cca-37a0702e44a5
# ╠═cd924f92-76eb-11eb-3275-eb4f55282d5d
# ╟─91c8aa1a-777c-11eb-0989-45fd3982ce8c
# ╠═0a8b0560-777d-11eb-2a90-25e0581a9407
# ╠═10ce3212-777d-11eb-13dd-2b02eaa97ef7
# ╟─22151cfa-777d-11eb-2769-75cc409a9902
# ╠═598ccad6-777d-11eb-3cac-3fcd0f0f3969
# ╠═6a6a3e0e-777d-11eb-09f7-99de6d4a20be
# ╠═6fe7f0c6-777d-11eb-0d78-dd28ca754de3
# ╠═802d53a4-777d-11eb-25cc-efc4d9dcf292
# ╠═49df45ee-7782-11eb-2363-d562aa91ca14
# ╠═3f7a0d3e-76ec-11eb-368e-17428098bc72
# ╠═52e23b14-76ec-11eb-1787-01b2cf9888eb
