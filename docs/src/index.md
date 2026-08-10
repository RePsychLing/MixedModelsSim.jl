# MixedModelsSim.jl

MixedModelsSim.jl provides utilities for designing and power-analysing experiments with mixed-effects models.
It is a companion to [MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/).

## What this package does

- **Generate crossed experimental designs** — subjects × items × conditions with arbitrary combinations of between-subject, between-item, and within-both factors ([`simdat_crossed`](@ref))
- **Specify random-effects covariance structure** — create the lower Cholesky factor of a covariance matrix from standard deviations and optional correlations ([`create_re`](@ref)), and convert it to the θ parameter vector MixedModels.jl expects ([`createθ`](@ref))
- **Run power analyses** — simulate data via `parametricbootstrap` from MixedModels.jl and summarise the results as a power table ([`power_table`](@ref))

## Installation

MixedModelsSim.jl is registered in the Julia package registry:

```julia
using Pkg
Pkg.add("MixedModelsSim")
```

## Getting started

| Guide | Description |
|:------|:------------|
| [Rapid Start](@ref "Simulating an Experiment from Scratch") | A concise worked example: build a design from scratch, specify parameters, simulate, and compute power. |
| [Tutorial](@ref "Power Analysis and Simulation Tutorial") | A step-by-step tutorial covering four scenarios: bootstrapping from existing data, adapting effect sizes, creating a design from scratch, and power curves over varying sample sizes. |
| [API Reference](@ref "API") | Complete documentation for all exported functions. |

