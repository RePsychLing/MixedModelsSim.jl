module MixedModelsSim

using DataFrames
using Distributions: Chisq, ccdf
using MixedModels
using PooledArrays
using Random
using Tables
using Statistics

export
    cyclicshift,
    factorproduct,
    nlevels,
    pooled!,
    simulate_waldtests,
    withinitem,
    power_table,
    simdat_crossed

include("columntable.jl")
include("power.jl")
include("simdat.jl")

end # module
