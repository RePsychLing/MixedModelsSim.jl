module MixedModelsSim

using DataFrames
using MixedModels
using PooledArrays
using Random
using StaticArrays
using Statistics


export
    cyclicshift,
    factorproduct,
    itemsubjdf,
    nlevels,
    pooled!,
    simulate_waldtests,
    withinitem,
    power_table,
    sim_to_df,
    simdat_crossed

include("columntable.jl")
include("power.jl")
include("simdat.jl")

end # module
