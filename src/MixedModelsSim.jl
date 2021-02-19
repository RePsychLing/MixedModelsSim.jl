module MixedModelsSim

using MixedModels
using PooledArrays
using Random
using Statistics
using Tables

using MixedModels: replicate

export
    cyclicshift,
    factorproduct,
    nlevels,
    #withinitem,
    pooled!,
    power_table,
    simdat_crossed

include("utilities.jl")
include("power.jl")
include("simdat.jl")

end # module
