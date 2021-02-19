module MixedModelsSim

using LinearAlgebra
using MixedModels
using PooledArrays
using Random
using Statistics
using Tables

using MixedModels: replicate

export
    cyclicshift,
    factorproduct,
    flatlowertri,
    nlevels,
    #withinitem,
    pooled!,
    power_table,
    simdat_crossed

include("utilities.jl")
include("power.jl")
include("simdat.jl")

end # module
