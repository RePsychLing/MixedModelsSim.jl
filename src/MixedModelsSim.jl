module MixedModelsSim

using LinearAlgebra
using MixedModels
using PooledArrays
using PrettyTables
using Random
using Statistics
using Tables

using MixedModels: replicate

export
       create_re,
       create_theta,
       createθ,
       cyclicshift,
       factorproduct,
       flatlowertri,
       nlevels,
       nlevstbl,
#withinitem,
       pooled!,
       power_table,
       simdat_crossed,
       update!

export pretty_table, @pt # re-exports

include("utilities.jl")
include("power.jl")
include("simdat.jl")

end # module
