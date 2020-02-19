module MixedModelsSim

using DataFrames
using PooledArrays
using Tables
using Statistics
using Pkg.Artifacts

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

function __init__()
    global TestData = artifact"TestData"
end

include("MMutils.jl")
include("columntable.jl")
include("power.jl")
include("simdat.jl")

end # module
