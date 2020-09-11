using MixedModelsSim
using Test

include("levels.jl")
include("product.jl")
include("sim.jl")
include("simdat.jl")
# there's no point in testing the big power stuff 
# until all the stuff it depepnds on has been tested
include("power.jl")

