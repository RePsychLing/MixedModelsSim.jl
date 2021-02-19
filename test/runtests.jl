using MixedModelsSim
using Test

include("utilities.jl")
include("simdat.jl")
# there's no point in testing the big power stuff
# until all the stuff it depends on has been tested
include("power.jl")
