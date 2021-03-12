using Documenter
using MixedModelsSim

makedocs(
    modules = [MixedModelsSim],
    root = joinpath(dirname(pathof(MixedModelsSim)), "..", "docs"),
    sitename = "MixedModelsSim.jl",
    doctest = true,
    authors = "Phillip Alday, Douglas Bates, Lisa DeBruine, Reinhold Kliegl",
    pages = ["Home" => "index.md"
             "Full Example" => "simulation.md"],
)

deploydocs(; repo = "github.com/RePsychLing/MixedModelsSim.jl", push_preview = true)
