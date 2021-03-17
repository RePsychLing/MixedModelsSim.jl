using Documenter
using MixedModelsSim

makedocs(
    root = joinpath(dirname(pathof(MixedModelsSim)), "..", "docs"),
    sitename = "MixedModelsSim.jl",
    doctest = true,
    authors = "Phillip Alday, Douglas Bates, Lisa DeBruine, Reinhold Kliegl",
    pages = [
        "index.md",
        "simulation.md",
        "simulation_tutorial.md",
    ],
)

deploydocs(; repo = "github.com/RePsychLing/MixedModelsSim.jl", push_preview = true)
