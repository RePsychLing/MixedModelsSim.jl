using Documenter
using MixedModelsSim

makedocs(; sitename="MixedModelsSim.jl",
         doctest=true,
         authors="Phillip Alday, Douglas Bates, Lisa DeBruine, Reinhold Kliegl",
         pages=["index.md",
                "Rapid Start" => "simulation.md",
                "Beginner Friendly Tutorial" => "simulation_tutorial.md"])

deploydocs(; repo="github.com/RePsychLing/MixedModelsSim.jl", push_preview=true)
