using Documenter, MixedModelsSim

makedocs(;
    modules=[MixedModelsSim],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/RePsychLing/MixedModelsSim.jl/blob/{commit}{path}#L{line}",
    sitename="MixedModelsSim.jl",
    authors="Phillip Alday, Douglas Bates, Shravan Vasishth, Reinhold Kliegl",
    assets=String[],
)

deploydocs(;
    repo="github.com/RePsychLing/MixedModelsSim.jl",
)
