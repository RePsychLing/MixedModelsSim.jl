using Documenter, MixedModelsSim

makedocs(;
    modules=[MixedModelsSim],
    format=Documenter.HTML(assets=String[]),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/RePsychLing/MixedModelsSim.jl/blob/{commit}{path}#L{line}",
    sitename="MixedModelsSim.jl",
    authors="Phillip Alday, Douglas Bates, Lisa DeBruine, Reinhold Kliegl",
)

deploydocs(;
    repo="github.com/RePsychLing/MixedModelsSim.jl",
)
