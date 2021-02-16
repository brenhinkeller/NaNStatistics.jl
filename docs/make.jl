using NaNStatistics
using Documenter

DocMeta.setdocmeta!(NaNStatistics, :DocTestSetup, :(using NaNStatistics); recursive=true)

makedocs(;
    modules=[NaNStatistics],
    authors="C. Brenhin Keller",
    repo="https://github.com/brenhinkeller/NaNStatistics.jl/blob/{commit}{path}#{line}",
    sitename="NaNStatistics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brenhinkeller.github.io/NaNStatistics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brenhinkeller/NaNStatistics.jl",
)
