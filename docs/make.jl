using Thermochron
using Documenter

DocMeta.setdocmeta!(Thermochron, :DocTestSetup, :(using Thermochron); recursive=true)

makedocs(;
    modules=[Thermochron],
    authors="C. Brenhin Keller",
    repo="https://github.com/brenhinkeller/Thermochron.jl/blob/{commit}{path}#{line}",
    sitename="Thermochron.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brenhinkeller.github.io/Thermochron.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brenhinkeller/Thermochron.jl",
)
