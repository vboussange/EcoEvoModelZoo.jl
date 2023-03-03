using EcoEvoModelZoo
using Documenter
import ParametricModels.AbstractModel

DocMeta.setdocmeta!(EcoEvoModelZoo, :DocTestSetup, :(using EcoEvoModelZoo); recursive=true)

makedocs(;
    modules=[EcoEvoModelZoo],
    authors="Victor Boussange",
    repo="https://github.com/vboussange/EcoEvoModelZoo.jl/blob/{commit}{path}#{line}",
    sitename="EcoEvoModelZoo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vboussange.github.io/EcoEvoModelZoo.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vboussange/EcoEvoModelZoo.jl",
    devbranch="main",
)
