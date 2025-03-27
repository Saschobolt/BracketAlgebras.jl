using BracketAlgebras
using Documenter

DocMeta.setdocmeta!(BracketAlgebras, :DocTestSetup, :(using BracketAlgebras, AbstractAlgebra); recursive=true)

makedocs(;
    modules=[BracketAlgebras],
    authors="Sascha Stüttgen <sascha.stuettgen@rwth-aachen.de> and contributors",
    sitename="BracketAlgebras.jl",
    format=Documenter.HTML(;
        canonical="https://Saschobolt.github.io/BracketAlgebras.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "The Bracket Algebra" => "bracket_algebra.md",
    ],
)

deploydocs(;
    repo="github.com/Saschobolt/BracketAlgebras.jl",
    devbranch="master",
)
