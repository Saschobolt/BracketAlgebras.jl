using BracketAlgebra
using Documenter

DocMeta.setdocmeta!(BracketAlgebra, :DocTestSetup, :(using BracketAlgebra); recursive=true)

makedocs(;
    modules=[BracketAlgebra],
    authors="Sascha Stüttgen <sascha.stuettgen@rwth-aachen.de> and contributors",
    sitename="BracketAlgebra.jl",
    format=Documenter.HTML(;
        canonical="https://Saschobolt.github.io/BracketAlgebra.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Saschobolt/BracketAlgebra.jl",
    devbranch="master",
)
