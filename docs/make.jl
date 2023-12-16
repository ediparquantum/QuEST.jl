using QuEST
using Documenter

DocMeta.setdocmeta!(QuEST, :DocTestSetup, :(using QuEST); recursive=true)

makedocs(;
    modules=[QuEST],
    authors="Jonathan Miller, jonathan.miller@fieldofnodes.com",
    repo="https://github.com/fieldofnodes/QuEST.jl/blob/{commit}{path}#{line}",
    sitename="QuEST.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fieldofnodes.github.io/QuEST.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fieldofnodes/QuEST.jl",
    devbranch="main",
)
