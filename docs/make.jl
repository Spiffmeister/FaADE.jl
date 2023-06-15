push!(LOAD_PATH,"../src/")

using Documenter, Literate,
    FaADE

LitPath1D = joinpath(@__DIR__,"..","tutorials","example1D.jl")
LitPath2D = joinpath(@__DIR__,"..","tutorials","example2D.jl")
LitPathBF = joinpath(@__DIR__,"..","tutorials","example3D.jl")
DocSrc = joinpath(@__DIR__,"src","tutorials") #.md creation path

Literate.markdown(LitPath1D,DocSrc)
Literate.markdown(LitPath2D,DocSrc)
Literate.markdown(LitPathBF,DocSrc)

makedocs(sitename="SBP Operators Documentation",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "tutorials/example1D.md"
            "tutorials/example2D.md"
            "tutorials/example3D.md"
        ],
        "Modules" => [
            "Derivative Operators" => "Derivatives.md",
            "SATs" => "SATs.md",
            "Solvers" => "solvers.md",
            "Helpers" => "Helpers.md"
        ]
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [FaADE,FaADE.Derivatives,FaADE.Helpers,FaADE.SATs,FaADE.solvers]
    )