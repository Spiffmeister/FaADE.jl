using Pkg
Pkg.activate("..")
push!(LOAD_PATH,"../src/")

using Documenter, Literate,
    SBP_operators,
    SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers

LitPath1D = joinpath(@__DIR__,"..","tutorials","example1D.jl")
LitPath2D = joinpath(@__DIR__,"..","tutorials","example2D.jl")
LitPathBF = joinpath(@__DIR__,"..","tutorials","example3D.jl")
DocSrc = joinpath(@__DIR__,"src")

Literate.markdown(LitPath1D,DocSrc)
Literate.markdown(LitPath2D,DocSrc)
Literate.markdown(LitPathBF,DocSrc)

makedocs(sitename="SBP Operators Documentation",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "1D Example" => "example1D.md"
            "2D Example" => "example2D.md"
            "3D Example" => "exampleBField.md"
        ],
        "Derivative Operators" => "Derivatives.md",
        "SATs" => "SATs.md",
        "Solvers" => "solvers.md",
        "Helpers" => "Helpers.md"
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [SBP_operators,SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers]
    )