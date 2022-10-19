using Pkg
Pkg.activate("..")
push!(LOAD_PATH,"../src/")

using Documenter, Literate,
    SBP_operators,
    SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers

LitPath1D = joinpath(@__DIR__,"..","tutorials","example1D.jl")
LitPath2D = joinpath(@__DIR__,"..","tutorials","example2D.jl")
DocSrc = joinpath(@__DIR__,"src")

Literate.markdown(LitPath1D,DocSrc)

makedocs(sitename="SBP Operators Documentation",
    pages = [
        "Home" => "index.md"
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [SBP_operators,SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers]
    )