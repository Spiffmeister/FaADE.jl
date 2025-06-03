push!(LOAD_PATH,"../src/")

using Documenter, Literate,
    FaADE

LitPath1D = joinpath(@__DIR__,"..","tutorials","example1D.jl")
LitPath2D = joinpath(@__DIR__,"..","tutorials","example2D.jl")
LitPathBF = joinpath(@__DIR__,"..","tutorials","example2Dparallel.jl")
LitPathCF = joinpath(@__DIR__,"..","tutorials","example2DCurvilinear.jl")
DocSrc = joinpath(@__DIR__,"src","tutorials") #.md creation path

Literate.markdown(LitPath1D,DocSrc)
Literate.markdown(LitPath2D,DocSrc)
Literate.markdown(LitPathBF,DocSrc)
Literate.markdown(LitPathCF,DocSrc)

makedocs(sitename="FaADE Documentation",
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "tutorials/example1D.md"
            "tutorials/example2D.md"
            "tutorials/example2Dparallel.md"
            "tutorials/example2DCurvilinear.md"
        ],
        "Modules" => [
            "Derivative Operators" => "Derivatives.md",
            "SATs" => "SATs.md",
            "Grids" => "Grid.md",
            "Solvers" => "solvers.md",
            "Solution" => "solution.md",
            "Problems" => "UserInteraction.md",
            "Parallel Operator" => "Parallel.md"
        ]
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [FaADE,FaADE.Derivatives,FaADE.Helpers,FaADE.SATs,FaADE.solvers],
    warnonly = Documenter.except(:linkcheck,:footnote)
    )