using Pkg
Pkg.activate("..")
push!(LOAD_PATH,"../src/")

using Documenter, SBP_operators,
    SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers


makedocs(sitename="SBP Operators Documentation",
    pages = [
        "Home" => "index.md"
    ],
    format=Documenter.HTML(prettyurls=false),
    modules = [SBP_operators,SBP_operators.Derivatives,SBP_operators.Helpers,SBP_operators.SATs,SBP_operators.solvers]
    )