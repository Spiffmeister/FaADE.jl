using Pkg
Pkg.activate("..")
push!(LOAD_PATH,"../src/")

using Documenter, SBP_operators


makedocs(sitename="SBP Operators Documentation",
    format=Documenter.HTML(prettyurls=false)
    )