using Pkg
using LiveServer

Pkg.activate("docs")
servedocs(skip_dirs=[joinpath("docs","src","tutorials")])


servedocs(literate=joinpath("docs","src","tutorials"))
