"""
    Grid

Submodule containing the grid data structure and functions.
"""
module Grid

    # using FaADE.Helpers
    # using Derivatives: DerivativeOrder, D₁!

    include("grids.jl")
    include("toroidal.jl")
    # include("data.jl")

    export GridType, LocalGridType
    export CartesianMetric, CurvilinearMetric

    export Grid1D, Grid2D, GridMultiBlock
    export Torus
    export meshgrid
    
    export eachgrid

end