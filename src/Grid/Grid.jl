"""
    Grid

Submodule containing the grid data structure and functions.
"""
module Grid

    using FaADE.Helpers
    # using Derivatives: DerivativeOrder, D₁!

    include("DerivativeFirst.jl")
    include("toroidal.jl")
    include("generating.jl")
    include("grids.jl")
    # include("data.jl")

    export GridType, LocalGridType
    export CartesianMetric, CurvilinearMetric
    export coordtype

    export Joint

    export Grid1D, Grid2D, GridMultiBlock
    export Torus
    export meshgrid
    
    export eachgrid

end