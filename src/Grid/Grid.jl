"""
    Grid

Submodule containing the grid data structure and functions.
"""
module Grid

    using FaADE.Helpers

    include("DerivativeFirst.jl")
    include("toroidal.jl")
    include("generating.jl")
    include("grids.jl")
    # include("data.jl")

    export GridType, LocalGridType
    export CartesianMetric, CurvilinearMetric
    export coordtype
    export GetMetricType

    export Joint

    export Grid1D, Grid2D, GridMultiBlock
    export Torus
    export meshgrid
    
    export GetBoundaryCoordinates
    export eachgrid

end