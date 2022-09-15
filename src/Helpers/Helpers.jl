module Helpers

    include("types.jl")
    include("nifty.jl")

    include("grid.jl")
    include("UserTypes.jl")

    include("DataStorage.jl")



    # types import
    export NodeType, Left, Internal, Right
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    export SATMode, DataMode, SolutionMode

    # nifty import
    export check_order, check_boundary,
        halforder, BoundaryNodeInput, BoundaryNodeOutput

    # grid import
    export GridType, Grid1D, Grid2D

    # Export UserTypes
    export BoundaryConditionData, Boundary, PeriodicBoundary,
    VariableCoefficientPDE1D, VariableCoefficientPDE2D
    
    # DataStorage import
    export DataBlockType, BoundaryStorage, 
        DataBlock, ConjGradBlock,
        BoundaryData1D, BoundaryData2D, 
        copyUtoSAT
    
end

