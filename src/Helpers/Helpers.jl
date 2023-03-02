"""
    Helpers
Module contains many of the data structures for users and internal use, as well as internal helper functions.
"""
module Helpers

    include("types.jl")
    include("nifty.jl")

    include("grid.jl")
    include("UserTypes.jl")

    include("ConjGradData.jl")
    include("DataStorage.jl")

    include("postnifty.jl")

    include("penaltyfn.jl")



    # types import
    export NodeType, Left, Internal, Right, Up, Down
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    export SATMode, DataMode, SolutionMode

    # nifty import
    export check_order, check_boundary,
        halforder, BoundaryNodeInput, BoundaryNodeOutput,
        SelectLoopDirection
    
    # grid import
    export GridType, Grid1D, Grid2D

    # Export UserTypes
    export BoundaryConditionData, Boundary, PeriodicBoundary,
    PDEProblem, VariableCoefficientPDE1D, VariableCoefficientPDE2D
    
    # DataStorage import
    export DataBlockType, BoundaryStorage, 
        DataBlock,
        BoundaryData1D, BoundaryData2D, 
        copyUtoSAT!, copySATtoU!, addSATtoU!,
        addSource!, setBoundary!
    
    export ConjGradBlock, build_H, ExplicitBlock
    
    # PostNifty import
    export GetAxis, GetDim,
        GetMinΔ
    
    # Parallel penalty functions
    export ParallelGrid, generate_parallel_penalty
    
end

