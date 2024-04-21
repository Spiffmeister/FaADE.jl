"""
    Helpers
Module contains many of the data structures for users and internal use, as well as nifty functions.
"""
module Helpers

    include("types.jl")
    include("nifty.jl")

    # include("gridconstruction.jl")
    # include("grid.jl")
    # include("UserTypes.jl")

    include("HMatrix.jl")
    # include("DataStorage.jl")


    # types import
    export NodeType, Left, Internal, Right, Up, Down
    export _flip
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    # export SATMode, DataMode, SolutionMode
    export SourceTerm, DiffusionCoefficient
    export BoundaryStorage

    # nifty import
    export check_order, check_boundary,
        halforder, BoundaryNodeInput, BoundaryNodeOutput,
        SelectLoopDirection, SATNodeOutput
    
    # grid import
    # export GridType, LocalGridType
    # export Grid1D, Grid2D, GridMultiBlock
    # export CartesianMetric, CurvilinearMetric
    # export eachgrid


    # # Export UserTypes
    # export BoundaryConditionData, Boundary, PeriodicBoundary,
    # PDEProblem, VariableCoefficientPDE1D, VariableCoefficientPDE2D,BoundaryConditions
    # # DataStorage import
    # export DataBlockType, BoundaryStorage, 
    #     DataBlock,
    #     BoundaryData1D, BoundaryData2D, 
    #     copyUtoSAT!, copySATtoU!, addSATtoU!,
    #     addSource!, BoundaryConditions
    
        
    # export build_H, innerH
    
    # PostNifty import
    export GetAxis, GetDim,
        GetMinΔ
    
end

