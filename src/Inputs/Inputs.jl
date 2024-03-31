module Inputs

    using FaADE.Helpers
    using FaADE.Grid
    using FaADE.SATs
    using FaADE.Derivatives
    using FaADE.ParallelOperator

    # include("oldUserTypes.jl")
    # include("oldDataStorage.jl")
    # include("SAT_Interface.jl")
    include("UserTypes.jl")

    export newPDEProblem
    export Problem1D, Problem2D

    export SATBoundaries

    # oldUserTypes
    # export BoundaryConditionData
    # export Boundary, PeriodicBoundary,
    # PDEProblem, VariableCoefficientPDE1D, VariableCoefficientPDE2D,BoundaryConditions
    # export SAT_Interface, SAT_Interface!
    # export DataBlockType, BoundaryStorage, 
    #     DataBlock,
    #     BoundaryData1D, BoundaryData2D, 
    #     copyUtoSAT!, copySATtoU!, addSATtoU!,
    #     addSource!, BoundaryConditions


end