"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using FaADE.Helpers

    using FaADE.Derivatives: generate_SecondDerivative, DerivativeOrder, GetOrder, DerivativeOperator, mul!,
    innerH

    using FaADE.Grid
    
    using FaADE.SATs: SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, SAT_Dirichlet_implicit_data!, SimultanousApproximationTerm, SAT_Interface, SAT_Interface!, SATMode, DataMode, SolutionMode#,
    # SATBoundaries, SAT, construct_SAT

    using FaADE.ParallelOperator: ParallelGrid, generate_parallel_penalty, ParallelData, applyParallelPenalty!

    using FaADE.Inputs: newPDEProblem, newProblem1D, newProblem2D, 
    BoundaryConditionData, Boundary, PeriodicBoundary,
    PDEProblem, VariableCoefficientPDE1D, VariableCoefficientPDE2D,BoundaryConditions,
    # old imports
    DataBlockType, BoundaryStorage, DataBlock, BoundaryData1D, BoundaryData2D, SAT#, copyUtoSAT!, copySATtoU!, addSATtoU!,addSource!, BoundaryConditions

    
    include("DataStorage.jl")
    include("IntegeratorData.jl")
    include("SettingValues.jl")
    include("solution.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end