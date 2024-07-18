"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using FaADE.Helpers

    using FaADE.Derivatives: mul!,
    innerH,
    DerivativeOperatorType, DiffusionOperator, DiffusionOperatorND

    using FaADE.Grid
    
    using FaADE.SATs: SAT_Periodic, SAT_Periodic!,
    SAT_Dirichlet, SAT_Dirichlet_explicit!, SAT_Dirichlet_solution!, SAT_Dirichlet_data!,
    SAT_Neumann, SAT_Neumann_data!, SAT_Neumann_solution!,
    SATpenalties, SimultanousApproximationTerm, SAT_Interface, SAT_Interface!, SAT_Interface_cache!,
    SATMode, DataMode, SolutionMode, ExplicitMode

    using FaADE.ParallelOperator: ParallelGrid, ParallelMultiBlock, ParallelData, applyParallelPenalty!, compute_parallel_operator, computeglobalw!

    using FaADE.Inputs: newPDEProblem, Problem1D, Problem2D#, 

    
    include("DataStorage.jl")
    include("IntegeratorData.jl")
    include("SettingValues.jl")
    include("solution.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end