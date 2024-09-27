"""
    solvers
Module containting the time and spatial solvers, data structures require dfor solving, and all functions for data passing during solve.
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
    SAT_Robin, SAT_Robin_data!, SAT_Robin_solution!,
    SAT_Interface, SAT_Interface!, SAT_Interface_cache!, SATpenalties, SimultanousApproximationTerm,  
    SATMode, DataMode, SolutionMode, ExplicitMode

    using FaADE.ParallelOperator: ParallelGrid, ParallelMultiBlock, ParallelData, applyParallelPenalty!, compute_parallel_operator, computeglobalw!

    using FaADE.Inputs: PDEProblem, Problem1D, Problem2D#, 

    
    include("DataStorage.jl")
    include("IntegeratorData.jl")
    include("SettingValues.jl")
    include("solution.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end