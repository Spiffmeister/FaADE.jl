"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using FaADE.Helpers

    using FaADE.Derivatives: DerivativeOrder, GetOrder, DerivativeOperator, mul!,
    innerH

    using FaADE.Grid
    
    using FaADE.SATs: SAT_Periodic, SAT_Periodic!,
    SAT_Dirichlet, SAT_Dirichlet_explicit!, SAT_Dirichlet_solution!, SAT_Dirichlet_data!,
    SAT_Neumann, SAT_Neumann_data!, SAT_Neumann_solution!,
    SATpenalties, SimultanousApproximationTerm, SAT_Interface, SAT_Interface!, 
    SATMode, DataMode, SolutionMode, ExplicitMode

    using FaADE.ParallelOperator: ParallelGrid, generate_parallel_penalty, ParallelData, applyParallelPenalty!, compute_parallel_operator

    using FaADE.Inputs: newPDEProblem, newProblem1D, newProblem2D#, 

    
    include("DataStorage.jl")
    include("IntegeratorData.jl")
    include("SettingValues.jl")
    include("solution.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end