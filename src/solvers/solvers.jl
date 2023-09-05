"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using FaADE.Helpers

    using FaADE.Derivatives: generate_SecondDerivative, DerivativeOrder, GetOrder, DerivativeOperator
    
    using FaADE.SATs: SAT, construct_SAT, 
    SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, SAT_Dirichlet_implicit_data!, SATBoundaries, SimultanousApproximationTerm, SAT_Interface

    using FaADE.ParallelOperator: ParallelGrid, generate_parallel_penalty, ParallelData, applyParallelPenalty!

    using FaADE.Inputs: newPDEProblem, newProblem1D, newProblem2D
    
    include("DataStorage.jl")
    include("IntegeratorData.jl")
    include("SettingValues.jl")
    include("solution.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end