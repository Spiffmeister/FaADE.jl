"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using FaADE.Helpers

    using FaADE.Derivatives: generate_SecondDerivative
    
    using FaADE.SATs: SAT, construct_SAT, 
    SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, SAT_Dirichlet_implicit_data!

    using FaADE.Parallel: ParallelGrid, generate_parallel_penalty
    

    include("solution.jl")
    include("IntegeratorData.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end