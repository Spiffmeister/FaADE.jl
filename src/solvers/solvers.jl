"""
    solvers
Module containting the time and spatial solvers.

WARNING: EXPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using SBP_operators.Helpers

    using SBP_operators.Derivatives: generate_SecondDerivative
    
    using SBP_operators.SATs: SAT, construct_SAT, 
    SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, SAT_Dirichlet_implicit_data!

    using SBP_operators.Parallel: ParallelGrid
    

    include("solution.jl")
    include("IntegeratorData.jl")
    include("integrators.jl")
    include("timesolvers.jl")

    export solve

end