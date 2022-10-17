"""
    solvers
Module containting the time and spatial solvers.

WARNING: IMPLICIT METHODS ARE DEPRECATED
"""
module solvers

    import LinearAlgebra: dot, norm

    using SBP_operators.Helpers

    using SBP_operators.Derivatives: Dₓₓ!, generate_Derivative
    
    using SBP_operators.SATs: SAT, construct_SAT, 
    SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, BDₓᵀ, SAT_Dirichlet_implicit_data!
    

    include("solution.jl")
    include("timesolvers.jl")
    include("integrators.jl")

    export time_solver, build_H, solve

end