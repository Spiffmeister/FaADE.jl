module solvers

    using Distributed 
    using SharedArrays
    using LinearAlgebra



    using SBP_operators.Helpers
    
    using SBP_operators.SATs: SAT_Periodic, SAT_Periodic!, SAT_Dirichlet_implicit!, SATpenalties, BDₓᵀ, SAT_Dirichlet_implicit_data!, 
    SimultanousApproximationTermContainer, construct_SATs

    include("solution.jl")
    include("timesolvers.jl")
    include("integrators.jl")

    export time_solver, build_H

end