module solvers

    using Distributed 
    using SharedArrays
    using LinearAlgebra



    using SBP_operators.types
    using SBP_operators.SATs: SAT, SAT!, SAT_Periodic, SAT_Periodic!, SATAdd!, SAT_Dirichlet_internal!, SATpenalties, BDₓᵀ, SAT_Dirichlet_internal_forcing!


    include("solution.jl")
    include("timesolvers.jl")
    include("integrators.jl")

    export time_solver, build_H

end