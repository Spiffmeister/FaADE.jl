module solvers

    using Distributed 
    using SharedArrays
    using LinearAlgebra



    using SBP_operators.types
    using SBP_operators.SATs: SAT!, SAT_Periodic!


    include("solution.jl")
    include("timesolvers.jl")
    include("integrators.jl")

    export time_solver, build_H

end