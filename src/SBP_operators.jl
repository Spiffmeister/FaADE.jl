module SBP_operators

    # Inbuild julia packages
    using LinearAlgebra


    # Include the files (add as completed)
    include("Op_Dx.jl")
    include("Op_Dxx.jl")
    include("BoundaryOp_Dxx.jl")
    include("TimeSolver.jl")


    # Export the functions for direct user interaction
    export Dₓ, Dₓₓ, Dₓ!, Dₓₓ!, SAT_left, SAT_right, SAT_Periodic, Split_domain, time_solver

end # module
