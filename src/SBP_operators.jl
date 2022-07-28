module SBP_operators

    # Inbuild julia packages
    using Distributed
    using LinearAlgebra


    # Include the files (add as completed)
    include("types.jl")
    include("Op_Dx.jl")
    include("Op_Dxx.jl")
    include("SAT_SecondOrder.jl")
    include("TimeSolver.jl")
    include("integrators.jl")


    # Export the functions for direct user interaction
    export NodeLeft, NodeInternal, NodeRight, 
        Dirichlet, Neumann, Robin, Periodic, Interface
    export Dₓ, Dₓₓ, Dₓ!, Dₓₓ!, SAT_left, SAT_right, SAT_Periodic, Split_domain, time_solver, build_H

end # module
