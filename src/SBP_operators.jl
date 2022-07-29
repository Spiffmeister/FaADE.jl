"""
Module for `SBP_operators.jl` package - for computing finite differences by Summation by Parts 
    with Simulatenous Approximation Terms.

# Public functions

- [`Dₓ`](@ref)
- [`Dₓₓ`](@ref)
- [`SAT`](@ref)
- [`SAT_Periodic`](@ref)
- [`Split_domain`](@ref)
- [`time_solver`](@ref)

"""
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
    export Left, Internal, Right, 
        Dirichlet, Neumann, Robin, Periodic, Interface
    export Dₓ, Dₓₓ, Dₓₓ!, SAT, SAT_Periodic, Split_domain, time_solver, build_H

end # module
