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
    using SharedArrays
    using LinearAlgebra
    using JLD2




    # Include the files (add as completed)
    include("types/types.jl")
    using SBP_operators.types: Dirichlet, Neumann, Robin, Periodic, Interface
    using SBP_operators.types: Left, Internal, Right

    include("Derivatives/Derivatives.jl")
    using SBP_operators.Derivatives: Dₓ, Dₓₓ, Dₓₓ!, Stencil2D
    
    include("SATs/SATs.jl")
    using SBP_operators.SATs: SAT, SAT_Periodic

    include("solvers/solvers.jl")
    using SBP_operators.solvers: time_solver, build_H


    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right

    export Dₓ, Dₓₓ, Dₓₓ!, Stencil2D

    export time_solver, build_H
    
    export SAT, SAT!, 
        SAT_Dirichlet, SAT_Dirichlet!, 
        SAT_Neumann, SAT_Neumann!, 
        SAT_Periodic, SAT_Periodic!, 
        SAT_Robin, 
        Split_domain

end # module
