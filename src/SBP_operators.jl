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


    # Includes
    include("Helpers/Helpers.jl")
    include("Derivatives/Derivatives.jl")
    include("SATs/SATs.jl")
    include("solvers/solvers.jl")


    # Include the files (add as completed)
    # Helpers export
    using SBP_operators.Helpers: Dirichlet, Neumann, Robin, Periodic, Interface
    using SBP_operators.Helpers: Left, Internal, Right
    using SBP_operators.Helpers: Grid1D, Grid2D
    using SBP_operators.Helpers: Boundary, PeriodicBoundary
    using SBP_operators.Helpers: VariableCoefficientPDE1D, VariableCoefficientPDE2D

    using SBP_operators.Derivatives: Dₓ, Dₓₓ, Dₓₓ!
    
    using SBP_operators.SATs: SAT_Periodic, SATDirichlet, SimultanousApproximationTerm

    using SBP_operators.solvers: solve, build_H


    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right
    export Grid1D, Grid2D, Boundary, PeriodicBoundary
    export VariableCoefficientPDE1D, VariableCoefficientPDE2D

    export Dₓ, Dₓₓ, Dₓₓ!

    export solve, build_H
    
    export SimultanousApproximationTerm,
        SAT_Dirichlet, SAT_Dirichlet!, SATDirichlet,
        SAT_Periodic, SAT_Periodic!, 
        Split_domain
        #SAT_Neumann, SAT_Neumann!, 
        #SAT_Robin, 

end # module
