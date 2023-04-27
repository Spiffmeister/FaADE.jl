"""
Module for `SBP_operators.jl` package - for computing finite differences by Summation by Parts 
    with Simulatenous Approximation Terms.

# Submodules

- [`Helpers`](@ref)
- [`Derivatives`](@ref)
- [`SATs`](@ref)
- [`solvers`](@ref)

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
    include("Parallel/Parallel.jl")
    include("solvers/solvers.jl")


    # Include the files (add as completed)
    # Helpers export
    using SBP_operators.Helpers: Dirichlet, Neumann, Robin, Periodic, Interface
    using SBP_operators.Helpers: Left, Internal, Right, Up, Down
    using SBP_operators.Helpers: Grid1D, Grid2D
    using SBP_operators.Helpers: Boundary, PeriodicBoundary
    using SBP_operators.Helpers: VariableCoefficientPDE1D, VariableCoefficientPDE2D

    using SBP_operators.Derivatives: Dₓ, Dₓₓ, Dₓₓ!
    
    using SBP_operators.SATs: SAT_Periodic, SAT_Dirichlet, SAT_Neumann, SimultanousApproximationTerm

    using SBP_operators.solvers: solve, build_H

    using SBP_operators.Parallel: construct_grid, generate_parallel_penalty

    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right, Up, Down
    export Grid1D, Grid2D, Boundary, PeriodicBoundary
    export VariableCoefficientPDE1D, VariableCoefficientPDE2D
    export construct_grid, generate_parallel_penalty

    export Dₓ, Dₓₓ, Dₓₓ!

    export solve, build_H
    
    export SimultanousApproximationTerm,
        SAT_Dirichlet,
        SAT_Neumann,
        SAT_Periodic
        #SAT_Robin, 
        #Split_domain

end # module
