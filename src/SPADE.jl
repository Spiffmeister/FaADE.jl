"""
Module for `SPADE.jl` package - for computing finite differences by Summation by Parts 
    with Simulatenous Approximation Terms.

# Submodules

- [`Helpers`](@ref)
- [`Derivatives`](@ref)
- [`SATs`](@ref)
- [`solvers`](@ref)

"""
module SPADE

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
    using SPADE.Helpers: Dirichlet, Neumann, Robin, Periodic, Interface
    using SPADE.Helpers: Left, Internal, Right, Up, Down
    using SPADE.Helpers: Grid1D, Grid2D
    using SPADE.Helpers: Boundary, PeriodicBoundary
    using SPADE.Helpers: VariableCoefficientPDE1D, VariableCoefficientPDE2D

    using SPADE.Derivatives: D₁, D₂, D₂!
    
    using SPADE.SATs: SAT_Periodic, SAT_Dirichlet, SAT_Neumann, SimultanousApproximationTerm

    using SPADE.solvers: solve, build_H

    using SPADE.Parallel: construct_grid, generate_parallel_penalty

    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right, Up, Down
    export Grid1D, Grid2D, Boundary, PeriodicBoundary
    export VariableCoefficientPDE1D, VariableCoefficientPDE2D
    export construct_grid, generate_parallel_penalty

    export D₁, D₂, D₂!

    export solve, build_H
    
    export SimultanousApproximationTerm,
        SAT_Dirichlet,
        SAT_Neumann,
        SAT_Periodic
        #SAT_Robin, 
        #Split_domain

end # module
