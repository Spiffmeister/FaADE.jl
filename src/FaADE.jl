"""
Module for `FaADE.jl` package - for computing finite differences by Summation by Parts 
    with Simulatenous Approximation Terms.

# Submodules

- [`Helpers`](@ref)
- [`Derivatives`](@ref)
- [`SATs`](@ref)
- [`solvers`](@ref)

"""
module FaADE

    # Inbuild julia packages
    using Distributed
    using SharedArrays
    using LinearAlgebra
    using JLD2


    # Includes
    include("Helpers/Helpers.jl")
    include("Derivatives/Derivatives.jl")
    include("Grid/Grid.jl")
    include("SATs/SATs.jl")
    include("ParallelOperator/ParallelOperator.jl")
    include("Inputs/Inputs.jl")
    include("solvers/solvers.jl")


    # Include the files (add as completed)
    # Helpers export
    using FaADE.Helpers: Dirichlet, Neumann, Robin, Periodic, Interface
    using FaADE.Helpers: Left, Internal, Right, Up, Down
    # using FaADE.Helpers: Grid1D, Grid2D, GridMultiBlock, CartesianMetric, CurvilinearMetric
    # using FaADE.Helpers: Boundary, PeriodicBoundary
    # using FaADE.Helpers: VariableCoefficientPDE1D, VariableCoefficientPDE2D

    using FaADE.Derivatives: D₁, D₂!

    using FaADE.Grid: Grid1D, Grid2D, GridMultiBlock, CartesianMetric, CurvilinearMetric, Joint
    
    using FaADE.SATs: SAT_Periodic, SAT_Dirichlet, SAT_Neumann, SimultanousApproximationTerm

    using FaADE.solvers: solve

    using FaADE.ParallelOperator: construct_grid, generate_parallel_penalty, ParallelData

    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right, Up, Down
    export Grid1D, Grid2D, GridMultiBlock, Joint
    export Boundary, PeriodicBoundary, CartesianMetric, CurvilinearMetric
    export VariableCoefficientPDE1D, VariableCoefficientPDE2D
    export construct_grid, generate_parallel_penalty

    export D₁, D₂, D₂!

    export solve#, build_H
    
    export SimultanousApproximationTerm,
        SAT_Dirichlet,
        SAT_Neumann,
        SAT_Periodic
        #SAT_Robin, 
        #Split_domain

    using FaADE.Inputs: newProblem1D, newProblem2D,
    #old using
    Boundary,PeriodicBoundary, VariableCoefficientPDE1D, VariableCoefficientPDE2D
    export newProblem1D, newProblem2D
    


end # module
