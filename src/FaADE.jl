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
    # using JLD2


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

    using FaADE.Derivatives: D₁, D₂, D₁!, D₂!

    using FaADE.Grid: Grid1D, Grid2D, GridMultiBlock, CartesianMetric, CurvilinearMetric, Joint
    
    using FaADE.SATs: SAT_Periodic, SAT_Dirichlet, SAT_Neumann, SimultanousApproximationTerm

    using FaADE.solvers: solve

    using FaADE.ParallelOperator: construct_grid, ParallelData

    # Export the functions for direct user interaction
    export Dirichlet, Neumann, Robin, Periodic, Interface
    export Left, Internal, Right, Up, Down
    export Grid1D, Grid2D, GridMultiBlock, Joint
    export CartesianMetric, CurvilinearMetric
    export construct_grid, ParallelData

    export D₁, D₂, D₂!, D₁!

    export solve
    
    export SimultanousApproximationTerm,
        SAT_Dirichlet,
        SAT_Neumann,
        SAT_Periodic
        #SAT_Robin, 
        #Split_domain

    using FaADE.Inputs: Problem1D, Problem2D, SATBoundaries



    export Problem1D, Problem2D, SATBoundaries
    


end # module
