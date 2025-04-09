"""
    ParallelOperator
Methods for constructing the parallel operator
"""
module ParallelOperator

    using Base.Threads

    using DifferentialEquations: ODEProblem, EnsembleProblem, EnsembleSerial, solve, Tsit5, remake
    # using Interpolations: LinearInterpolation
    using Interpolations
    using BasicInterpolators: BicubicInterpolator, LinearInterpolator
    using CubicHermiteSpline

    using LinearAlgebra: norm

    using FaADE.Helpers
    using FaADE.Derivatives: build_H
    using FaADE.Derivatives: MassMatrix, DiagonalH, CompositeH
    using FaADE.Derivatives: D₁!

    using FaADE.Grid

    include("types.jl")
    include("ParallelData.jl")
    include("ConstructParallelGrid.jl")
    include("penaltyfn.jl")

    export ParallelGrid
    export ParallelData, ParallelMultiBlock

    export construct_grid
    export applyParallelPenalty!, computeglobalw!
    export compute_parallel_operator
    export MagneticField

end