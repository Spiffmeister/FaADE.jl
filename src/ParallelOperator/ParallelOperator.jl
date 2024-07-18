"""
    Parallel
Methods for constructing the parallel operator
"""
module ParallelOperator

    using DifferentialEquations: ODEProblem, EnsembleProblem, EnsembleSerial, solve, Tsit5, remake
    # using Interpolations: LinearInterpolation
    using Interpolations
    using BasicInterpolators: BicubicInterpolator, LinearInterpolator

    using LinearAlgebra: norm

    using FaADE.Helpers
    using FaADE.Derivatives: build_H
    using FaADE.Derivatives: MassMatrix, DiagonalH, CompositeH

    using FaADE.Grid

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