"""
    Parallel
Methods for constructing the parallel operator
"""
module ParallelOperator

    using DifferentialEquations: ODEProblem, EnsembleProblem, EnsembleSerial, solve, Tsit5, remake
    # using Interpolations: LinearInterpolation
    using Interpolations
    using LinearAlgebra: norm

    using FaADE.Helpers
    using FaADE.Derivatives: build_H

    using FaADE.Grid

    include("ParallelData.jl")
    include("ConstructParallelGrid.jl")
    include("penaltyfn.jl")

    export ParallelGrid
    export ParallelData

    export construct_grid
    export generate_parallel_penalty
    export applyParallelPenalty!

end