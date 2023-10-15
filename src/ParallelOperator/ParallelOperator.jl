"""
    Parallel
Methods for constructing the parallel operator
"""
module ParallelOperator

    using DifferentialEquations: ODEProblem, EnsembleProblem, EnsembleSerial, solve, Tsit5, remake
    # using Interpolations: LinearInterpolation
    using Interpolations

    using FaADE.Helpers

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