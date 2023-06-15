"""
    Parallel
Methods for constructing the parallel operator
"""
module Parallel

    using DifferentialEquations: ODEProblem, EnsembleProblem, EnsembleSerial, solve, Tsit5, remake
    # using Interpolations: LinearInterpolation
    using Interpolations

    using FaADE.Helpers


    include("PGrid.jl")
    include("penaltyfn.jl")

    export ParallelGrid

    export construct_grid
    export generate_parallel_penalty

end