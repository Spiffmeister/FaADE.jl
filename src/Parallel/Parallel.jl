"""
    Parallel
Methods for constructing the parallel operator
"""
module Parallel

    using DifferentialEquations

    using SBP_operators.Helpers


    include("PGrid.jl")

    export ParallelGrid

end