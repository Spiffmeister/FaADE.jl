"""
    Inputs
Module for interfacing with `FaADE.jl` and creating problem objects to run the code on.
"""
module Inputs

    using FaADE.Helpers
    using FaADE.Grid
    using FaADE.SATs
    using FaADE.Derivatives
    using FaADE.ParallelOperator

    include("parsing.jl")
    include("UserTypes.jl")
    

    export PDEProblem
    export Problem1D, Problem2D



end