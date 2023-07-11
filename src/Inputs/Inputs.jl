module Inputs

    using FaADE.Helpers
    using FaADE.SATs
    using FaADE.Derivatives
    using FaADE.ParallelOperator

    include("UserTypes.jl")

    export newProblem1D, newProblem2D

end