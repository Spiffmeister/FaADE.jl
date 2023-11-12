"""
    Derivatives
Module containing first and second derivative variable coefficient SBP operators.
"""
module Derivatives

    using Distributed
    using FaADE.Helpers

    include("types.jl")

    include("DerivativeFirst.jl")
    include("DerivativeSecond.jl")

    include("Op_Dx.jl")
    include("Op_Dxx.jl")

    include("HMatrix.jl")

    include("Diff.jl")

    # export SecondDerivative
    export DerivativeOrder, GetOrder, DerivativeOperator
    export D₁, D₂
    export D₁!, D₂!
    export mul!

    export build_H, innerH

end