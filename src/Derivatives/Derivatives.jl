"""
    Derivatives
Module containing first and second derivative variable coefficient SBP operators.
"""
module Derivatives

    using Distributed
    using SBP_operators.Helpers

    include("DerivativeFirst.jl")
    include("DerivativeSecond.jl")

    include("Op_Dx.jl")
    include("Op_Dxx.jl")

    include("Diff.jl")

    export SecondDerivative
    export Dₓ, Dₓ!, Dₓₓ, Dₓₓ!, generate_SecondDerivative

end