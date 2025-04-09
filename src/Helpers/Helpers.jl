"""
    Helpers
Module contains many of the data structures for users and internal use, as well as nifty functions. Some functions are included here and commented out in higher level
    packages because they are required in multiple packages. In theory they could be uncommented in the higher level packages so the package can be used by itself.
"""
module Helpers

    include("types.jl")
    include("nifty.jl")

    include("HMatrix.jl")

    export NodeType, Left, Internal, Right, Up, Down
    export _flipside
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    export SourceTerm, DiffusionCoefficient
    export BoundaryStorage

    export check_order, check_boundary,
        halforder, BoundaryNodeInput, BoundaryNodeOutput,
        _SelectLoopDirection, SATNodeOutput


    export GetAxis, #GetDim,
        GetMinΔ

end
