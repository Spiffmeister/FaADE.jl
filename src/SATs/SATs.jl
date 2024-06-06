"""
    SATs

Submodule containing the simulatenous approximation term constructor and functions.
"""
module SATs


    # import Base: +
    using LinearAlgebra: dot
    using FaADE.Helpers
    using FaADE.Derivatives: FirstDerivativeTranspose!, D₁!
    # using StaticArrays

    include("types.jl")

    # include("SAT_Interface.jl")

    include("Boundary_Operators.jl")

    include("Dirichlet.jl")
    include("Neumann.jl")
    # include("Robin.jl")
    include("Periodic.jl")
    include("Split.jl")
    # include("Helpers.jl")

    export SATMode, DataMode, SolutionMode, ExplicitMode
    export SAT_Dirichlet, SAT_Dirichlet_explicit!, SAT_Dirichlet_solution!, SAT_Dirichlet_data!
    export SAT_Neumann, SAT_Neumann_data!, SAT_Neumann_solution!
    export SAT_Periodic, SAT_Periodic!
    # export SAT_Robin
    # export SAT_Interface, SAT_Interface!


    export SimultanousApproximationTerm
    # export construct_SAT
    # export SATBoundaries
    # export applySAT!
    
end
