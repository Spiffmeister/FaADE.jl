"""
    SATs

Submodule containing the simulatenous approximation term constructor and functions.
"""
module SATs


    # import Base: +
    using LinearAlgebra: dot
    using SPADE.Helpers
    # using StaticArrays


    include("SAT_Interface.jl")

    include("Boundary_Operators.jl")

    include("Dirichlet.jl")
    include("Neumann.jl")
    # include("Robin.jl")
    include("Periodic.jl")
    # include("Split.jl")


    export SAT_Dirichlet, SAT_Dirichlet_implicit_data!, SAT_Dirichlet_implicit!, SAT_Dirichlet_explicit!
    export SAT_Neumann
    export SAT_Periodic, SAT_Periodic!
    # export SAT_Robin
    # export Split_domain


    export SimultanousApproximationTerm
    export construct_SAT
    
end
