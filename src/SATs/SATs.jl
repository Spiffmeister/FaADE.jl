module SATs


    # import Base: +

    using SBP_operators.types
    using StaticArrays

    include("SAT_Interface.jl")

    include("Boundary_Operators.jl")

    include("Dirichlet.jl")
    include("Neumann.jl")
    include("Robin.jl")
    include("Periodic.jl")
    include("Split.jl")


    export SAT, SAT!, SATAdd!
    export SAT_Dirichlet, SAT_Dirichlet!, SAT_Dirichlet_internal!, SATpenalties, BDₓᵀ, SAT_Dirichlet_internal_forcing!
    export SAT_Neumann, SAT_Neumann!
    export SAT_Periodic, SAT_Periodic!
    export SAT_Robin
    export Split_domain
    
end

# get_node_coordinates