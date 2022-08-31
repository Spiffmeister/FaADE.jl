module SATs

    using SBP_operators.types
    using StaticArrays


    include("Boundary_Operators.jl")

    include("Dirichlet.jl")
    include("Neumann.jl")
    include("Robin.jl")
    include("Periodic.jl")
    include("Split.jl")

    include("SAT_Interface.jl")

    export SAT, SAT!, SATAdd!
    export SAT_Dirichlet, SAT_Dirichlet!, SAT_Dirichlet_internal!, SATpenalties, BDₓᵀ
    export SAT_Neumann, SAT_Neumann!
    export SAT_Periodic, SAT_Periodic!
    export SAT_Robin
    export Split_domain
    
end

# get_node_coordinates