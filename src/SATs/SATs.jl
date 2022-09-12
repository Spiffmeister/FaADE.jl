module SATs


    # import Base: +

    using SBP_operators.Helpers
    # using StaticArrays


    include("SAT_Interface.jl")

    include("Boundary_Operators.jl")

    include("Dirichlet.jl")
    # include("Neumann.jl")
    # include("Robin.jl")
    include("Periodic.jl")
    # include("Split.jl")


    export SAT_Dirichlet_implicit_data!, SAT_Dirichlet_implicit!
    # export SAT_Neumann, SAT_Neumann!
    export SAT_Periodic, SAT_Periodic!
    # export SAT_Robin
    # export Split_domain


    export SimultanousApproximationTerm
    export SATDirichlet
    export construct_SAT
    
end
