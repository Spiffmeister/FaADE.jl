module Helpers

    include("types.jl")
    include("checks.jl")
    include("problem.jl")


    # types import
    export NodeType, Left, Internal, Right
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    export SATMode, DataMode, SolutionMode
    
    # problem import
    export GridType, Grid1D


end

