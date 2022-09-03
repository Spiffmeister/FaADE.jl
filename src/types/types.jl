module types


    include("Boundaries.jl")


    struct NodeType{T} end
    const Left = NodeType{:Left}()
    const Internal = NodeType{:Internal}()
    const Right = NodeType{:Right}()


    struct SATMode{T} end
    const DataMode = SATMode{:DataMode}()
    const SolutionMode = SATMode{:SolutionMode}()


    export NodeType, Left, Internal, Right
    export BoundaryConditionType, Dirichlet, Neumann, Robin, Periodic, Interface
    # export Boundary, PeriodicBoundary
    export SATMode, DataMode, SolutionMode

    
end

