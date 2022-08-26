module types

    struct NodeType{T} end
    const Left = NodeType{:Left}()
    const Internal = NodeType{:Internal}()
    const Right = NodeType{:Right}()

    struct BoundaryCondition{T} end
    const Dirichlet = BoundaryCondition{:Dirichlet}()
    const Neumann = BoundaryCondition{:Neumann}()
    const Robin = BoundaryCondition{:Robin}()
    const Periodic = BoundaryCondition{:Periodic}()
    const Interface = BoundaryCondition{:Interface}()

    export NodeType, Left, Internal, Right
    export BoundaryCondition, Dirichlet, Neumann, Robin, Periodic, Interface
    
end

