"""
    NodeType

Used in `SecondDerivative` to indicate which node type should be used.
- NodeLeft
- NodeInternal
- NodeRight
"""
struct NodeType{T} end
const NodeLeft = NodeType{:Left}()
const NodeInternal = NodeType{:Internal}()
const NodeRight = NodeType{:Right}()


struct BoundaryCondition{T} end
const Dirichlet = BoundaryCondition{:Dirichlet}()
const Neumann = BoundaryCondition{:Neumann}()
const Robin = BoundaryCondition{:Robin}()
const Periodic = BoundaryCondition{:Periodic}()
const Interface = BoundaryCondition{:Interface}()