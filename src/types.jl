"""
    NodeType

- Used by `FirstDerivative()` and `SecondDerivative()` to indicate which node type should be used.
- `Left` and `Right` are used by the `SAT()` function to determine which side of the SAT to call.

Node types are:
- Left
- Internal
- Right
"""
struct NodeType{T} end
const Left = NodeType{:Left}()
const Internal = NodeType{:Internal}()
const Right = NodeType{:Right}()

"""
    BoundaryCondition

Used in `SAT()` function to determine which boundary conditions to call.
- Dirichlet
- Neumann
- Robin

- Periodic
- Interface
"""
struct BoundaryCondition{T} end
const Dirichlet = BoundaryCondition{:Dirichlet}()
const Neumann = BoundaryCondition{:Neumann}()
const Robin = BoundaryCondition{:Robin}()
const Periodic = BoundaryCondition{:Periodic}()
const Interface = BoundaryCondition{:Interface}()