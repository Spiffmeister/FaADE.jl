
"""
    NodeType

- Used by `FirstDerivative()` and `SecondDerivative()` to indicate which node type should be used.
- `Left` and `Right` are used by the `SAT()` function to determine which side of the SAT to call.

Node types are:
- Left
- Internal
- Right
- Up
- Down
"""
struct NodeType{T,D} end
const Left = NodeType{:Left,1}()
const Internal = NodeType{:Internal,0}()
const Right = NodeType{:Right,1}()
const Up = NodeType{:Right,2}()
const Down = NodeType{:Left,2}()

function _flipside(NT::NodeType)
    if NT == Left
        return Right
    elseif NT == Right
        return Left
    elseif NT == Up
        return Down
    elseif NT == Down
        return Up
    else
        return NT
    end
end


"""
    BoundaryCondition

Defines the boundary conditions, options are:
- Dirichlet
- Neumann
- Robin
- Periodic
- Interface
"""
struct BoundaryConditionType{T} end
const Dirichlet = BoundaryConditionType{:Dirichlet}()
const Neumann = BoundaryConditionType{:Neumann}()
const Robin = BoundaryConditionType{:Robin}()
const Periodic = BoundaryConditionType{:Periodic}()
const Interface = BoundaryConditionType{:Interface}()



"""
    SATMode

Used when the conjugate gradient solver is being used so the solver knows which part of the SAT to call, since ``F`` in ``Ab=F`` contains all data (the right hand side of ``u(0)=f(t)``, but not the solution ``u``). 
"""
struct SATMode{T} end
const DataMode = SATMode{:DataMode}()
const SolutionMode = SATMode{:SolutionMode}()
const ExplicitMode = SATMode{:Explicit}()


struct DerivativeOrder{InternalOrder,BoundarySize} end


abstract type DataBlockType{dtype<:AbstractFloat,N, atype<:AbstractArray{dtype}} end
# abstract type BoundaryStorage{dtype<:AbstractFloat,N, atype<:AbstractArray{dtype}} end
# abstract type GridType{dtype<:AbstractFloat,N} end
abstract type ParallelGridStorage{dtype<:AbstractFloat,N} end



struct SourceTerm{ST}
    source :: ST
end


Base.show(io::IO,S::SourceTerm{TT}) where TT = S.source

