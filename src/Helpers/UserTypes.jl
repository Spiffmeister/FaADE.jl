

"""
    BoundaryConditionData
Abstract type for boundary condition input structs
"""
abstract type BoundaryConditionData end
"""
    PDEProblem
Abstract type for the problem struct
"""
abstract type PDEProblem end


"""
    Boundary{BT<:BoundaryConditionType,BF<:Function,BN<:NodeType,BA<:Integer}
User defined boundary conditions for Dirichlet or Neumann.
TODO: Robin boundaries

Constructors:
```julia
Boundary(type::BoundaryConditionType,RHS::Function,side::NodeType,axis::Int)
Boundary(type::BoundaryConditionType,RHS::Function,side::NodeType)
```

Inputs:
- `Dirichlet` or `Neumann` (see [`BoundaryConditionType`](@ref)).
- `f(x,t)` or `f(x,y,t)` depending on the function which gives the boundary condition.
- `Left`, `Right`, `Up`, or `Down` (see [`NodeType`](@ref))
- `1` if along the ``x`` axis or `2` along the ``y`` axis. _TODO: Remove this requirement._ 

Returns:
- Struct required for `SAT` construction.
"""
struct Boundary{BT<:BoundaryConditionType,BF<:Function,BN<:NodeType,BA<:Integer} <: BoundaryConditionData
    type    :: BT
    RHS     :: BF
    side    :: BN
    axis    :: BA
    function Boundary(BCtype::BoundaryConditionType,RHS::Function,side::NodeType,axis::Int)
        BCtype ∈ [Dirichlet,Neumann,Robin] ? nothing : error("Input 1 must be Dirichlet, Neumann or Robin")
        new{typeof(BCtype),typeof(RHS),NodeType,typeof(axis)}(BCtype,RHS,side,axis)
    end
end
function Boundary(BCtype::BoundaryConditionType,RHS::Function,side::NodeType)
    BCtype ∈ [Dirichlet,Neumann,Robin] ? nothing : error("Input 1 must be Dirichlet, Neumann or Robin")
    Boundary(BCtype,RHS,side,GetAxis(side))
end
"""
    PeriodicBoundary(axis::Int)
User defined periodic boundary conditions for Periodic boundaries.

Inputs: 
- `1` if periodic in ``x`` and `2` if periodic in ``y``.

Returns:
- Struct required for `SAT` construction.
"""
struct PeriodicBoundary <: BoundaryConditionData
    type    :: BoundaryConditionType
    axis    :: Int
    function PeriodicBoundary(axis::Int)
        new(Periodic,axis)
    end
end

struct BoundaryConditions{L<:BoundaryConditionData,R<:BoundaryConditionData,U<:BoundaryConditionData,D<:BoundaryConditionData}
    Left    :: L
    Right   :: R
    Up      :: U
    Down    :: D
    function BoundaryConditions(BC::NamedTuple)
        new{typeof(BC.Left),typeof(BC.Right),typeof(BC.Up),typeof(BC.Down)}(BC.Left,BC.Right,BC.Up,BC.Down)
    end
end


"""
    VariableCoefficientPDE1D(u₀::Function,K::Function,order::Int,BCs::BoundaryConditionData...)
PDE Problem type for user input

Inputs:
- Initial condition
- Diffusion coefficient function
- Solution order
- Boundary conditions given by [`Boundary`](@ref) or [`PeriodicBoundary`](@ref)

Returns
- Struct of data required for `solver`
"""
struct VariableCoefficientPDE1D <: PDEProblem
    InitialCondition    :: Function
    K                   :: Function
    order               :: Int
    BoundaryConditions  :: NamedTuple
    function VariableCoefficientPDE1D(u₀::Function,K::Function,order::Int,BCs::BoundaryConditionData...)
        Bounds = NamedTuple()
        for BC in BCs
            Bounds = merge(Bounds,BC)
        end
        new(u₀,K,order,Bounds)
    end
end
"""
    VariableCoefficientPDE2D
PDE Problem type for user input

Inputs:
- Initial condition
- Diffusion coefficient function in ``x``
- Diffusion coefficient function in ``y``
- Solution order
- Boundary conditions given by [`Boundary`](@ref) or [`PeriodicBoundary`](@ref)
    
Returns
- Struct of data required for `solver`
"""
struct VariableCoefficientPDE2D{FC,FX,FY} <: PDEProblem
    InitialCondition    :: FC
    Kx                  :: FX
    Ky                  :: FY
    order               :: Int
    BoundaryConditions  :: NamedTuple
    function VariableCoefficientPDE2D(u₀,Kx::Function,Ky::Function,order::Int,BCs::BoundaryConditionData...)

        Bounds = NamedTuple()
        for BC in BCs
            Bounds = merge(Bounds,BC)
        end

        # new(u₀,Kx,Ky,order,BoundaryConditions(Bounds))
        new{typeof(u₀),typeof(Kx),typeof(Ky)}(u₀,Kx,Ky,order,Bounds)
    end
end



"""
    merge
Extension of Base.merge to add user boundary conditions to the PDEProblem
"""
function Base.merge(Bound::NamedTuple,BC::PeriodicBoundary)
    if BC.axis == 1
        Bound = merge(Bound, (Left=BC,Right=BC))
    elseif BC.axis == 2
        Bound = merge(Bound, (Up=BC,Down=BC))
    end
    return Bound
end
function Base.merge(Bound::NamedTuple,BC::Boundary)
    if BC.side == Left
        Bound = merge(Bound, (Left=BC,))
    elseif BC.side == Right
        Bound = merge(Bound, (Right=BC,))
    elseif BC.side == Up
        Bound = merge(Bound, (Up=BC,))
    elseif BC.side == Down
        Bound = merge(Bound, (Down=BC,))
    end
    return Bound
end



    

# Base.show(io::IO, B::Boundary) = print(io, )



