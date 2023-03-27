


abstract type BoundaryConditionData end
abstract type PDEProblem end


"""
    Boundary
User defined boundary conditions.
"""
struct Boundary <: BoundaryConditionData
    type    :: BoundaryConditionType
    RHS     :: Function
    side    :: NodeType
    axis    :: Int
    function Boundary(type::BoundaryConditionType,RHS::Function,side::NodeType,axis::Int)
        type ∈ [Dirichlet,Neumann,Robin] ? nothing : error("Input 1 must be Dirichlet, Neumann or Robin")
        new(type,RHS,side,axis)
    end
end
"""
    PeriodicBoundary
User defined periodic boundary conditions.
"""
struct PeriodicBoundary <: BoundaryConditionData
    type    :: BoundaryConditionType
    axis    :: Int
    function PeriodicBoundary(axis::Int)
        new(Periodic,axis)
    end
end




"""
    VariableCoefficientPDE1D
PDE Problem type for user input
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
"""
struct VariableCoefficientPDE2D <: PDEProblem
    InitialCondition    :: Function
    Kx                  :: Function
    Ky                  :: Function
    order               :: Int
    BoundaryConditions  :: NamedTuple
    function VariableCoefficientPDE2D(u₀,Kx::Function,Ky::Function,order::Int,BCs::BoundaryConditionData...)

        Bounds = NamedTuple()
        for BC in BCs
            Bounds = merge(Bounds,BC)
        end
        new(u₀,Kx,Ky,order,Bounds)
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



    




