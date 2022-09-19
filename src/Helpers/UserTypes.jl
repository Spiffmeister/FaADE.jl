


abstract type BoundaryConditionData end
abstract type PDEProblem end


"""
    User defined boundary types
"""
struct Boundary <: BoundaryConditionData
    type    :: BoundaryConditionType
    RHS     :: Function
    side    :: NodeType
    axis    :: Int
    function Boundary(type::BoundaryConditionType,RHS::Function,side::NodeType,axis::Int)
        type ∈ [Dirichlet,Neumann,Robin] ? nothing : error("Input 1 must be Dirichlet, Neumann or Robin")
        side ∈ [Left,Right] ? nothing : error("Input 3 must be Left or Right")
        new(type,RHS,side,axis)
    end
end
struct PeriodicBoundary <: BoundaryConditionData
    type    :: BoundaryConditionType
    axis    :: Int
    function PeriodicBoundary(axis::Int)
        new(Periodic,axis)
    end
end




"""
    PDE Problem type for user input
"""
struct VariableCoefficientPDE1D{T} <: PDEProblem
    InitialCondition    :: Function
    K                   :: AbstractArray{T}
    order               :: Int
    BoundaryConditions  :: NamedTuple
    function VariableCoefficientPDE1D{T}(u₀,K,order,BC) where T
        new(u₀,K,order,BC)
    end
end
struct VariableCoefficientPDE2D{T,D} <: PDEProblem
    InitialCondition    :: Function
    Kx                  :: AbstractArray{T}
    Ky                  :: AbstractArray{T}
    order               :: Vector{Int}
    BoundaryConditions  :: Vector{D}
    function VariableCoefficientPDE2D{T,D}(u₀,Kx::AbstractArray{T},Ky::AbstractArray{T},order::Int,BCs::D...) where {T,D<:BoundaryConditionData}
        new(u₀,Kx,Ky,order,[bc for bc in BCs])
    end
end



function VariableCoefficientPDE1D(u₀::Function,K::AbstractVector{T},order::Int,BCs::D...) where {T,D<:BoundaryConditionData}
    BCs = (Left=BCs[1],Right=BCs[2])
    VariableCoefficientPDE1D{T}(u₀,K,order,BCs)
end






