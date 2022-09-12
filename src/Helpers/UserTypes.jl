


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
struct VariableCoefficientPDE1D{T,D} <: PDEProblem
    grid                :: GridType{T}
    K                   :: AbstractArray{T}
    order               :: Int
    BoundaryConditions  :: Vector{D}
    function VariableCoefficientPDE1D{T,D}(grid,K,order,BC) where {T,D}
        new(grid,K,order,BC)
    end
end
struct VariableCoefficientPDE2D{T} <: PDEProblem
    grid                :: GridType{T}
    Kx                  :: AbstractArray{T}
    Ky                  :: AbstractArray{T}
    order               :: Vector{Int}
    BoundaryConditions  :: Vector{BoundaryConditionData}
    function VariableCoefficientPDE2D{T,D}(grid::GridType{T},Kx::AbstractArray{T},Ky::AbstractArray{T},order::Int,BCs::D...) where {T,D}
        new(grid,K,order,[bc for bc in BCs])
    end
end



function VariableCoefficientPDE1D(grid::GridType{T},K::AbstractVector{T},order::Int,BCs::D...) where {T,D<:BoundaryConditionData}
    VariableCoefficientPDE1D{T,D}(grid,K,order,[bc for bc in BCs])
end

