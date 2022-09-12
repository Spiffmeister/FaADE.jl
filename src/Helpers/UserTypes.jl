


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
    grid                :: GridType
    K                   :: AbstractArray{T}
    order               :: Vector{Int}
    BoundaryConditions  :: Vector{BoundaryConditionData}
    function VariableCoefficientPDE1D(grid,K,order,BCs...) where T
        new{T}(grid,K,order,BCs)
    end
end
struct VariableCoefficientPDE2D{T} <: PDEProblem
    grid                :: GridType
    Kx                  :: AbstractArray{T}
    Ky                  :: AbstractArray{T}
    order               :: Vector{Int}
    BoundaryConditions  :: Vector{BoundaryConditionData}
    function VariableCoefficientPDE2D(grid,K,order,BCs...) where T
        new{T}(grid,K,order,BCs)
    end
end




