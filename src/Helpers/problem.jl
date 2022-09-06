


abstract type GridType end



struct Grid1D{T} <: GridType
    domain  :: Vector{T}
    Δ       :: T
    n       :: Int
    order   :: Int
end



function Grid1D(𝒟::Vector{T},n::Int,order::Int) where T
    check_order(order)

    Δx = (𝒟[2]-𝒟[1])/(n-1)
    x = collect(range(𝒟[1],𝒟[2],length=n))

    Grid1D{T}(x,Δx,n,order)
end




struct SquareGrid <: GridType
    grids   :: Tuple{Grid1D}
    axis    :: Vector{Int}
end


function SquareGrid(gridax::Tuple{Grid1D,Int}...)
    grids = Tuple([g[1] for g in gridax])
    axis = [ax[2] for ax in gridax]
    SquareGrid(grids,axis)
end

# CartesianIndices((x_internal.indices[1],y_internal.indices[1]))
# y_internal = CartesianIndices(y)
# x_internal = CartesianIndices(x)

#=
struct PDE
    PDE     :: Function
    grid    :: NamedTuple{grid}
    order   :: Int
end


struct VariableCoefficientPDE
    grid    :: Vector{grid}
    K       :: Vector{AbstractArray}
    order   :: Vector{Int}
end
=#




