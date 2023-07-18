

# abstract type GridType{DIM,COORD,dtype<:Real} end
abstract type GridType{dtype<:AbstractFloat,DIM,COORD} end
abstract type LocalGridType{dtype,DIM,COORD} <: GridType{dtype,DIM,COORD} end

abstract type MetricType{MType} end




struct CartesianMetric <: MetricType{:cartesian} end
struct CurvilinearMetric{
        F1<:Function,
        F2<:Function,
        F3<:Function,
        F4<:Function} <: MetricType{:curvilinear}
    dxdr    ::  Function
    dxdq    ::  Function
    dydr    ::  Function
    dydq    ::  Function
end

"""
    Grid1D
Grid data structure for 1 dimensional problems.
    
    `Grid1D{T}(𝒟::Vector,n::Integer)`
Inputs:
- Vector of domain boundaries `[x_0,x_N]`
- Number of grid points.

Returns:
- Struct for 1D grid object containing vector of grid points, ``\\Delta x`` and ``n``.
"""
struct Grid1D{TT,MET} <: LocalGridType{TT,1,MET}
    grid    :: Vector{TT}
    Δx      :: TT
    n       :: Integer
    M       :: MET
    function Grid1D(𝒟::Vector{TT},n::Integer,MType) where TT
        Δx = (𝒟[2]-𝒟[1])/(n-1)
        x = collect(range(𝒟[1],𝒟[2],length=n))

        new{TT,typeof(MType)}(x,Δx,n,MType)
    end
end
Grid1D(𝒟,n) = Grid1D(𝒟,n,CartesianMetric())
"""
    Grid2D
Grid data structure for 2 dimensional problems.

    `Grid2D{T}(𝒟x::Vector,𝒟y::Vector,nx::Integer,ny::Integer)`
Inputs:
- Domain boundaries in ``x``
- Domain boundaries in ``y``
- Number of nodes in ``x``
- Number of nodes in ``y``

Returns:
- Struct for 2D grid object containing grid points in ``x`` and ``y``, ``\\Delta x`` and ``\\Delta y``, and ``n_x`` and ``n_y``.
"""
struct Grid2D{TT,MET} <: LocalGridType{TT,2,MET}
    gridx   :: Vector{TT}
    gridy   :: Vector{TT}
    Δx      :: TT
    Δy      :: TT
    nx      :: Integer
    ny      :: Integer
    M       :: MET
    function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer,MType) where TT
        gx = Grid1D(𝒟x,nx)
        gy = Grid1D(𝒟y,ny)

        new{TT,typeof(MType)}(gx.grid,gy.grid, gx.Δx,gy.Δx, gx.n,gy.n,MType)
    end
end
Grid2D(𝒟x,𝒟y,nx,ny) = Grid2D(𝒟x,𝒟y,nx,ny,CartesianMetric())



"""
    GridMultiBlock
Grid data for SAT boundary problems

```
D1 = Grid1D([0.0,0.5],11)
D2 = Grid1D([0.5,1.0],6)

FaADE.Helpers.GridMultiBlock([D1,D2])
```

```
D1 = Grid1D([0.0,0.5],11)
D2 = Grid1D([0.5,1.0],6)

FaADE.Helpers.GridMultiBlock
```
"""
struct GridMultiBlock{TT  <: Real,
        DIM,
        ST,
        MET,
        TG,
        TJ } <: GridType{TT,DIM,MET}
    
    Grids           :: TG
    Joint           :: TJ

    function GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}

        # new{TT, ndims(grids[1]), MET, typeof(TG)}(grids)
        J = [(i,i+1,Right) for i in 1:length(grids)]

        new{TT, 1, :structured, MET, typeof(grids), typeof(J)}(grids,J)
    end
end
"""
    GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}}) where {TT,MET}
Assumes the grid layout is stitched together in `x`, starting at the lower bound in `x`
"""
function GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}}) where {TT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)]
    new{TT,2,:structured,MET,typeof(grids),typeof(J)}(grids,J)
end
"""
    GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}},J::Vector{Tuple{Int,Int,NodeType}}) where {TT,MET}

`J` is a vector of tuples which tells the code what the layout of the local grid blocks in `grids`
"""
function GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}},J::Vector{Tuple{Int,Int,NodeType}}) where {TT,MET}
    new{TT,2,:unstructured,MET,typeof(grids),typeof(J)}(grids,J)
end

"""
    GetMinΔ
Return the miniumum grid size between ``\\Delta x`` and ``\\Delta y``
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)


Base.ndims(G::GridType{TT,DIM,MET}) where {TT,DIM,MET} = DIM

# Base.typeof(M::MetricType{MType}) where MType = MType


