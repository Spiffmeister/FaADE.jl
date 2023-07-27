

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

Base.show(M::CartesianMetric) = print("Cartesian Metric")
Base.show(M::CurvilinearMetric) = print("Curvilinear Metric")


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
        NG,
        MET,
        TG,
        TJ,
        IT <: AbstractArray{Int}} <: GridType{TT,DIM,MET}
    
    Grids   :: TG
    Joint   :: TJ
    inds    :: IT

    function GridMultiBlock(grids::Vector{Grid1D{TT,MET}},joints) where {TT,MET}

        inds = [sum([grids[j].n for j in 1:i]) for i in 1:length(grids)]

        new{TT, 1, length(grids), MET, typeof(grids), typeof(joints),typeof(inds)}(grids,joints,inds)
    end
    function GridMultiBlock(grids::Vector{Grid2D{TT,MET}},joints) where {TT,MET}
        
        inds = [sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]

        new{TT,2, length(grids),MET,typeof(grids),typeof(joints),typeof(inds)}(grids,joints,[1])
    end
end
"""
    GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)]
    GridMultiBlock(grids,J)
end
"""
    GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}},J::Vector{Tuple{Int,Int,NodeType}}) where {TT,MET}

Assumes blocks are stacked in `x`
"""
function GridMultiBlock(grids::AbstractArray{Grid2D{TT,MET}}) where {TT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)]
    GridMultiBlock(grids,J)
end

"""
    GetMinΔ
Return the miniumum grid size between ``\\Delta x`` and ``\\Delta y``
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)







function Base.getindex(G::GridMultiBlock{TT,1},i::Integer) where TT
    ii = findfirst(x->x ≥ i, G.inds)
    ii == 1 ? iii = i : iii = i - G.inds[ii-1]
    return G.Grids[ii].grid[iii]
end
function Base.getindex(G::GridMultiBlock{TT,2},i::Integer,j::Integer) where TT
    ii = findfirst(x->x ≥ i, G.inds[1,:])
    ii == 1 ? iii = i : iii = i - G.inds[1,ii-1]
    
    jj = findfirst(x->x ≥ j, G.inds)
    jj == 1 ? jjj = j : jjj = j - G.inds[jj-1]
    return G.Grids[ii,jj].grid[iii,jjj]
end



"""
    size(G::GridType)
"""
Base.size(G::Grid1D) = (G.n,)
Base.size(G::Grid2D) = (G.nx,G.ny)
Base.size(G::GridMultiBlock{TT,1}) where {TT} = (G.inds[end]-length(G.inds)+1,)
Base.size(G::GridMultiBlock{TT,2}) where {TT} = (G.inds[end,1],G.inds[end,2])
"""
    length(G::GridType)
"""
Base.length(G::GridType) = prod(size(G))

Base.ndims(G::GridType{TT,DIM,AT}) where {TT,DIM,AT} = DIM

Base.eachindex(G::GridMultiBlock{TT,1}) where {TT} = Base.OneTo(length(G))


# Base.typeof(M::MetricType{MType}) where MType = MType

# Base.getindex(G::GridMultiBlock{},i)

