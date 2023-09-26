

# abstract type GridType{DIM,COORD,dtype<:Real} end
abstract type GridType{dtype<:AbstractFloat,DIM,COORD} end
abstract type LocalGridType{dtype,DIM,COORD} <: GridType{dtype,DIM,COORD} end

abstract type MetricType{MType} end


struct CartesianMetric <: MetricType{:cartesian} end
struct CurvilinearMetric <: MetricType{:curvilinear} end

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
struct Grid1D{TT<:Real,
        DT<:Union{Real,Vector{TT}},
        MET<:MetricType} <: LocalGridType{TT,1,MET}

    grid    :: Vector{TT}
    Δx      :: DT
    n       :: Int64

end
function Grid1D(𝒟::Vector{TT},n::Integer) where TT
    Δx = (𝒟[2]-𝒟[1])/(n-1)
    x = collect(range(𝒟[1],𝒟[2],length=n))

    # new{TT,TT,CartesianMetric}(x,Δx,n)
    return Grid1D{TT,typeof(Δx),CartesianMetric}(x,Δx,n)
end
function Grid1D(𝒟::Vector{TT}) where TT
    Δx = diff(𝒟)
    # new{TT,Vector{TT},CurvilinearMetric}(𝒟,Δx,length(𝒟))
    return Grid1D{TT,typeof(Δx),CurvilinearMetric}(𝒟,Δx,length(𝒟))
end


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
struct Grid2D{TT,
        DT<:Union{Real,Vector{TT}},
        MET<:MetricType} <: LocalGridType{TT,2,MET}
    gridx   :: Vector{TT}
    gridy   :: Vector{TT}
    Δx      :: DT
    Δy      :: DT
    nx      :: Integer
    ny      :: Integer
end
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer) where TT
    gx = Grid1D(𝒟x,nx)
    gy = Grid1D(𝒟y,ny)

    return Grid2D{TT,typeof(gx.Δx),CartesianMetric}(gx.grid, gy.grid, gx.Δx, gy.Δx, gx.n, gy.n)
end
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT}) where TT
    Δx = diff(𝒟x)
    Δy = diff(𝒟y)
    return Grid2D{TT,typeof(Δx),CurvilinearMetric}(𝒟x,𝒟y,Δx,Δy,length(𝒟x),length(𝒟y))
end



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
        MET,
        TG,
        TJ,
        IT} <: GridType{TT,DIM,MET}
    
    Grids   :: TG
    Joint   :: TJ
    inds    :: IT
end

function GridMultiBlock(grids::Vector{Grid1D{TT,DT,MET}},joints) where {TT,DT,MET}
    inds = [sum([grids[j].n for j in 1:i]) for i in 1:length(grids)]
    return GridMultiBlock{TT, 1, MET, typeof(grids), typeof(joints),typeof(inds)}(grids,joints,inds)
end
function GridMultiBlock(grids::Vector{Grid2D{TT,DT,MET}},joints) where {TT,DT,MET}
    inds = [sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]
    return GridMultiBlock{TT,2, MET,typeof(grids),typeof(joints),typeof(inds)}(grids,joints,inds)
end
"""
    GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::Vector{Grid1D{TT,DT,MET}}) where {TT,DT,MET}
    J = [(i,i+1,Right) for i in 1:length(grids)-1]
    # J = [(1,2,Right)], [(i,i-1,Left),]
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

eachgrid(G::GridMultiBlock) = Base.OneTo(length(G.Grids))


# Base.typeof(M::MetricType{MType}) where MType = MType

# Base.getindex(G::GridMultiBlock{},i)
Base.getindex(G::Grid1D,i...) = G.grid[i...]
Base.getindex(G::Grid2D,i::Integer,j::Integer) = (G.gridx[i],G.gridy[j])



Base.lastindex(G::Grid1D) = G.n
Base.lastindex(G::Grid2D) = size(G)


Base.eltype(G::GridType{TT}) where TT = TT