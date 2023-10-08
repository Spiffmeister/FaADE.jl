

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
        MET <:MetricType,
        DT  <:Union{Real,Vector{TT}},
        GT  <:Vector{TT}} <: LocalGridType{TT,1,MET}

    grid    :: GT
    Δx      :: DT
    n       :: Int64

end
function Grid1D(𝒟::Vector{TT},n::Integer) where TT
    Δx = (𝒟[2]-𝒟[1])/(n-1)
    x = collect(range(𝒟[1],𝒟[2],length=n))

    # new{TT,TT,CartesianMetric}(x,Δx,n)
    return Grid1D{TT,CartesianMetric,TT,typeof(x)}(x,Δx,n)
end
function Grid1D(𝒟::Vector{TT}) where TT
    Δx = diff(𝒟)
    # new{TT,Vector{TT},CurvilinearMetric}(𝒟,Δx,length(𝒟))
    return Grid1D{TT,CurvilinearMetric,typeof(Δx),typeof(x)}(𝒟,Δx,length(𝒟))
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
        MET<:MetricType,
        DT<:Union{Real,Vector{TT}},
        GT<:AbstractArray{TT}
            } <: LocalGridType{TT,2,MET}
    gridx   :: GT
    gridy   :: GT
    Δx      :: DT
    Δy      :: DT
    nx      :: Integer
    ny      :: Integer
end
"""
    Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer)
2D cartesian grid with domain boundaries ``𝒟x`` and ``𝒟y`` and ``nx`` and ``ny`` nodes in ``x`` and ``y`` respectively.
"""
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer) where TT
    gx = Grid1D(𝒟x,nx)
    gy = Grid1D(𝒟y,ny)
    return Grid2D{TT,CartesianMetric,typeof(gx.grid),typeof(gx.Δx)}(gx.grid, gy.grid, gx.Δx, gy.Δx, gx.n, gy.n)
end
"""
    Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer) where TT
2D curvilinear grid using transfinite interpolation between the four boundary functions.

See [`meshgrid`](@ref) for more details.
"""
function Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer) where TT
    gridx,gridy = meshgrid(cbottom,cleft,cright,ctop,nx,ny)
    Δx = diff(gridx)
    Δy = diff(gridy)
    return Grid2D{TT,CurvilinearMetric,typeof(gridx),typeof(Δx)}(gridx,gridy,Δx,Δy,nx,ny)
end



"""
    GridMultiBlock
Grid data for SAT boundary problems
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

function GridMultiBlock(grids::Vector{Grid1D{TT,MET,DT}},joints) where {TT,MET,DT}
    inds = [sum([grids[j].n for j in 1:i]) for i in 1:length(grids)]
    return GridMultiBlock{TT, 1, MET, typeof(grids), typeof(joints),typeof(inds)}(grids,joints,inds)
end
function GridMultiBlock(grids::Vector{Grid2D{TT,MET,DT}},joints) where {TT,MET,DT}
    inds = [sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]
    return GridMultiBlock{TT,2, MET,typeof(grids),typeof(joints),typeof(inds)}(grids,joints,inds)
end
"""
    GridMultiBlock(grids::Vector{Grid1D{TT,MET}}) where {TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::Vector{Grid1D{TT,MET,DT}}) where {TT,MET,DT}
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