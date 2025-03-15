

# abstract type GridType{DIM,COORD,dtype<:Real} end
abstract type GridType{dtype<:AbstractFloat,DIM,COORD} end
abstract type LocalGridType{dtype,DIM,COORD} <: GridType{dtype,DIM,COORD} end

abstract type MetricType{MType} end


struct CartesianMetric <: MetricType{:cartesian} end
struct CurvilinearMetric <: MetricType{:curvilinear} end
struct StretchedMetric <: MetricType{:stretched} end

Base.show(M::CartesianMetric) = print("Cartesian metric")
Base.show(M::CurvilinearMetric) = print("Curvilinear metric")
Base.show(M::StretchedMetric) = print("Stretched grid metric")

"""
    Grid1D{TT<:Real,MET<:MetricType,GT<:Vector{TT},DT<:Union{Real,Vector{TT}}} <: LocalGridType{TT,1,MET}
Grid data structure for 1 dimensional problems.

Constructors:
```
Grid1D(𝒟::Vector{TT},n::Integer) where TT
Grid1D(𝒟::Vector{TT}) where TT
```

Grid points can be accessed as if the grid object is an array;
```julia
G = Grid2D([0.0,0.5],[0.0,1.0],5,5)
G[1]
```

Inputs:
- Vector of domain boundaries `[x_0,x_N]`
- Number of grid points.

Returns:
- Struct for 1D grid object containing vector of grid points, ``\\Delta x`` and ``n``.
"""
struct Grid1D{TT<:Real,
        MET <:MetricType,
        GT  <:Vector{TT},
        DT  <:Union{Real,Vector{TT}}
            } <: LocalGridType{TT,1,MET}

    grid    :: GT
    Δx      :: DT
    n       :: Int64

end
function Grid1D(𝒟::Vector{TT},n::Integer) where TT
    Δx = (𝒟[2]-𝒟[1])/(n-1)
    x = collect(range(𝒟[1],𝒟[2],length=n))
    return Grid1D{TT,CartesianMetric,typeof(x),typeof(Δx)}(x,Δx,n)
end
function Grid1D(𝒟::Vector{TT}) where TT
    Δx = diff(𝒟)
    return Grid1D{TT,CurvilinearMetric,typeof(x),typeof(Δx)}(𝒟,Δx,length(𝒟))
end


"""
    Grid2D{TT,MET<:MetricType,GT<:AbstractArray{TT}} <: LocalGridType{TT,2,MET}
Grid data structure for 2 dimensional problems.

Grid points can be accessed as if the grid object is an array;
```julia
G = Grid2D([0.0,0.5],[0.0,1.0],5,5)
G[1,1]
```

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
        GT<:AbstractArray{TT},
            } <: LocalGridType{TT,2,MET}
    gridx   :: GT
    gridy   :: GT
    Δx      :: TT
    Δy      :: TT
    nx      :: Integer
    ny      :: Integer

    J       :: GT
    qx      :: GT
    qy      :: GT
    rx      :: GT
    ry      :: GT
end
"""
    Grid2D(𝒟x::Vector,𝒟y::Vector,nx::Integer,ny::Integer)
Construct a 2D grid from the domain boundaries in ``x`` and ``y`` and the number of nodes in ``x`` and ``y``.
"""
function Grid2D(𝒟x::Vector{TT},𝒟y::Vector{TT},nx::Integer,ny::Integer;coord=CartesianMetric) where TT
    gx = Grid1D(𝒟x,nx)
    gy = Grid1D(𝒟y,ny)

    X = repeat(gx.grid,1,ny)
    Y = repeat(gy.grid',nx,1)

    if coord == CartesianMetric
        J = qx = qy = rx = ry = ones(eltype(gx.grid),(1,1))
    else
        J = qx = ry = ones(TT,size(X))
        qy = rx = zeros(TT,size(X))
    end

    return Grid2D{TT,coord,typeof(X)}(X, Y, gx.Δx, gy.Δx, gx.n, gy.n,
        J, qx, qy, rx, ry)
end
"""
    Grid2D(𝒟x::Vector,𝒟y::Vector)
Construct a 2D grid from vectors in ``x`` and ``y`` for curvilinear ``x,y``.

Inputs:
- Grid in ``x``
- Grid in ``y``
Optional:
- ``order=nothing``: Order of the derivative, will use highest order not specified
- ``periodicx=false``: Domain is periodic in ``x``
- ``periodicy=false``: Domain is periodic in ``y``
- ``xmap=true``: Map the grid in ``x`` to the computational domain, set false if the grid is already in logical coordinates
- ``ymap=true``: Map the grid in ``y`` to the computational domain, set false if the grid is already in logical coordinates
"""
function Grid2D(𝒟x::Matrix{TT},𝒟y::Matrix{TT};order=nothing,periodicx=false,periodicy=false,xmap=true,ymap=true) where TT
    
    nx, ny = size(𝒟x)
    if xmap
        Δx = TT(1)/TT(nx-1)
    else
        Δx = (𝒟x[end,1]-𝒟x[1,1])/TT(nx-1)
    end
    if ymap
        Δy = TT(1)/TT(ny-1)
    else
        Δy = (𝒟y[1,end]-𝒟y[1,1])/TT(ny-1)
    end
    
    if isnothing(order)
        if (nx ≥ 16) & (ny ≥ 16)
            order = 4
        else
            order = 2
        end
    end

    xq = zeros(eltype(𝒟x),size(𝒟x))
    xr = zeros(eltype(𝒟x),size(𝒟x))
    yq = zeros(eltype(𝒟y),size(𝒟y))
    yr = zeros(eltype(𝒟y),size(𝒟y))

    # Derivatives of x,y wrt q
    if periodicx
        PeriodicD₁!(xq,𝒟x,nx,Δx,order,1)
        PeriodicD₁!(yq,𝒟y,nx,Δx,order,1)
    else
        D₁!(xq,𝒟x,nx,Δx,order,1)
        D₁!(yq,𝒟y,nx,Δx,order,1)
    end
    
    # Derivatives of x,y wrt r
    if periodicy
        PeriodicD₁!(xr,𝒟x,ny,Δy,order,2)
        PeriodicD₁!(yr,𝒟y,ny,Δy,order,2)
    else
        D₁!(xr,𝒟x,ny,Δy,order,2)
        D₁!(yr,𝒟y,ny,Δy,order,2)
    end

    # Jacobian
    J = zeros(eltype(𝒟x),size(𝒟x))
    for i = 1:nx
        for j = 1:ny
            J[i,j] = (xq[i,j]*yr[i,j] - xr[i,j]*yq[i,j])
        end
    end

    # Computational coordinate derivatives
    qx = yr./J # yr -> qx
    qy = -xr./J # xr -> qy
    rx = -yq./J # yr -> rx
    ry = xq./J # xq -> ry
    
    return Grid2D{TT,CurvilinearMetric,typeof(𝒟x)}(𝒟x, 𝒟y, TT(1)/TT(nx-1), TT(1)/TT(ny-1), nx, ny,
        J, qx, qy, rx, ry)
end
"""
    Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer)
Construct a 2D grid from the boundary functions in ``x`` and ``y`` and the number of nodes in ``x`` and ``y``.

Curves ``c`` are parameterised by ``u`` and ``v`` where ``u`` is the coordinate in the ``x`` direction and ``v`` is the coordinate in the ``y`` direction and where ``u`` and ``v`` are in the range ``[0,1]``.
"""
function Grid2D(cbottom::Function,cleft::Function,cright::Function,ctop::Function,nx::Integer,ny::Integer;order=nothing,stretchu=u->u,stretchv=v->v)
    X,Y = meshgrid(cbottom,cleft,cright,ctop,nx,ny,stretchu=stretchu,stretchv=stretchv)
    Grid2D(X,Y;order=order)
end




"""
    Joint
Specifies the joining surface between two grids.
"""
struct Joint
    index   :: Int64
    side    :: NodeType
    # pair    :: NodeType

    # Joint(index::Int64,side::NodeType;pair::NodeType) = new(index,side,pair)
end




"""
    GridMultiBlock{TT,DIM,MET,TG,TJ,IT} <: GridType{TT,DIM,MET}
Grid data for multiblock problems

Grabbing a particular subgrid can be done by `G.Grids[i]` which indexes in the order the grids are given. 
Indexing can be performed by `G[i]` for 1D or `G[i,j]` for 2D multiblock problems.
`GridMultiBlock.Joint` contains the information on how to connect grids. If periodic boundary conditions are being used, do not specify the joint across that boundary.

Example grid creation,
```julia
    D1  = Grid2D([0.0,0.5],[0.0,1.0],5,5)
    D2  = Grid2D([0.5,1.0],[0.0,1.0],5,5)
    D3  = Grid2D([1.0,1.5],[0.0,1.0],5,5)

    glayout = ([(2,Right)],
                [(1,Left),(3,Right)],
                [(2,Left)])

    G = GridMultiBlock((D1,D2,D3),glayout)
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
    ngrids  :: Int
end
"""
    GridMultiBlock(grids::LocalGridType{TT,1,MET}...) where {TT,MET}
Multiblock grid for 1D grids, assumes the grids are stacked one after the other from left to right
"""
function GridMultiBlock(grids::LocalGridType{TT,1,MET}...) where {TT,MET}
    if length(grids) == 2
        J = ((Joint(2,Right),),(Joint(1,Left),))
    elseif length(grids) > 2
        J = [(Joint(i-1,Left),Joint(i+1,Right)) for i = 2:length(grids)-1]
        J = tuple((Joint(2,Right),),J...,(Joint(length(grids)-1,Left),))
    end
    inds = [sum([grids[j].n for j in 1:i]) for i in 1:length(grids)]
    return GridMultiBlock{TT,1,MET,typeof(grids),typeof(J),typeof(inds)}(grids,J,inds,length(inds))
end
"""
    GridMultiBlock(grids::Tuple{Vararg{Grid2D{TT,MET},N}}) where {N,TT,MET}
Multiblock grid for 2D grids

For a stacked grid (i.e. a series of grids in a row) one can run,
```julia
    D1  = Grid2D([0.0,0.5],[0.0,1.0],5,5)
    D2  = Grid2D([0.5,1.0],[0.0,1.0],5,5)
    D3  = Grid2D([1.0,1.5],[0.0,1.0],5,5)

    glayout = (((2,Right),),
                ((1,Left),(3,Right)),
                ((2,Left),))

    G = GridMultiBlock((D1,D2,D3),glayout)
```

TODO: Add example for non-stacked grids
TODO: Add checking that there are no hanging nodes
"""
function GridMultiBlock(grids::Tuple{Vararg{Grid2D{TT,MET,GT},N}},joints) where {N,TT,GT,MET}
    # inds = [sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]
    # inds = [sum([grids[j].nx])]
    indx = vcat([sum([grids[j].nx] for j in 1:i) for i in 1:length(grids)]...)
    indy = vcat([sum([grids[j].ny] for j in 1:i) for i in 1:length(grids)]...)
    # _checkjoints(joints)
    inds = [indx,indy]
    return GridMultiBlock{TT,2, MET,typeof(grids),typeof(joints),typeof(inds)}(grids,joints,inds,length(inds))
end

#============ Functions ============#

"""
    GetMinΔ
Return the miniumum grid size between ``\\Delta x`` and ``\\Delta y``
"""
function GetMinΔ end
GetMinΔ(grid::Grid1D) = grid.Δx
GetMinΔ(grid::Grid2D) = min(grid.Δx,grid.Δy)

"""
    GetMetricType
"""
GetMetricType(G::GridType{TT,DIM,COORD}) where {TT,DIM,COORD} = COORD





function GetBoundaryCoordinates(grid::Grid2D,side::NodeType)
    if side == Left
        return [(grid.gridx[1,i],grid.gridy[1,i]) for i in 1:grid.ny]
    elseif side == Right
        return [(grid.gridx[end,i],grid.gridy[end,i]) for i in 1:grid.ny]
    elseif side == Down
        return [(grid.gridx[i,1],grid.gridy[i,1]) for i in 1:grid.nx]
    elseif side == Up
        return [(grid.gridx[i,end],grid.gridy[i,end]) for i in 1:grid.nx]
    end
end


function _checkjoints(G::GridMultiBlock)
    joints = G.Joint
    warnprint = false
    for (I,BLOCK) in enumerate(joints)
        for JOINT in BLOCK
            for joint in joints[JOINT.index]
                if (joint.index == I) & (joint.side != Helpers._flipside(JOINT.side))
                    @warn "Blocks $I and $(JOINT.index) have side mismatch. Block $I has block $(JOINT.index) as $(JOINT.side) but block $(JOINT.index) has block $I as $(joint.side)"
                    warnprint = true
                end
            end
        end
    end
    if warnprint
        @warn "The above issues may be caused by certain grid layouts, solver will attempt to work out the correct normal vectors (typically normals defined by blocks with higher index take higher precidence)."
    end
end


"""
Used to index the side of a matrix for multi-block problems using a `joint` from the grid.
"""
function Base.getindex(x::AT,::NodeType{SIDE,AX}) where {AT<:AbstractArray,SIDE,AX}
    SIDE == :Left ? I = 1 : I = size(x)[AX]
    if AX == 1
		return x[I,:]
	elseif AX == 2
		return x[:,I]
	else
		error("Does not correspond to a boundary.")
	end
end



"""
    Base.getindex(G::GridType,i::Integer)
"""
Base.getindex(G::Grid1D,i::Integer) = G.grid[i]
Base.getindex(G::Grid2D,i::Integer) = (G.gridx[i],G.gridy[i])
Base.getindex(G::Grid2D,i::Integer,j::Integer) = (G.gridx[i,j],G.gridy[i,j])



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

"""
    size(G::GridMultiBlock{TT,1})
"""
function Base.size(G::GridMultiBlock{TT,1}) where {TT}
    sz = (0,)
    for i = 1:G.ngrids
        sz = sz .+ size(G.Grids[i])
    end
    return sz
end

function Base.size(G::GridMultiBlock{TT,2}) where {TT}
    sz = (0,0)
    for i = 1:G.ngrids
        sz = sz .+ size(G.Grids[i])
    end
    return sz
end


"""
    Base.length(G::GridType)
"""
Base.length(G::GridType) = prod(size(G)) :: Int

"""
    Base.ndims(G::GridType{TT,DIM,AT}) where {TT,DIM,AT}
"""
Base.ndims(G::GridType{TT,DIM,AT}) where {TT,DIM,AT} = DIM

"""
    Base.eachindex(G::GridType)
"""
# Base.eachindex(G::GridType) = Base.OneTo(length(G))
Base.eachindex(G::GridType) = eachindex(1:length(G))

"""
    eachgrid(G::GridMultiBlock)
"""
eachgrid(G::GridMultiBlock) = Base.OneTo(length(G.Grids))

"""
    eachjoint(G::GridMultiBlock)
"""
eachjoint(G::GridMultiBlock) = Base.OneTo(length(G.Joint))

"""
    Base.lastindex(G::GridType)
"""
Base.lastindex(G::Grid1D) = G.n
Base.lastindex(G::Grid2D) = size(G)

"""
    Base.lastindex(G::GridMultiBlock)
"""
Base.eltype(G::GridType{TT}) where TT = TT


"""
    coordtype(G::GridType{TT,DIM,MET}) where {TT,DIM,MET}
"""
coordtype(G::GridType{TT,DIM,MET}) where {TT,DIM,MET} = MET