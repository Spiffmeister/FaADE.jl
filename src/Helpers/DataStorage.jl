

# abstract type Solution{T<:AbstractFloat} end


#========== WHOLE PROBLEM DATA ==========#
"""
    DataBlock{T}
Passed around internally between functions. Only contains data required for current timestep.
"""
mutable struct DataBlock{T,N} <: DataBlockType{T,N}
    dim         :: Int
    grid        :: GridType
    u           :: AbstractArray{T}
    uₙ₊₁        :: AbstractArray{T}
    K           :: Union{Vector,Vector{AbstractArray{T}}}
    boundary    :: BoundaryStorage
    t           :: T
    Δt          :: T
    function DataBlock{T}(
            boundaries::NamedTuple,
            grid::GridType,
            Δt::T,
            order::Int,
            K::AbstractArray{T}...) where {T}
    
        # If grid is 1D or 2D construct the right DataBlock
        if typeof(grid) <: Grid1D
            u   = zeros(T,grid.n)
            uₓₓ = zeros(T,grid.n)
            BStor = BoundaryData1D{T}(boundaries,order)
            DiffCoeff = K[1]
            dim = 1
        elseif typeof(grid) <: Grid2D
            u   = zeros(T,(grid.nx,grid.ny))
            uₓₓ = zeros(T,(grid.nx,grid.ny))
            BStor = BoundaryData2D{T}(boundaries,grid,order)
            DiffCoeff = [K[1],K[2]]
            dim = 2
        end
        new{T,dim}(dim,grid,u,uₓₓ,DiffCoeff,BStor,0,Δt)
    end
end


#========== WHOLE PROBLEM DATA ==========#


#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D{T}
Data structure for storage of SATs in 1 dimensional problems
"""
mutable struct BoundaryData1D{T} <: BoundaryStorage{T,1}
    Type_Left   :: BoundaryConditionType
    Type_Right  :: BoundaryConditionType

    SAT_Left    :: AbstractArray{T}
    SAT_Right   :: AbstractArray{T}

    u_Left      :: AbstractArray{T}
    u_Right     :: AbstractArray{T}

    RHS_Left    :: AbstractArray{T}
    RHS_Right   :: AbstractArray{T}

    function BoundaryData1D{T}(BC::NamedTuple,order::Int) where {T}

        nnodes = SATNodeOutput(order)
        
        SAT_Left    = zeros(T,nnodes)
        SAT_Right   = zeros(T,nnodes)

        u_Left      = zeros(T,nnodes)
        u_Right     = zeros(T,nnodes)

        if length(BC) == 2
            BCL = BC.Left.type
            BCR = BC.Right.type
        elseif length(BC) == 1
            BCL = BCR = Periodic
        end

        new{T}(BCL,BCR,SAT_Left,SAT_Right,u_Left,u_Right,[0.0],[0.0])

    end
end

"""
    BoundaryData2D{T}
Data structure for storage of SATs in 2 dimensional problems
"""
struct BoundaryData2D{T} <: BoundaryStorage{T,2}

    Type_Left    :: BoundaryConditionType
    Type_Right   :: BoundaryConditionType
    Type_Up      :: BoundaryConditionType
    Type_Down    :: BoundaryConditionType

    SAT_Left    :: AbstractArray{T} #Same as u_ but for SAT
    SAT_Right   :: AbstractArray{T} #Same as u_ but for SAT
    SAT_Up      :: AbstractArray{T} #Same as u_ but for SAT
    SAT_Down    :: AbstractArray{T} #Same as u_ but for SAT

    u_Left      :: AbstractArray{T} #Solution along boundary with size determined by derivative order
    u_Right     :: AbstractArray{T} #Solution along boundary with size determined by derivative order
    u_Up        :: AbstractArray{T} #Solution along boundary with size determined by derivative order
    u_Down      :: AbstractArray{T} #Solution along boundary with size determined by derivative order

    RHS_Left    :: AbstractArray{T}
    RHS_Right   :: AbstractArray{T}
    RHS_Up      :: AbstractArray{T}
    RHS_Down    :: AbstractArray{T}

    function BoundaryData2D{T}(BC::NamedTuple,grid::Grid2D,order::Int) where {T}

        nnodes = SATNodeOutput(order)
        ny = grid.ny
        nx = grid.nx

        SAT_Left =  zeros(T,(nnodes,ny)) #x
        SAT_Right = zeros(T,(nnodes,ny)) #x
        SAT_Up =    zeros(T,(nx,nnodes)) #y 
        SAT_Down =  zeros(T,(nx,nnodes)) #y

        u_Left =    zeros(T,(nnodes,ny)) #x
        u_Right =   zeros(T,(nnodes,ny)) #x
        u_Up =      zeros(T,(nx,nnodes)) #y 
        u_Down =    zeros(T,(nx,nnodes)) #y

        RHS_Left =    zeros(T,(1,ny)) #x
        RHS_Right =   zeros(T,(1,ny)) #x
        RHS_Up =      zeros(T,(nx,1)) #y 
        RHS_Down =    zeros(T,(nx,1)) #y

        new{T}(BC.Left.type,BC.Right.type,BC.Up.type,BC.Down.type,
            SAT_Left,SAT_Right,SAT_Up,SAT_Down,
            u_Left,u_Right,u_Up,u_Down,
            RHS_Left,RHS_Right,RHS_Up,RHS_Down)

    end
end




function BuildStorage(Prob,grid,Δt)
    DB = DataBlock{Float64}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.K)
    CG = ConjGradBlock{Float64}(grid.n)
    return DB,CG
end





"""
    copyUtoSAT!
Moves data from the solution `u` at a given boundary to the `SAT_` field in `BoundaryStorage` structs. Or moves all data to `SAT_` fields.
"""
function copyUtoSAT!(SAT::AbstractArray,u::AbstractArray,side::NodeType,order::Int)
    nnodes = SATNodeOutput(order)
    if side == Left
        SAT .= u[1:nnodes,:]
    elseif side == Right
        SAT .= u[end-nnodes+1:end,:]
    elseif side == Up
        SAT .= u[:,1:nnodes]
    elseif side == Down
        SAT .= u[:,end-nnodes+1:end]
    end
end
function copyUtoSAT!(Bound::BoundaryStorage,u::AbstractArray,order::Int)
    copyUtoSAT!(Bound.SAT_Left,u,Left,order)
    copyUtoSAT!(Bound.SAT_Right,u,Right,order)
    if typeof(Bound) <: BoundaryStorage{T,2} where T
        copyUtoSAT!(Bound.SAT_Up,u,Up,order)
        copyUtoSAT!(Bound.SAT_Down,u,Down,order)
    end
end

"""
    copySATtoU!
Moves data from the `SAT_` field  on the given side in `BoundaryStorage` to `u`. Or moves all data from `u` to `SAT_` fields.
"""
function copySATtoU!(u::AbstractArray,SAT::AbstractArray,side::NodeType,order::Int)
    nnodes = SATNodeOutput(order)
    if side == Left
        u[1:nnodes,:] .= SAT
    elseif side == Right
        u[end-nnodes+1:end,:] .= SAT
    elseif side == Up
        u[:,1:nnodes] .= SAT
    elseif side == Down
        u[:,end-nnodes+1:end] .= SAT
    end
end
function copySATtoU!(u::AbstractArray,Bound::BoundaryStorage,order::Int)
    copySATtoU!(u,Bound.SAT_Left,Left,order)
    copySATtoU!(u,Bound.SAT_Right,Right,order)
    if typeof(Bound) <: BoundaryStorage{T,2} where T
        copySATtoU!(u,Bound.SAT_Up,Up,order)
        copySATtoU!(u,Bound.SAT_Down,Down,order)
    end
end

"""
    addSATtoU!
Add data from the `SAT_` field  on the given side in `BoundaryStorage` to `u`. Or add all data from `u` to `SAT_` fields.
"""
function addSATtoU!(u::AbstractArray,SAT::AbstractArray,side::NodeType,order::Int)
    nnodes = SATNodeOutput(order)
    if side == Left
        u[1:nnodes,:] .+= SAT
    elseif side == Right
        u[end-nnodes+1:end,:] .+= SAT
    elseif side == Up
        u[:,1:nnodes] .+= SAT
    elseif side == Down
        u[:,end-nnodes+1:end] .+= SAT
    end
end
function addSATtoU!(u::AbstractArray,Bound::BoundaryStorage,order::Int)
    addSATtoU!(u,Bound.SAT_Left,Left,order)
    addSATtoU!(u,Bound.SAT_Right,Right,order)
    if typeof(Bound) <: BoundaryStorage{T,2} where T
        addSATtoU!(u,Bound.SAT_Up,Up,order)
        addSATtoU!(u,Bound.SAT_Down,Down,order)
    end
end

"""
    addSource!
Add the source term `F(x,y,t)` to the array `u`.
"""
function addSource! end
function addSource!(F::Function,u::AbstractArray{T},grid::Grid2D{T},t::T,Δt) where T
    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] += Δt*F(grid.gridx[i],grid.gridy[j],t)
        end
    end
end
function addSource!(F::Function,u::AbstractArray{T},grid::Grid1D{T},t::T,Δt::T) where T
    u .+= Δt*F.(grid.grid,t)
end

"""
    setBoundary!
"""
function setBoundary!(RHS::Function,Bound::AbstractArray{T},grid::AbstractArray{T},n::Int,t::T,Δt::T) where T
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end
function setBoundary!(RHS::Function,Bound::AbstractArray{T},t::T,Δt::T) where T
        Bound[1] = Δt*RHS(t)
end

# function setBoundary!(BC::NamedTuple,BoundStore::BoundaryStorage{T,N},grid::AbstractArray{T},t::T,Δt::T,side::NodeType)
#     if side == Left
#         setBoundary!(BC.Left.RHS,BoundStore.RHS_Left,grid,n,t,Δt)
#     elseif side == Right
#     elseif side == Up
#     elseif side == Down
#     end
# end
# function setBoundaries!(BCs::NamedTuple,BoundaryStore::BoundaryStorage{T,N},grid::GridType,t::T,Δt::T,sides::NodeType...)
#     if GetDim(BoundaryStore) == 2
#     end
# end