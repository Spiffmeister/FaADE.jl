

# abstract type Solution{T<:AbstractFloat} end
abstract type DataBlockType{T<:AbstractFloat} end
abstract type BoundaryStorage{T<:AbstractFloat} end


#========== WHOLE PROBLEM DATA ==========#
"""
    DataBlock{T}
Passed around internally between functions. Only contains the current timestep
"""
struct DataBlock{T} <: DataBlockType{T}
    grid        :: GridType
    u           :: AbstractArray{T}
    uₓₓ         :: AbstractArray{T}
    boundary    :: BoundaryStorage
    Δt          :: T
    function DataBlock{T}(grid,u,uₓₓ,boundary,Δt)
        new(grid,u,uₓₓ,boundary,Δt)
    end
end

"""
    DataBlock
External constructor for `DataBlock{T}` for 1D and 2D problems.
"""
function DataBlock(grid::Grid,Δt::T,order::Int,boundaries...) where T
    # Build tuple of boundary types, ensure that Periodic boundaries do not result in too few types
    BTypes = []
    for bound in boundaries
        if bound.type != Periodic
            push!(BTypes,bound.type)
        else
            push!(BTypes,Periodic)
            push!(BTypes,Periodic)
        end
    end
    BTypes = Tuple(BTypes)

    # If grid is 1D or 2D construct the right DataBlock
    if typeof(grid) <: Grid1D
        u   = zeros(T,grid.n)
        uₓₓ = zeros(T,grid.n)
        BStor = BoundaryData1D{T}(BTypes,grid.n,order)
    elseif typeof(grid) <: Grid2D
        u   = zeros(T,(grid.nx,grid.ny))
        uₓₓ = zeros(T,(grid.nx,grid.ny))
        BStor = BoundaryData2D{T}(BTypes,grid.nx,grid.ny,order)
    end
    
    DBlock = DataBlock{T}(grid,u,uₓₓ,Bstor,Δt)
    return DBlock
end



#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D{T}
Data structure for storage of SATs in 1 dimensional problems
"""
struct BoundaryData1D{T} <: BoundaryStorage{T}
    Type_Left   :: BoundaryConditionType
    Type_Right  :: BoundaryConditionType

    SAT_Left    :: AbstractArray{T}
    SAT_Right   :: AbstractArray{T}

    u_Left      :: AbstractArray{T}
    u_Right     :: AbstractArray{T}

    function BoundaryData1D{T}(type::Tuple{BoundaryConditionType,BoundaryConditionType},n::Int,order::Int) where T

        nnodes = SATNodeOutput(order)

        SAT_Left    = zeros(T,nnodes)
        SAT_Right   = zeros(T,nnodes)

        u_Left      = zeros(T,nnodes)
        u_Right     = zeros(T,nnodes)

        new{T}(type[1],type[2],SAT_Left,SAT_Right,u_Left,u_Right)

    end
end

"""
    BoundaryData2D{T}
Data structure for storage of SATs in 2 dimensional problems
"""
struct BoundaryData2D{T} <: BoundaryStorage{T}

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

    function BoundaryData2D{T}(type::Tuple{BoundaryConditionType,BoundaryConditionType,BoundaryConditionType,BoundaryConditionType},nx::Int,ny::Int,order::Int) where T

        nnodes = SATNodeOutput(order)

        SAT_Left =  zeros(T,(nnodes,ny)) #x
        SAT_Right = zeros(T,(nnodes,ny)) #x
        SAT_Up =    zeros(T,(nx,nnodes)) #y 
        SAT_Down =  zeros(T,(nx,nnodes)) #y

        u_Left =    zeros(T,(nnodes,ny)) #x
        u_Right =   zeros(T,(nnodes,ny)) #x
        u_Up =      zeros(T,(nx,nnodes)) #y 
        u_Down =    zeros(T,(nx,nnodes)) #y

        new{T}(type[1],type[2],type[3],type[4],
            SAT_Left,SAT_Right,SAT_Up,SAT_Down,
            u_Left,u_Right,u_Up,u_Down)

    end
end







"""
    copyUtoSAT
Moves data from the solution `u` at a given boundary to the `SAT_` field in `BoundaryStorage` structs.
"""
function copyUtoSAT(SAT::AbstractArray,u::AbstractArray,side::NodeType,nnodes::Int)
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



