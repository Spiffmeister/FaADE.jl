

# abstract type Solution{T<:AbstractFloat} end
abstract type DataBlockType{T<:AbstractFloat,N} end
abstract type BoundaryStorage{T<:AbstractFloat} end


#========== WHOLE PROBLEM DATA ==========#
"""
    DataBlock{T}
Passed around internally between functions. Only contains data required for current timestep.
"""
struct DataBlock{T,N} <: DataBlockType{T,N}
    dim         :: Int
    grid        :: GridType
    u           :: AbstractArray{T}
    uₙ₊₁        :: AbstractArray{T}
    K           :: Union{Vector,Vector{AbstractArray{T}}}
    boundary    :: BoundaryStorage
    Δt          :: T
    function DataBlock{T}(
            boundaries::NamedTuple,
            grid::GridType,
            Δt::T,
            order::Int,
            K::AbstractArray{T}) where {T}
    
        # If grid is 1D or 2D construct the right DataBlock
        if typeof(grid) <: Grid1D
            u   = zeros(T,grid.n)
            uₓₓ = zeros(T,grid.n)
            BStor = BoundaryData1D{T}(boundaries,order)
            DiffCoeff = K
            dim = 1
        elseif typeof(grid) <: Grid2D
            u   = zeros(T,(grid.nx,grid.ny))
            uₓₓ = zeros(T,(grid.nx,grid.ny))
            BStor = BoundaryData2D{T}(boundaries,grid.nx,grid.ny,order)
            DiffCoeff = [K]
            dim = 2
        end
        new{T,dim}(dim,grid,u,uₓₓ,DiffCoeff,BStor,Δt)
    end
end


#========== WHOLE PROBLEM DATA ==========#
struct ConjGradBlock{T,N} <: DataBlockType{T,N}
    b   :: AbstractArray{T,N} # b = uⁿ⁺¹ + F
    rₖ  :: AbstractArray{T,N} # (uⁿ⁺¹ - Δt*uₓₓⁿ⁺¹) - b
    Adₖ :: AbstractArray{T,N} # Adₖ = dₖ - Δt*D(dₖ)
    Drₖ :: AbstractArray{T,N} # Drₖ = rₖ - Δt*D(rₖ)
    dₖ  :: AbstractArray{T,N} # dₖ = -rₖ, -rₖ .+ βₖ*dₖ

    converged :: Bool

    function ConjGradBlock{T}(n::Int...) where T
        b   = zeros(T,n)
        rₖ  = zeros(T,n)
        Adₖ = zeros(T,n)
        Arₖ = zeros(T,n)
        dₖ  = zeros(T,n)

        dims = length(n)

        new{T,dims}(b, rₖ, Adₖ, Arₖ, dₖ, true)
    end
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

    RHS_Left    :: Float64
    RHS_Right   :: Float64

    function BoundaryData1D{T}(BC::NamedTuple,order::Int) where {T}

        nnodes = SATNodeOutput(order)

        SAT_Left    = zeros(T,nnodes)
        SAT_Right   = zeros(T,nnodes)

        u_Left      = zeros(T,nnodes)
        u_Right     = zeros(T,nnodes)

        new{T}(BC.Left.type,BC.Right.type,SAT_Left,SAT_Right,u_Left,u_Right,0.0,0.0)

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

    function BoundaryData2D{T}(type::Tuple{D},nx::Int,ny::Int,order::Int) where {T,D<:BoundaryConditionData}

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


# function copySATtoU!(DBlock::DataBlock{T,1},order::Int) where T
#     copySATtoU!(u,SAT_Left)
#     copySATtoU!(u,SAT_Left)
# end
# function copySATtoU!(DBlock::DataBlock{T,2},order::Int) where T
#     copySATtoU!(u,SAT_Left)
#     copySATtoU!(u,SAT_Left)
#     copySATtoU!(u,SAT_Left)
#     copySATtoU!(u,SAT_Left)
# end