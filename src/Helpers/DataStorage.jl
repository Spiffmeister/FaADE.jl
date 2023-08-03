

abstract type LocalDataBlockType{dtype,DIM,atype} <: DataBlockType{dtype,DIM,atype} end


# abstract type Solution{T<:AbstractFloat} end


#========== WHOLE PROBLEM DATA ==========#
"""
    DataBlock
Passed around internally between functions. Only contains data required for current timestep.
"""
mutable struct DataBlock{T,N,AT,KT<:Union{AT,Vector{AT}}} <: LocalDataBlockType{T,N,AT}
    grid        :: GridType
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BoundaryStorage{T,N,AT}
    t           :: T
    Δt          :: T
    Δu          :: T
    function DataBlock{T}(
            PDE::PDEProblem,
            grid::GridType{T,D},
            Δt::T,
            K::Function...) where {T,D}
    
        # If grid is 1D or 2D construct the right DataBlock
        if typeof(grid) <: Grid1D
            u   = zeros(T,grid.n)
            uₙ₊₁ = zeros(T,grid.n)
            BStor = BoundaryData1D{T}(PDE.BoundaryConditions,PDE.order)

            DiffCoeff = zeros(T,grid.n)
            setCoefficient!(PDE.K,DiffCoeff,grid)

            dim = 1
            
        elseif typeof(grid) <: Grid2D
            u   = zeros(T,(grid.nx,grid.ny))
            uₙ₊₁ = zeros(T,(grid.nx,grid.ny))
            BStor = BoundaryData2D{T}(PDE.BoundaryConditions,grid,PDE.order)

            DiffCoeff = [zeros(T,(grid.nx,grid.ny)),zeros(T,(grid.nx,grid.ny))]
            setCoefficient!(PDE.Kx,DiffCoeff[1],grid)
            setCoefficient!(PDE.Ky,DiffCoeff[2],grid)

            dim = 2

        end
        new{T,dim,typeof(u),typeof(DiffCoeff)}(grid,u,uₙ₊₁,DiffCoeff,BStor,0,Δt,0.0)
    end
end



#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D
Data structure for storage of SATs in 1 dimensional problems
"""
mutable struct BoundaryData1D{T,AT} <: BoundaryStorage{T,1, AT}
    Type_Left   :: BoundaryConditionType
    Type_Right  :: BoundaryConditionType

    SAT_Left    :: AT
    SAT_Right   :: AT

    u_Left      :: AT
    u_Right     :: AT

    RHS_Left    :: AT
    RHS_Right   :: AT

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

        new{T,typeof(u_Left)}(BCL,BCR,SAT_Left,SAT_Right,u_Left,u_Right,[0.0],[0.0])

    end
end

"""
    BoundaryData2D
Data structure for storage of SATs in 2 dimensional problems
"""
struct BoundaryData2D{T,AT} <: BoundaryStorage{T,2, AT}

    Type_Left    :: BoundaryConditionType
    Type_Right   :: BoundaryConditionType
    Type_Up      :: BoundaryConditionType
    Type_Down    :: BoundaryConditionType

    SAT_Left    :: AT #Same as u_ but for SAT
    SAT_Right   :: AT #Same as u_ but for SAT
    SAT_Up      :: AT #Same as u_ but for SAT
    SAT_Down    :: AT #Same as u_ but for SAT

    u_Left      :: AT #Solution along boundary with size determined by derivative order
    u_Right     :: AT #Solution along boundary with size determined by derivative order
    u_Up        :: AT #Solution along boundary with size determined by derivative order
    u_Down      :: AT #Solution along boundary with size determined by derivative order

    RHS_Left    :: AT
    RHS_Right   :: AT
    RHS_Up      :: AT
    RHS_Down    :: AT

    function BoundaryData2D{T}(BC::NamedTuple,grid::Grid2D,order::Int) where T

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

        new{T,typeof(u_Left)}(BC.Left.type,BC.Right.type,BC.Up.type,BC.Down.type,
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




#=
"""
    copyUtoSAT!
Moves data from the solution `u` at a given boundary to the `SAT_` field in `BoundaryStorage` structs. Or moves all data to `SAT_` fields.
"""
function copyUtoSAT! end
function copyUtoSAT!(SAT::AT,u::AT,side::NodeType,order::Int) where AT
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
function copyUtoSAT!(Bound::BoundaryStorage{T,N,AT},u::AT,order::Int) where {T,N,AT}
    copyUtoSAT!(Bound.SAT_Left,u,Left,order)
    copyUtoSAT!(Bound.SAT_Right,u,Right,order)
    if typeof(Bound) <: BoundaryStorage{T,2,AT}
        copyUtoSAT!(Bound.SAT_Up,u,Up,order)
        copyUtoSAT!(Bound.SAT_Down,u,Down,order)
    end
end

"""
    copySATtoU!
Moves data from the `SAT_` field  on the given side in `BoundaryStorage` to `u`. Or moves all data from `u` to `SAT_` fields.
"""
function copySATtoU! end
function copySATtoU!(u::AT,SAT::AT,side::NodeType,order::Int) where AT
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
function copySATtoU!(u::AT,Bound::BoundaryStorage{T,N,AT},order::Int) where {T,N,AT}
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
function addSATtoU! end
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
function addSource!(S::Nothing,u::AbstractArray{TT},grid::LocalGridType,t::TT,Δt::TT) where TT end
function addSource!(F::Function,u::AbstractArray{TT},grid::Grid1D{TT},t::TT,Δt::TT) where TT
    u .+= Δt*F.(grid.grid,t)
end
function addSource!(F::Function,u::AbstractArray{TT},grid::Grid2D{TT},t::TT,Δt::TT) where TT
    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] += Δt*F(grid.gridx[i],grid.gridy[j],t)
        end
    end
end


"""
    setBoundary!
Sets the value of the boundary.
"""
function setBoundary! end
function setBoundary!(RHS::F,Bound::AT,t::T,Δt::T) where {F,AT,T}
    Bound[1] = Δt*RHS(t)
end
function setBoundary!(RHS::F,Bound::AT,grid::Vector{T},n::Int,t::T,Δt::T) where {F,AT,T}
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end



"""
    setCoefficient!
Sets the diffusion coefficient
"""
function setCoefficient! end
function setCoefficient!(K::Function,κ::AbstractArray,grid::Grid1D)
    for i = 1:grid.n
        κ[i] = K(grid.grid[i])
    end
end
function setCoefficient!(K::Function,κ::AbstractArray,grid::Grid2D)
    for i = 1:grid.nx
        for j = 1:grid.ny
            κ[i,j] = K(grid.gridx[i],grid.gridy[j])
        end
    end
end
function setCoefficient!(DC::DiffusionCoefficient{F},κ::AbstractArray,grid::Grid2D) where {F<:Function}
    setCoefficient!(DC.coeff,κ,grid)
end
function setCoefficient!(DC::DiffusionCoefficient{TT},κ::AbstractArray{TT},grid::GridType) where {TT}
    κ .= DC.coeff
end

# function setCoefficient!(P::newPDEProblem,D::LocalDataBlockType,grid::GridType)
#     setCoefficient!(P.Kx,D.Kx,grid)
#     setCoefficient!(P.Ky,D.Ky,grid)
# end
=#