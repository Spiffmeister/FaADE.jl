
abstract type newDataBlockType{dtype,DIM} end
abstract type newLocalDataBlockType{dtype,DIM} <: newDataBlockType{dtype,DIM} end



"""
    newDataBlock
Data structure for multiblock problems
"""
struct newDataBlock{TT<:Real,
    DIM,
    DT <: Union{newLocalDataBlockType,Vector{newLocalDataBlockType}},
        } <: newDataBlockType{TT,DIM}
    Data    :: DT
    t       :: TT
    Δt      :: TT
    Δu      :: TT

    function newDataBlock(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = [LocalDataBlock(P,Grid,Δt) for Grid in G.Grids]
        new{TT,DIM,typeof(DTA)}(DTA,t,Δt,0.0)
    end
    function newDataBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = LocalDataBlock(P,G)
        new{TT,DIM,typeof(DTA)}(DTA,t,Δt,0.0)
    end
end




#========== NEW  DATA ==============#
"""
    LocalDataBlock
"""
mutable struct LocalDataBlock{TT<:Real,
        DIMS,
        AT  <: AbstractArray{TT},
        KT  <: Union{Vector{TT},Vector{Matrix{TT}}},
        BT  <: BoundaryStorage{TT,DIMS,AT}
        } <: newLocalDataBlockType{TT,DIMS}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BT
    Δu          :: TT
    function LocalDataBlock(
            P::newPDEProblem{TT,DIMS},
            G::LocalGridType) where {TT,DIMS}

        u   = zeros(TT,size(G))
        uₙ₊₁ = zeros(TT,size(G))


        if typeof(G) <: Grid1D
            BStor = newBoundaryData1D(G,P.order)
            
            DiffCoeff = zeros(TT,size(G))
            setCoefficient!(P.K,DiffCoeff,G)

        elseif typeof(G) <: Grid2D
            
            BStor = newBoundaryData2D(G,P.order)

            DiffCoeff = [zeros(TT,size(G)),zeros(TT,size(G))]
            setCoefficient!(P.Kx,DiffCoeff[1],G)
            setCoefficient!(P.Ky,DiffCoeff[2],G)

        end

        new{TT,DIMS,typeof(u),typeof(DiffCoeff),typeof(BStor)}(u,uₙ₊₁,DiffCoeff,BStor,0.0)
    end

end



#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D
Data structure for storage of SATs in 1 dimensional problems
"""
mutable struct newBoundaryData1D{TT,
        AT} <: BoundaryStorage{TT,1,AT}

    SAT_Left    :: AT
    SAT_Right   :: AT

    u_Left      :: AT
    u_Right     :: AT

    RHS_Left    :: AT
    RHS_Right   :: AT

    function newBoundaryData1D(G::Grid1D{TT,MET},order::DerivativeOrder{O}) where {TT,O,MET}

        nnodes = SATNodeOutput(O)

        SAT_Left    = zeros(TT,nnodes)
        SAT_Right   = zeros(TT,nnodes)

        u_Left      = zeros(TT,nnodes)
        u_Right     = zeros(TT,nnodes)

        new{TT,typeof(u_Left)}(SAT_Left,SAT_Right,u_Left,u_Right,[0.0],[0.0])
    end
end

"""
    BoundaryData2D
Data structure for storage of SATs in 2 dimensional problems
"""
struct newBoundaryData2D{TT,AT} <: BoundaryStorage{TT,2,AT}

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

    function newBoundaryData2D(grid::Grid2D{TT,MET},order::DerivativeOrder{O}) where {TT,O,MET}

        nnodes = SATNodeOutput(order)
        ny = grid.ny
        nx = grid.nx

        SAT_Left =  zeros(TT,(nnodes,ny)) #x
        SAT_Right = zeros(TT,(nnodes,ny)) #x
        SAT_Up =    zeros(TT,(nx,nnodes)) #y 
        SAT_Down =  zeros(TT,(nx,nnodes)) #y

        u_Left =    zeros(TT,(nnodes,ny)) #x
        u_Right =   zeros(TT,(nnodes,ny)) #x
        u_Up =      zeros(TT,(nx,nnodes)) #y 
        u_Down =    zeros(TT,(nx,nnodes)) #y

        RHS_Left =    zeros(TT,(1,ny)) #x
        RHS_Right =   zeros(TT,(1,ny)) #x
        RHS_Up =      zeros(TT,(nx,1)) #y 
        RHS_Down =    zeros(TT,(nx,1)) #y

        new{TT,typeof(u_Left)}(SAT_Left,SAT_Right,SAT_Up,SAT_Down,
            u_Left,u_Right,u_Up,u_Down,
            RHS_Left,RHS_Right,RHS_Up,RHS_Down)
    end
end







# function joining()
#     for joint in joints
#         if i ∈ j
#             ii = findfirst(x->x==i,joint)
#             if ii == 2
#                 -1*joint
#             end

#             direction = (1-i)*joint[end]
#         end
#     end
# end



#=
function Base.mul(i::Integer,Node::NodeType{T,D}) where {T,D}
    if ispositive(i)
        return Node
    else
        if NodeType == Left
            return Right
        end
        if NodeType == Right
            return Left
        end
        if NodeType == Up
            return Down
        end
        if NodeType == Down
            return Up
        end
    end
end
=#












"""
    setBoundary!
Sets the value of the boundary.
"""
function setBoundary! end
function setBoundary!(RHS::F,Bound::AT,grid::Vector{T},n::Int,t::T,Δt::T) where {F,AT,T}
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end
function setBoundary!(RHS::F,Bound::AT,t::T,Δt::T) where {F,AT,T}
    Bound[1] = Δt*RHS(t)
end
function setBoundaries(D::newDataBlockType,B::SATBoundaries,G::GridType,t::TT,Δt::TT) where TT
    # setBoundary!(B.Boundary.RHS,D.Data.RHS_Left,G,t,Δt)
end


# function fetchBoudnary! end
# function fetchBoudnary!(D1,D2) where {F,AT,T}
# end



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




"""
    CommBoundary(D,J)
"""
function CommBoundary(D::newDataBlock,G::GridMultiBlock)
    for J in G.Joint
        CommBoundary(D.Data[J[1]],D.Data[J[2]],J[3])
    end
end
function CommBoundary(D1::LocalDataBlock,D2::LocalDataBlock,N::NodeType)
    D1.boundary .= D2.u()
    D2.boundary .= D1.u()
end