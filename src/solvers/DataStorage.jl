
abstract type newDataBlockType{dtype,DIM} end
abstract type newLocalDataBlockType{dtype,DIM} <: newDataBlockType{dtype,DIM} end


struct DerivativeInfo{TT}
    order   :: DerivativeOrder
    nx      :: Int64
    ny      :: Int64
    Δx      :: TT
    Δy      :: TT
end

mutable struct StepConfig{TT}
    t   :: TT
    Δt  :: TT
    Δu  :: TT
end

"""
    newDataBlock
Data structure for multiblock problems
"""
struct newDataBlock{TT<:Real,
    DIM,
    NBLOCK,
    TDBLOCK <: Union{newLocalDataBlockType,Vector{newLocalDataBlockType}},
        } <: newDataBlockType{TT,DIM}
    Data    :: TDBLOCK
    SC      :: StepConfig{TT}

    function newDataBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = LocalDataBlock(P,G)
        new{TT,DIM,1,typeof(DTA)}(DTA,StepConfig(t,Δt,0.0))
    end
    function newDataBlock(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = [LocalDataBlock(P,Grid,Δt) for Grid in G.Grids]
        new{TT,DIM,NG,typeof(DTA)}(DTA,t,Δt,0.0)
    end
end




#========== NEW  DATA ==============#
"""
    LocalDataBlock
"""
mutable struct LocalDataBlock{TT<:Real,
        DIM,
        AT  <: AbstractArray{TT},
        KT  <: Union{Vector{TT},Vector{Matrix{TT}}},
        BT  <: BoundaryStorage{TT,DIM,AT},
        DT  <: DerivativeOperator
        } <: newLocalDataBlockType{TT,DIM}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BT
    Derivative  :: DT
    Δu          :: TT

    function LocalDataBlock(P::newPDEProblem{TT,DIM},G::LocalGridType) where {TT,DIM}

        u   = zeros(TT,size(G))
        uₙ₊₁ = zeros(TT,size(G))

        if DIM == 1
            BStor = newBoundaryData1D(G,P.order,P.BoundaryConditions.BoundaryLeft,P.BoundaryConditions.BoundaryRight)
            DiffCoeff = zeros(TT,size(G))
            setCoefficient!(P.K,DiffCoeff,G)
            D = DerivativeOperator{TT,1,true,false,false}(P.order,G.n,0,G.Δx,TT(0))
        elseif DIM == 2
            BStor = newBoundaryData2D(G,P.order)
            DiffCoeff = [zeros(TT,size(G)),zeros(TT,size(G))]
            setCoefficient!(P.Kx,DiffCoeff[1],G)
            setCoefficient!(P.Ky,DiffCoeff[2],G)
            D = DerivativeOperator{TT,2,true,false,false}(P.order,G.nx,G.ny,G.Δx,G.Δy)
        end

        new{TT,DIM,typeof(u),typeof(DiffCoeff),typeof(BStor),typeof(D)}(u,uₙ₊₁,DiffCoeff,BStor,D,0.0)
    end
    function LocalDataBlock(P::newPDEProblem{TT,DIM},G::GridMultiBlock{TT,DIM,NG},ig::Integer) where {TT,DIM,NG}
        # (My index, neighbours index, Boundary)

        u   = zeros(TT,size(G.Grids[ig]))
        uₙ₊₁ = zeros(TT,size(G.Grids[ig]))

        
        if DIM == 1
            D = DerivativeInfo(P.order,G.nx,0,G.Δx,TT(0))
            if ig != 1
                BCL = SAT_Split()
            else
                BCL = P.Boundaries.BoundaryLeft
            end
            if ig != NG
                BCL = SAT_Split()
            else
                BCR = P.Boundaries.BoundaryRight
            end
            BStor = newBoundaryData1D(G,P.order,BCL,BCR)
        elseif DIM == 2
        end
        new{TT,DIM,typeof(u),typeof(DiffCoeff),typeof(BStor),typeof(D)}(u,uₙ₊₁,DiffCoeff,Bstor,D,0.0)
    end

end

"""
    newBoundaryData
"""
struct newBoundaryData{TT,SAT<:SimultanousApproximationTerm,AT}
    BC      :: SAT
    u       :: AT
    RHS     :: AT
    function newBoundaryData(G::Grid1D{TT},order::DerivativeOrder{O},BC) where {TT,O}
        nnodes = SATNodeOutput(O)
        u = zeros(TT,nnodes)
        new{TT,typeof(BC),typeof(u)}(BC,u,[0.0])
    end
end

#========== BOUNDARY DATA ==========#
"""
    BoundaryData1D
Data structure for storage of SATs in 1 dimensional problems
"""
struct newBoundaryData1D{TT,
        SATL <: SimultanousApproximationTerm,
        SATR <: SimultanousApproximationTerm,
        AT} <: BoundaryStorage{TT,1,AT}

    Left        :: SATL
    Right       :: SATR

    u_Left      :: AT
    u_Right     :: AT

    RHS_Left    :: AT
    RHS_Right   :: AT

    function newBoundaryData1D(G::Grid1D{TT,MET},order::DerivativeOrder{O},BCL,BCR) where {TT,O,MET}

        nnodes = SATNodeOutput(O)

        u_Left      = zeros(TT,nnodes)
        u_Right     = zeros(TT,nnodes)

        new{TT,typeof(BCL),typeof(BCR),typeof(u_Left)}(BCL,BCR,u_Left,u_Right,[0.0],[0.0])
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

    # Left    :: B1
    # Right   :: B1
    # Up      :: B1
    # Down    :: B1

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
function setBoundaries(D::newDataBlockType,G::Grid1D,t::TT,Δt::TT) where TT
    setBoundary!(D.Data.boundary.Left.RHS,D.Data.boundary.RHS_Left,t,Δt)
    setBoundary!(D.Data.boundary.Right.RHS,D.Data.boundary.RHS_Right,t,Δt)
end
function setBoundaries(D::DataBlockType,B::SATBoundaries,G::Grid2D,t::TT,Δt::TT) where TT
    setBoundary!(B.BoundaryLeft.RHS,D.Data.RHS_Left,G.nx,G.gridx,t,Δt)
    setBoundary!(B.BoundaryRight.RHS,D.Data.RHS_Right,G.nx,G.gridx,t,Δt)
    setBoundary!(B.BoundaryUp.RHS,D.Data.RHS_Up,G.ny,G.gridy,t,Δt)
    setBoundary!(B.BoundaryDown.RHS,D.Data.RHS_Down,G.ny,G.gridy,t,Δt)
end

function setBoundaries(D::newDataBlockType,B::SATBoundaries,G::GridType,t::TT,Δt::TT) where TT
    for DB in D.Data
        DB.Boundary
    end
end

function setBoundary!()
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
        κ[i] = K(grid[i])
    end
end
function setCoefficient!(K::Function,κ::AbstractArray,grid::Grid2D)
    for i = 1:grid.nx
        for j = 1:grid.ny
            κ[i,j] = K(grid[i,j]...)
        end
    end
end
# function setCoefficient!(K::Function,κ::AbstractArray,grid::GridType)
#     for i in eachindex(grid)
#         κ[i] = K(grid[i]...)
#     end
# end
function setCoefficient!(DC::DiffusionCoefficient{F},κ::AbstractArray,grid::LocalGridType) where {F<:Function}
    setCoefficient!(DC.coeff,κ,grid)
end
function setCoefficient!(DC::DiffusionCoefficient{TT},κ::AbstractArray{TT},grid::LocalGridType) where {TT<:Real}
    κ .= DC.coeff
end


"""
    applySATs
"""
function applySATs(CG,D::newDataBlock{TT,1},mode::SATMode) where {TT}
    applySAT!(D.Data.boundary.Left, CG.b, D.Data.boundary.RHS_Left, D.Data.K, mode)
    applySAT!(D.Data.boundary.Right, CG.b, D.Data.boundary.RHS_Right, D.Data.K, mode)
end

#== INTER-BLOCK COMMUNICATION ==#
"""
    CommBoundary(dest,source,Joint,inds)
Given a `Joint` between grid regions, send the data between blocks
"""
function CommBoundaries(DB::DataBlock{TT,DIM,NBLOCK}) where {TT,DIM,NBLOCK}
    # Iterate over blocks
    # Foreach boundary
    # If boundary Periodic or Split
    #   Communicate boundary data
    # else
    #   do nothing
    # end
    for i in 1:NBLOCK
        for joint in DB[i].boundary.joints
                CommBoundary(D,DB[i].boundary)
        end
    end
end
function CommBoundary(B::BoundaryStorage,DB{TT,DIM,NBLOCK})
end
function CommBoundary(dest::AT,source::AT,inds::Tuple) where {AT}
    dest .= source[inds]
end
function CommBoundary(D1::LocalDataBlock,D2::LocalDataBlock,N::NodeType)
    D1.boundary .= D2.u()
    D2.boundary .= D1.u()
end