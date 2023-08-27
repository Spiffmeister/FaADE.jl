
abstract type newDataBlockType{dtype,DIM} end
abstract type newLocalDataBlockType{dtype,DIM} <: newDataBlockType{dtype,DIM} end

#=
"""
    BlockRelations
Object for communication handling between blocks
"""
struct BlockRelations{TJ<:Union{Nothing,Tuple}}
    joints

    function BlockRelations(G::GridMultiBlock{TT,1}) where TT

        new{Tuple}(G.Joint)
    end
end
=#

"""
    StepConfig{TT}
Configuration and stats for the current step for the current block
"""
mutable struct StepConfig{TT}
    t   :: TT
    Δt  :: TT
    Δu  :: TT
end


mutable struct DiffCoeff
    K   :: Union{Function,Real}
    KA  :: Array
end


#========== NEW  DATA ==============#



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
    LocalDataBlock
"""
mutable struct LocalDataBlock{TT<:Real,
        DIM,
        AT  <: AbstractArray{TT},
        KT  <: Union{Vector{TT},Vector{Matrix{TT}}},
        BT  <: BoundaryStorage{TT,DIM,AT},
        DT  <: DerivativeOperator,
        } <: newLocalDataBlockType{TT,DIM}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BT
    Derivative  :: DT



    SC          :: StepConfig{TT}
end
function LocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}

    u   = zeros(TT,size(G))
    uₙ₊₁ = zeros(TT,size(G))

    BStor = newBoundaryData1D(G,P.order,P.BoundaryConditions.BoundaryLeft,P.BoundaryConditions.BoundaryRight)
    DiffCoeff = zeros(TT,size(G))
    setCoefficient!(P.K,DiffCoeff,G)
    D = DerivativeOperator{TT,1,true,false,false}(P.order,G.n,0,G.Δx,TT(0))

    SC = StepConfig{TT}(TT(0),TT(0),TT(0))

    LocalDataBlock{TT,1,typeof(u),typeof(DiffCoeff),typeof(BStor),typeof(D)}(u,uₙ₊₁,DiffCoeff,BStor,D,SC)
end
function LocalDataBlock(P::newPDEProblem{TT,2},G::LocalGridType) where {TT}

    u   = zeros(TT,size(G))
    uₙ₊₁ = zeros(TT,size(G))

    BStor = newBoundaryData2D(G,P.order)
    DiffCoeff = [zeros(TT,size(G)),zeros(TT,size(G))]
    setCoefficient!(P.Kx,DiffCoeff[1],G)
    setCoefficient!(P.Ky,DiffCoeff[2],G)
    D = DerivativeOperator{TT,2,true,false,false}(P.order,G.nx,G.ny,G.Δx,G.Δy)

    SC = StepConfig{TT}(TT(0),TT(0),TT(0))

    LocalDataBlock{TT,2,typeof(u),typeof(DiffCoeff),typeof(BStor),typeof(D)}(u,uₙ₊₁,DiffCoeff,BStor,D,SC)
end



mutable struct newLocalDataBlock{TT<:Real,
        DIM,
        AT  <: AbstractArray{TT},
        KT  <: Union{Vector{TT},Vector{Matrix{TT}}},
        BT  <: BoundaryStorage{TT,DIM,AT},
        DT  <: DerivativeOperator,
        } <: newLocalDataBlockType{TT,DIM}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT
    boundary    :: BT
    Derivative  :: DT

    innerprod :: innerH{TT,DIM,Vector{TT}}
    cache     :: AT
    rₖ        :: AT
    dₖ        :: AT
    b         :: AT

    SC          :: StepConfig{TT}

    function newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}

        u       = zeros(TT,size(G))
        uₙ₊₁    = similar(u)
        cache   = similar(u)
        rₖ      = similar(u)
        dₖ      = similar(u)
        b       = similar(u)

        BStor = newBoundaryData1D(G,P.order,P.BoundaryConditions.BoundaryLeft,P.BoundaryConditions.BoundaryRight)
        DiffCoeff   = zeros(TT,size(G))

        IP = innerH(G,GetOrder(P.order))
        
        D = DerivativeOperator{TT,1,true,false,false}(P.order,G.n,0,G.Δx,TT(0))

        SC = StepConfig{TT}(TT(0),TT(0),TT(0))

        new{TT,1,typeof(u),typeof(DiffCoeff),typeof(BStor),typeof(D)}(u,uₙ₊₁,DiffCoeff,BStor,D,IP,cache,rₖ,dₖ,b,SC)
    end
end


"""
    DataMultiBlock
Data structure for multiblock problems
"""
struct DataMultiBlock{TT<:Real,
        DIM,
        TDBLOCK} <: newDataBlockType{TT,DIM}
        
    Block   :: TDBLOCK
    SC      :: StepConfig{TT}
    Joints  :: Nothing
    nblock  :: Int64


    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = [newLocalDataBlock(P,G)]
        new{TT,DIM,typeof(DTA)}(DTA,StepConfig{TT}(t,Δt,TT(0)),nothing,length(DTA))
    end
    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
        DTA = (LocalDataBlock(P,Grid,Δt) for Grid in G.Grids)
        new{TT,DIM,typeof(DTA)}(DTA,t,Δt,nothing,0.0)
    end
end



Base.getindex(DB::DataMultiBlock,i::Integer) = DB.Block[i]
Base.length(DB::DataMultiBlock) = DB.nblock
@inline eachblock(DB::DataMultiBlock) = Base.OneTo(length(DB))







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
function CommBoundary(B::BoundaryStorage,DB::DataMultiBlock{TT,DIM,NBLOCK}) where {TT,DIM,NBLOCK}
end
function CommBoundary(dest::AT,source::AT,inds::Tuple) where {AT}
    dest .= source[inds]
end
function CommBoundary(D1::LocalDataBlock,D2::LocalDataBlock,N::NodeType)
    D1.boundary .= D2.u()
    D2.boundary .= D1.u()
end