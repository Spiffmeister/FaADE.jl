
abstract type newDataBlockType{dtype,DIM} end
abstract type newLocalDataBlockType{dtype,DIM} <: newDataBlockType{dtype,DIM} end



"""
    StepConfig{TT}
Configuration and stats for the current step for the current block
"""
mutable struct StepConfig{TT}
    t   :: TT
    Δt  :: TT
    Δu  :: TT
    converged :: Bool

    function StepConfig{TT}(t=TT(0),Δt=TT(0)) where TT
        new{TT}(t,Δt,TT(0),true)
    end
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
        CLR <: CartesianIndices,
        AT} <: BoundaryStorage{TT,1,AT}

    Left        :: SATL
    Right       :: SATR

    JointLeft   :: Int64
    JointRight  :: Int64

    LeftIndexFrom   :: CLR
    RightIndexFrom  :: CLR
    LRStorageIndex  :: CLR

    u_Left      :: AT
    u_Right     :: AT

    RHS_Left    :: AT
    RHS_Right   :: AT

    function newBoundaryData1D(G::Grid1D{TT,MET},order::DerivativeOrder{O},BCL,BCR) where {TT,O,MET}

        nnodes = SATNodeOutput(O)

        u_Left      = zeros(TT,nnodes)
        u_Right     = zeros(TT,nnodes)

        LeftIndex   = CartesianIndices((1:GetOrder(order),))
        RightIndex  = CartesianIndices((lastindex(G)-O+1:lastindex(G),))

        new{TT,typeof(BCL),typeof(BCR),typeof(LeftIndex),typeof(u_Left)}(BCL,BCR,1,1,LeftIndex,RightIndex,LeftIndex,u_Left,u_Right,[0.0],[0.0])
    end
    function newBoundaryData1D(G::GridMultiBlock{TT,DIM},order::DerivativeOrder{O},SATB,i::Int64) where {TT,DIM,O}

        nnodes = SATNodeOutput(O)

        u_Left = zeros(TT,nnodes)
        u_Right = zeros(TT,nnodes)

        LeftIndex = CartesianIndices((1:O,))
        RightIndex = CartesianIndices((lastindex(G[i])-O+1:lastindex(G[i]),))

        J = filter(t -> (t[1] == i || t[2] == i), G.Joint)

        # Decide which indicies need to be written to
        J = filter(t -> t[1] == i, G.Joint)
        println(J)
        if isempty(J)
            JL = i
            IL = LeftIndex
        else
            J = J[1]
            JL = J[2]
            if J[3] == Left
                IL = RightIndex
            elseif J[3] == Right
                IL = LeftIndex
            end
            # SATL = SATB.BoundaryLeft
        end
        J = filter(t -> t[2] == i, G.Joint)
        if isempty(J)
            JR = i
            IR = RightIndex
        else
            J = J[1]
            JR = J[1]
            if J[3] == Left
                IR = RightIndex
            elseif J[3] == Right
                IR = LeftIndex
            end
            # SATR = SATB.BoundaryRight
        end
        SATL = SATB.BoundaryLeft
        SATR = SATB.BoundaryRight

        new{TT,typeof(SATL),typeof(SATR),typeof(LeftIndex),typeof(u_Left)}(SATL,SATR,JL,JR,IL,IR,LeftIndex,u_Left,u_Right,[0.0],[0.0])
    end
end


"""
    BoundaryData2D
Data structure for storage of SATs in 2 dimensional problems
"""
struct newBoundaryData2D{TT,
        SATL <: SimultanousApproximationTerm,
        SATR <: SimultanousApproximationTerm,
        SATD <: SimultanousApproximationTerm,
        SATU <: SimultanousApproximationTerm,
        AT} <: BoundaryStorage{TT,2,AT}

    Left        :: SATL
    Right       :: SATR
    Up          :: SATU
    Down        :: SATD

    JointLeft   :: Int64
    JointRight  :: Int64
    JointUp     :: Int64
    JointDown   :: Int64

    LeftIndex   :: CartesianIndices
    RightIndex  :: CartesianIndices
    UpIndex     :: CartesianIndices
    DownIndex   :: CartesianIndices

    u_Left      :: AT #Solution along boundary with size determined by derivative order
    u_Right     :: AT #Solution along boundary with size determined by derivative order
    u_Up        :: AT #Solution along boundary with size determined by derivative order
    u_Down      :: AT #Solution along boundary with size determined by derivative order

    RHS_Left    :: AT
    RHS_Right   :: AT
    RHS_Up      :: AT
    RHS_Down    :: AT

    function newBoundaryData2D(grid::Grid2D{TT,MET},order::DerivativeOrder{O},BCL,BCR,BCU,BCD) where {TT,O,MET}

        nnodes = SATNodeOutput(order)
        ny = grid.ny
        nx = grid.nx

        u_Left =    zeros(TT,(nnodes,ny)) #x
        u_Right =   zeros(TT,(nnodes,ny)) #x
        u_Up =      zeros(TT,(nx,nnodes)) #y
        u_Down =    zeros(TT,(nx,nnodes)) #y

        RHS_Left =    zeros(TT,(1,ny)) #x
        RHS_Right =   zeros(TT,(1,ny)) #x
        RHS_Up =      zeros(TT,(nx,1)) #y
        RHS_Down =    zeros(TT,(nx,1)) #y

        LeftIndex   = CartesianIndices((1:O,1:ny))
        RightIndex  = CartesianIndices((lastindex(G)-O+1:lastindex(G),1:ny))
        UpIndex     = CartesianIndices((1:nx,1:O))
        DownIndex   = CartesianIndices((1:nx,lastindex(G)-O+1:lastindex(G)))

        new{TT,typeof(BCL),typeof(BCR),typeof(BCU),typeof(BCD),typeof(u_Left)}(BCL,BCR,BCU,BCD,
            1,1,1,1,
            LeftIndex,RightIndex,UpIndex,DownIndex,
            u_Left,u_Right,u_Up,u_Down,
            RHS_Left,RHS_Right,RHS_Up,RHS_Down)
    end
end



function _newLocalDataBlockBlocks(G::LocalGridType{TT}) where {TT}
    u       = zeros(TT,size(G))
    uₙ₊₁    = zeros(TT,size(G))
    K       = zeros(TT,size(G))
    cache   = zeros(TT,size(G))
    rₖ      = zeros(TT,size(G))
    dₖ      = zeros(TT,size(G))
    b       = zeros(TT,size(G))

    return u, uₙ₊₁, K , cache, rₖ, dₖ, b
end
mutable struct newLocalDataBlock{TT<:Real,
        DIM,
        AT  <: AbstractArray{TT},
        KT,
        BT  <: BoundaryStorage{TT,DIM,AT},
        DT  <: DerivativeOperator,
        } <: newLocalDataBlockType{TT,DIM}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: AT

    κ          :: KT

    boundary    :: BT
    Derivative  :: DT

    innerprod :: innerH{TT,DIM,Vector{TT}}
    cache     :: AT
    rₖ        :: AT
    dₖ        :: AT
    b         :: AT

    SC          :: StepConfig{TT}

    function newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}

        u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(G)

        BStor = newBoundaryData1D(G,P.order,P.BoundaryConditions.BoundaryLeft,P.BoundaryConditions.BoundaryRight)

        IP = innerH(G,GetOrder(P.order))
        
        D = DerivativeOperator{TT,1,typeof(P.order),true,false,false}(P.order,G.n,0,G.Δx,TT(0))

        SC = StepConfig{TT}()

        new{TT,1,typeof(u),typeof(P.K),typeof(BStor),typeof(D)}(u,uₙ₊₁,K,P.K,BStor,D,IP,cache,rₖ,dₖ,b,SC)
    end
    function newLocalDataBlock(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
        u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(G.Grids[I])

        BStor = newBoundaryData1D(G,P.order,P.BoundaryConditions,I)

        IP = innerH(G.Grids[I],GetOrder(P.order))

        D = DerivativeOperator{TT,1,typeof(P.order),true,false,false}(P.order,G.Grids[I].n,0,G.Grids[I].Δx,TT(0))

        SC = StepConfig{TT}()

        new{TT,1,typeof(u),typeof(P.K),typeof(BStor),typeof(D)}(u,uₙ₊₁,K, P.K, BStor, D, IP, cache,rₖ,dₖ,b,SC)

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


    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT},Δt::TT,t::TT) where {TT,DIM}
        DTA = [newLocalDataBlock(P,G)]
        SC = StepConfig{TT}(t,Δt)
        new{TT,DIM,typeof(DTA)}(DTA,SC,nothing,length(DTA))
    end
    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridMultiBlock{TT,DIM},Δt::TT,t::TT) where {TT,DIM}
        DTA = [newLocalDataBlock(P,G,I) for I in eachgrid(G)]
        SC = StepConfig{TT}(t,Δt)
        new{TT,DIM,typeof(DTA)}(DTA,SC,nothing,length(DTA))
    end
    # function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
    #     DTA = [newLocalDataBlock(P,Grid,Δt) for Grid in G.Grids]
    #     new{TT,DIM,typeof(DTA)}(DTA,StepConfig{TT}(t,Δt,TT(0),true),nothing,0.0)
    # end
end



Base.getindex(DB::DataMultiBlock,i::Integer) = DB.Block[i]
Base.length(DB::DataMultiBlock) = DB.nblock
@inline eachblock(DB::DataMultiBlock) = Base.OneTo(length(DB))

Base.ndims(DB::DataMultiBlock{TT,DIM}) where {TT,DIM} = DIM

# Base.getindex(BB::newBoundaryData1D,i::Integer) = BB.
@inline eachjoint(BB::newBoundaryData1D) = Base.OneTo(2)





