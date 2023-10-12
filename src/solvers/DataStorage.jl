
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

"""
    newBoundaryData
Container for Dirichlet, Neumann and Robin boundary conditions.
"""
struct newBoundaryData{
        TT<:Real,
        DIM,
        F1<:Union{Real,Function},
        BCT,
        AT} <: BoundaryStorage{TT,0,AT}

    Boundary    :: BCT
    RHS         :: F1
    BufferRHS   :: AT
    X           :: Vector{TT}   # Grid points along boundary
    n           :: Int64        # Length of boundary
    DIM         :: Int64

    function newBoundaryData(G::Grid1D{TT},BC,Joint,order::DerivativeOrder{O}) where {TT,O}

        BufferRHS = zeros(TT,1)

        new{TT,1,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,[TT(0)],1,1)
    end
    function newBoundaryData(G::Grid2D{TT},BC,Joint,order::DerivativeOrder{O}) where {TT,O}

        if BC.side ∈ [Left,Right]
            n = G.ny
            BufferRHS = zeros(TT,(1,n))
            if BC.side == Left
                X = ones(TT,n) * G.gridy[1]
            else
                X = ones(TT,n) * G.gridy[n]
            end
        elseif BC.side ∈ [Up,Down]
            n = G.nx
            BufferRHS = zeros(TT,(n,1))
            if BC.side == Up
                X = ones(TT,n) * G.gridx[1]
            else
                X = ones(TT,n) * G.gridx[n]
            end
        end

        new{TT,2,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,X,n,2)
    end
end


"""
    newInterfaceBoundaryData
Container for Periodic and Interface boundary conditions.
"""
struct newInterfaceBoundaryData{
        TT<:Real,
        DIM,
        BCT,
        AT} <: BoundaryStorage{TT,0,AT}

    Boundary    :: BCT
    BufferOut   :: AT
    BufferIn    :: AT
    Joint       :: Int64
    DIM         :: Int64
    I           :: CartesianIndices

    function newInterfaceBoundaryData{TT}(G1::Grid1D,G2::Grid1D,BC,Joint,order::DerivativeOrder{O}) where {TT,O}

        BufferIn    = zeros(TT,O)
        BufferOut   = zeros(TT,O)
        if BC.side == Right
            I = CartesianIndices((1:O,))
        elseif BC.side == Left
            I = CartesianIndices((length(G2)-O+1:length(G2),))
        end

        new{TT,1,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint,1,I)
    end
    function newInterfaceBoundaryData{TT}(G1::Grid2D,G2::Grid2D,BC,Joint,order::DerivativeOrder{O}) where {TT,O}

        if BC.side ∈ [Left,Right]
            n = G2.ny
            BufferIn    = zeros(TT,(O,n))
            BufferOut   = zeros(TT,(O,n))

            if BC.side == Left
                I = CartesianIndices((1:O,1:n))
            else
                I = CartesianIndices((G2.nx-O+1:G2.nx,1:n))
            end

        elseif BC.side ∈ [Up,Down]
            n = G.nx
            BufferIn    = zeros(TT,(n,O))
            BufferOut   = zeros(TT,(n,O))

            if BC.side == Up
                I = CartesianIndices((1:n,1:O))
            else
                I = CartesianIndices((1:n,G2.ny-O+1:G2.ny))
            end
        end

        new{TT,2,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint,2,I)
    end
end
"""
    newBoundaryConditions
Container for all boundary conditions
"""
struct newBoundaryConditions{DIM,
        BCLT,
        BCRT,
        BCUT,
        BCDT}
    
    BC_Left     :: BCLT
    BC_Right    :: BCRT
    BC_Up       :: BCUT
    BC_Down     :: BCDT

    # BC      :: NTuple{DIM-1,BoundaryStorage}

    """
        newBoundaryConditions(P::newPDEProblem{TT,1},G::LocalGridType)
    """
    function newBoundaryConditions(P::newPDEProblem{TT,DIM},G::LocalGridType) where {TT,DIM}
        BCL = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryLeft,1,P.order)
        BCR = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryRight,1,P.order)
        if DIM == 1
            # NB = 2
            BCU = nothing
            BCD = nothing
            # BC = (BCL,BCR)
        elseif DIM ==2
            # NB = 4
            BCU = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryUp,1,P.order)
            BCD = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryDown,1,P.order)
            # BC = (BCL,BCR,BCU,BCD)
        end
        new{1,typeof(BCL),typeof(BCR),typeof(BCU),typeof(BCD)}(BCL,BCR,BCU,BCD)
    end
    """
        newBoundaryConditions(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Int64)
    """
    function newBoundaryConditions(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Int64) where {TT}

        J = G.Joint[I]
        if length(J) == 1
            if J[1][2] == Left
                BCLt = SAT_Interface(G.Grids[J[1][1]].Δx,G.Grids[I].Δx,J[1][2],1,GetOrder(P.order))
# println("type",typeof(P.BoundaryConditions.BoundaryLeft))
                BCL = _newBoundaryCondition(G.Grids[I],G.Grids[J[1][1]],BCLt,J[1][1],P.order)
                BCR = _newBoundaryCondition(G.Grids[I],P.BoundaryConditions.BoundaryRight,J[1][1],P.order)
            elseif J[1][2] == Right
                BCRt = SAT_Interface(G.Grids[I].Δx,G.Grids[J[1][1]].Δx,J[1][2],1,GetOrder(P.order))
                
                BCL = _newBoundaryCondition(G.Grids[I],P.BoundaryConditions.BoundaryLeft,J[1][1],P.order)
                BCR = _newBoundaryCondition(G.Grids[I],G.Grids[J[1][1]],BCRt,J[1][1],P.order)
            end
        elseif length(J) == 2

            JL = G.Joint[J[1][1]]
            JR = G.Joint[J[2][1]][1]
            length(JL) == 1 ? JL = JL[1] : JL = JL[2] # Correct for second node where `JL = ((1,Right),)`

            println(J)
            println(JL)
            println(JR)

            BCLt = SAT_Interface(G.Grids[J[1][1]].Δx,G.Grids[JL[1]].Δx,Left,1,GetOrder(P.order))
            BCRt = SAT_Interface(G.Grids[J[2][1]].Δx,G.Grids[JR[1]].Δx,Right,1,GetOrder(P.order))

            BCL = _newBoundaryCondition(G.Grids[I],G.Grids[J[1][1]],BCLt,J[1][1],P.order)
            BCR = _newBoundaryCondition(G.Grids[I],G.Grids[J[2][1]],BCRt,J[2][1],P.order)
        end

        BCU = nothing
        BCD = nothing

        new{1,typeof(BCL),typeof(BCR),Nothing,Nothing}(BCL,BCR,BCU,BCD)
    end
    function newBoundaryConditions(P::newPDEProblem{TT,2},G::GridMultiBlock,I::Int64) where {TT}
        J = G.Joint[I]

        faces = [Left,Right,Up,Down]
        Pfaces = [:BoundaryLeft,:BoundaryRight,:BoundaryUp,:BoundaryDown]

        BCS = []
        if length(G.Joint[I]) != 4
            error("Insufficient information to set up block boundary conditions")
        end

        # Loop over joints
        for i in 1:length(G.Joint[I])
        # for J in G.Joint[I]
            J = G.Joint[I][i]
            if J[1] == I
                BCt = getfield(P.BoundaryConditions,Pfaces[i])
                BC = _newBoundaryCondition(G.Grids[I],BCt,J[1],P.order)
                push!(BCS,BC)
            else
                BCt = SAT_Interface(G.Grids[J[1]].Δx,G.Grids[I].Δx,J[2],GetAxis(J[2]),GetOrder(P.order))
                BC = _newBoundaryCondition(G.Grids[I],G.Grids[J[1]],BCt,J[1],P.order)
                push!(BCS,BC)
            end
        end
        new{2,typeof(BCS[1]),typeof(BCS[2]),typeof(BCS[3]),typeof(BCS[4])}(BCS[1],BCS[2],BCS[3],BCS[4])
    end
end


function getjoint(BC::newBoundaryData) end
getjoint(BC::newInterfaceBoundaryData) = BC.Joint



function Base.iterate(B::newBoundaryConditions{DIM,BCLT,BCRT,BCUT,BCDT},state=0) where {DIM,BCLT,BCRT,BCUT,BCDT}
    state >= 2DIM && return
    return Base.getfield(B,state+1), state+1
end
Base.getindex(B::newBoundaryConditions,N::NodeType{:Left})  = B.BC_Left
Base.getindex(B::newBoundaryConditions,N::NodeType{:Right}) = B.BC_Right
Base.getindex(B::newBoundaryConditions,N::NodeType{:Up})    = B.BC_Up
Base.getindex(B::newBoundaryConditions,N::NodeType{:Down})  = B.BC_Down

_newBoundaryCondition(G1::GridType{TT},G2::GridType{TT},BC::SimultanousApproximationTerm{:Periodic},J,order) where TT  = newInterfaceBoundaryData{TT}(G1,G2,BC,J,order)
_newBoundaryCondition(G1::GridType{TT},G2::GridType{TT},BC::SimultanousApproximationTerm{:Interface},J,order) where TT = newInterfaceBoundaryData{TT}(G1,G2,BC,J,order)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Dirichlet},J,order) where TT = newBoundaryData(G,BC,J,order)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Neumann},J,order) where TT   = newBoundaryData(G,BC,J,order)

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



function _newLocalDataBlockBlocks(G::LocalGridType{TT,DIM}) where {TT,DIM}
    u       = zeros(TT,size(G))
    uₙ₊₁    = zeros(TT,size(G))
    if DIM == 1
        K       = zeros(TT,size(G))
        # K = (zeros(TT,G),)
    elseif DIM == 2
        K       = [zeros(TT,size(G)), zeros(TT,size(G))]
        # K = (zeros(TT,size(G)),zeros(TT,size(G)))
    end
    cache   = zeros(TT,size(G))
    rₖ      = zeros(TT,size(G))
    dₖ      = zeros(TT,size(G))
    b       = zeros(TT,size(G))

    return u, uₙ₊₁, K , cache, rₖ, dₖ, b
end
"""
    newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT}
Contains the data for a local block
"""
mutable struct newLocalDataBlock{TT<:Real,
        DIM,
        AT  <: AbstractArray{TT},
        KT,
        DCT,
        GT <: GridType,
        BT,
        DT  <: DerivativeOperator,
        ST  <: SourceTerm,
        PT, # Parallel map or Nothing
        } <: newLocalDataBlockType{TT,DIM}
    u           :: AT
    uₙ₊₁        :: AT
    K           :: KT

    κ          :: DCT
    
    grid        :: GT
    
    boundary    :: BT

    Derivative  :: DT # :: DerivativeOperator{TT,DIM,DerivativeOrder}

    source      :: ST

    Parallel    :: PT

    innerprod :: innerH{TT,DIM,Vector{TT}}
    cache     :: AT
    rₖ        :: AT
    dₖ        :: AT
    b         :: AT

    SC          :: StepConfig{TT}
end
"""
    newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}
Initialise a data block for a 1D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}

    u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(G)

    if typeof(P.K) <: Real
        PK = P.K
        K .= P.K
    elseif typeof(P.K) <: Function
        PK = P.K
        for i in eachindex(G)
            K[i] .= P.K(G[i])
        end
    end

    BStor = newBoundaryConditions(P,G)
    BS = (BStor.BC_Left,BStor.BC_Right)
    IP = innerH(G,GetOrder(P.order))
    D = DerivativeOperator{TT,1,typeof(P.order)}(P.order,G.n,0,G.Δx,TT(0))
    PMap = P.Parallel
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,1,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end
"""
    newLocalDataBlock(P::newPDEProblem{TT,2},G::LocalGridType) where {TT}
Initialise a data block for a 2D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function newLocalDataBlock(P::newPDEProblem{TT,2},G::LocalGridType) where {TT}

    u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(G)

    if typeof(P.Kx) <: Real
        PK = [P.Kx,P.Ky]
        K[1] .= P.Kx
        K[2] .= P.Ky
    elseif typeof(P.Kx) <: Function
        PK = [P.Kx,P.Ky]
        for i in eachindex(G)
            K[1][i] .= P.Ky(G[i])
            K[2][i] .= P.Kx(G[i])
        end
    end

    BStor = newBoundaryConditions(P,G)
    BS = (BStor.BC_Left,BStor.BC_Right,BStor.BC_Up,BStor.BC_Down)
    IP = innerH(G,GetOrder(P.order))
    D = DerivativeOperator{TT,2,typeof(P.order)}(P.order,G.nx,G.ny,G.Δx,G.Δy)
    PMap = P.Parallel
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,2,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end
"""
    newLocalDataBlock(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
Initialise a data block for a 1D multiblock problem
"""
function newLocalDataBlock(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
    u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(G.Grids[I])

    if typeof(P.K) <: Real
        K .= P.K
    elseif typeof(P.K) <: Function
        for i in eachindex(G)
            K[i] .= P.K(G[i])
        end
    end
    
    if typeof(P.Parallel) <: Vector
        PMap = P.Parallel[I]
    else
        PMap = P.Parallel
    end

    BStor = newBoundaryConditions(P,G,I)
    BS = (BStor.BC_Left,BStor.BC_Right)
    IP = innerH(G.Grids[I],GetOrder(P.order))
    D = DerivativeOperator{TT,1,typeof(P.order)}(P.order,G.Grids[I].n,0,G.Grids[I].Δx,TT(0))
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,1,typeof(u),typeof(K),typeof(P.K),typeof(G.Grids[I]),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K, P.K, G.Grids[I], BS, D,source,PMap, IP, cache,rₖ,dₖ,b,SC)
end
"""
    newLocalDataBlock(P::newPDEProblem{TT,2},G::GridMultiBlock,I::Integer)
Initialise a data block for a 2D multiblock problem
"""
function newLocalDataBlock(P::newPDEProblem{TT,2},G::GridMultiBlock,I::Integer) where TT
    LG = G.Grids[I]

    u, uₙ₊₁, K , cache, rₖ, dₖ, b = _newLocalDataBlockBlocks(LG)
    
    if typeof(P.Kx) <: Real
        K[1] .= P.Kx
        K[2] .= P.Ky
    elseif typeof(P.K) <: Function
        for i in eachindex(G[I])
            K[1][i] .= P.Kx(G[I][i])
            K[2][i] .= P.Ky(G[I][i])
        end
    end
    PK = [P.Kx,P.Ky]

    if typeof(P.Parallel) <: Vector
        PMap = P.Parallel[I]
    else
        PMap = P.Parallel
    end

    BStor = newBoundaryConditions(P,G,I)
    BS = (BStor.BC_Left,BStor.BC_Right,BStor.BC_Up,BStor.BC_Down)
    IP = innerH(LG,GetOrder(P.order))
    D = DerivativeOperator{TT,2,typeof(P.order)}(P.order,LG.nx,LG.ny,LG.Δx,LG.Δy)
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,2,typeof(u),typeof(K),typeof(PK),typeof(LG),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,LG,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end


"""
    getproperty(D::newLocalDataBlock,s::Symbol)
"""
@inline function Base.getproperty(D::newLocalDataBlock,s::Symbol)
    rt = getfield(D,s)
    return rt :: typeof(rt)
end
"""
    getarray(D::newLocalDataBlock,s::Symbol)
Type stable getfield for arrays
"""
@inline function getarray(D::newLocalDataBlock{TT,DIM,AT},s::Symbol) where {TT,DIM,AT}
    rt = getfield(D,s)
    return rt :: AT
end

@inline eachboundary(D::newLocalDataBlock{TT,DIM}) where {TT,DIM} = Base.OneTo(2DIM)

"""
    DataMultiBlock
Data structure for multiblock problems
"""
struct DataMultiBlock{TT<:Real,
        DIM,
        NB,
        TDBLOCK <: NTuple{NB,newLocalDataBlockType}
            } <: newDataBlockType{TT,DIM}
    
    Block   :: TDBLOCK
    SC      :: StepConfig{TT}
    nblock  :: Int64
    parallel:: Bool


    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT},Δt::TT,t::TT) where {TT,DIM}
        # DTA = [newLocalDataBlock(P,G)]
        DTA = (newLocalDataBlock(P,G),)
        SC = StepConfig{TT}(t,Δt)
        new{TT,DIM,1,typeof(DTA)}(DTA,SC,length(DTA),false)
    end
    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridMultiBlock{TT,DIM},Δt::TT,t::TT) where {TT,DIM}
        # DTA = [newLocalDataBlock(P,G,I)]
        # DTA = []
        # for I in eachgrid(G)
        #     push!(DTA, newLocalDataBlock(P,G,I))
        # end
        DTA = map((x)->newLocalDataBlock(P,G,x),eachgrid(G))
        DTA = tuple(DTA...)
        SC = StepConfig{TT}(t,Δt)
        new{TT,DIM,length(DTA),typeof(DTA)}(DTA,SC,length(DTA),false)
    end
    # function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM,MET},Δt::TT,t::TT) where {TT,DIM,MET}
    #     DTA = [newLocalDataBlock(P,Grid,Δt) for Grid in G.Grids]
    #     new{TT,DIM,typeof(DTA)}(DTA,StepConfig{TT}(t,Δt,TT(0),true),nothing,0.0)
    # end
end





function Base.iterate(DMB::DataMultiBlock,state=0)
    state >= length(DMB) && return
    return DMB.DTA[state], state+1
end

Base.getindex(DB::DataMultiBlock,i::Integer) = Base.getindex(DB.Block,i)
Base.length(DB::DataMultiBlock) = DB.nblock
@inline eachblock(DB::DataMultiBlock) = Base.OneTo(DB.nblock)

Base.ndims(DB::DataMultiBlock{TT,DIM}) where {TT,DIM} = DIM

@inline eachjoint(BB::newBoundaryData1D) = Base.OneTo(2)



"""
    collectΔu(DB::DataMultiBlock)
Collect the change in solution from each block and store it in the main block
"""
function collectΔu(DB::DataMultiBlock{TT}) where TT
    DB.SC.Δu = TT(0)
    for i = 1:length(DB)
        DB.SC.Δu += DB[i].SC.Δu
    end
end
function collectΔu(D::newLocalDataBlock) end
