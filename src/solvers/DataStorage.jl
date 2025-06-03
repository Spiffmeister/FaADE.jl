



"""
    StepConfig{TT}
Configuration and stats for the current step for the current block
"""
mutable struct StepConfig{TT}
    t   :: TT
    Δt  :: TT
    Δu  :: TT
    θ   :: TT
    converged :: Bool

    function StepConfig{TT}(t=TT(0),Δt=TT(0),θ=TT(1)) where TT
        new{TT}(t,Δt,TT(0),θ,true)
    end
end

#========== NEW  DATA ==============#

function _setKoefficient!(K,P::Problem2D,G::LocalGridType{TT,2,CurvilinearMetric},Para::PT) where {TT,PT}
    # Used to compute I-BB^T/|B|^2
    function MagField(X)
        if (PT <: Nothing)
            return zeros(TT,3)
        else
            return Para.MagneticField.B(X,TT(0.0))
        end
    end

    cache = zeros(TT,size(G)) # Temporary storage

    if typeof(P.Kx) <: Real
        for i in eachindex(G)
            B = MagField(G[I])
            NB = norm(B,2)^2
            if NB == TT(0)
                NB = TT(1)
            end
            K[1][i] = P.Kx * (TT(1) - B[1]^2/NB) * G.J[i] * (G.qx[i]^2 + G.qy[i]^2)
            K[2][i] = P.Ky * (TT(1) - B[2]^2/NB) * G.J[i] * (G.rx[i]^2 + G.ry[i]^2)

            K[3][i] = P.Kx * B[1]*B[2]/NB * G.J[i] * (G.qx[i]*G.rx[i] + G.qy[i]*G.ry[i])
        end
    elseif typeof(P.Kx) <: Function

        # D₁!(K[3],K[1],G.nx,G.Δx,2,TT(0),1) # Cross derivative term in x
        # D₁!(K[4],K[2],G.ny,G.Δy,2,TT(0),2) # Cross derivative term in y

        for i in eachindex(G)
            B = MagField(G[i])
            NB = norm(B,2)^2
            if NB == TT(0)
                NB = TT(1)
            end
            K[1][i] = P.Kx(G[i]) * J[i] * (G.qx[i]^2 * (TT(1) - B[1]^2/NB)
                    + G.qy[i]^2 * (TT(1) - B[2]^2/NB) + 2*G.qx[i]*G.qy[i]*B[1]*B[2]/NB)
            K[2][i] = P.Ky(G[i]) * J[i] * (G.rx[i]^2 * (TT(1) - B[1]^2/NB)
                    + G.ry[i]^2 * (TT(1) - B[2]^2/NB) + 2*G.rx[i]*G.ry[i]*B[1]*B[2]/NB)
            K[3][i] = P.Kx(G[i]) * G.J[i] * (G.qx[i]*G.rx[i] * (TT(1) - B[1]^2/NB)
                    + G.qy[i]*G.ry[i] * (TT(1) - B[2]^2/NB) + (G.qx[i]*G.ry[i] + G.qy[i]*G.rx[i])*B[1]*B[2]/NB)
        end
    end
    K
end

function _setKoefficient!(K,P::Problem2D,G::LocalGridType{TT,2,CartesianMetric},Para::PT) where {TT,PT}
    function MagField(X)
        if PT <: Nothing
            return zeros(TT,3)
        else
            return Para.MagneticField.B(X,TT(0.0))
        end
    end
    if typeof(P.Kx) <: Real
        for I in eachindex(G)
            B = MagField(G[I])
            NB = norm(B,2)^2
            # if NB ≤ TT(2e-16)
            # if NB ≤ 1e-14
            if NB == 0
                NB = TT(1)
                B[1] = B[2] = TT(0)
                # println(I)
            end
            K[1][I] = P.Kx * (TT(1) - B[1]^2/NB)
            K[2][I] = P.Ky * (TT(1) - B[2]^2/NB)
            if !(PT<:Nothing)
                K[3][I] = -P.Kx * B[1]*B[2]/NB
            end
            # if isnan(K[1][I]) || isnan(K[2][I]) || isnan(K[3][I])
            #     println("NaN ",I," ",B[1]," ",B[2]," ",NB)
            # end
            x,y = G[I]
            # if (abs(x) == 0.5) && (abs(y) == 0.5) #NIMROD TEST
                # K[1][I] = P.Kx * 0.5
                # K[2][I] = P.Ky * 0.5
                # if !(PT<:Nothing)
                #     K[3][I] = -P.Kx*0.5
                # end
            # end
            if (x==0) && (y==0) #SINGLE ISLAND TEST
                K[1][I] = P.Kx
                K[2][I] = P.Ky
                if !(PT<:Nothing)
                    K[3][I] = 0.0
                end
            end
        end
        # for r in eachrow(K[1])
        #     println(r)
        # end
    elseif typeof(P.Kx) <: Function
        tmp = zeros(eltype(G),size(G))
        for I in eachindex(G)
            B = MagField(G[I])
            NB = norm(B,2)^2
            if NB == TT(0)
                NB = TT(1)
            end
            K[1][I] = P.Kx(G[I]) * (TT(1) - B[1]^2/NB)
            K[2][I] = P.Ky(G[I]) * (TT(1) - B[1]^2/NB)
            K[3][I] = -P.Ky(G[I]) * B[1]*B[2]/NB
        end
    end
    K
end




function _BuildGenericLocalDataBlocks(G::LocalGridType{TT}) where {TT}
    u = zeros(TT,size(G))
    uₙ₊₁ = zeros(TT,size(G))
    cache = zeros(TT,size(G))
    rₖ = zeros(TT,size(G))
    dₖ = zeros(TT,size(G))
    b = zeros(TT,size(G))

    return u, uₙ₊₁, cache, rₖ, dₖ, b
end

function _BuildDiffusionMatrix end
function _BuildDiffusionMatrix(G::LocalGridType{TT,1,MET},Para::PT) where {TT,MET,PT}
    K = zeros(TT,size(G))
    return K
end
function _BuildDiffusionMatrix(G::LocalGridType{TT,2,CurvilinearMetric},P,Para::PT) where {TT,PT<:ParallelData}

    existb = true
    try
        existb = !(typeof(Para.MagneticField).parameters[1] == Nothing)
    catch
        existb = false
    end

    function MagField(X)
        if existb
            return Para.MagneticField.B(X,TT(0))
        else
            return zeros(TT,3)
        end
    end

    K = [zeros(TT,size(G)), zeros(TT,size(G)), zeros(TT,size(G))]

    if typeof(P.Kx) <: Real
        for i in eachindex(G)
            if existb
                B = MagField(G[i])
                NB = norm(B,2)^2
                if NB == TT(0) # Ensure there are no divisions by zero
                    @warn "The field is singular at $(G[i])"
                    NB = TT(1)
                end
                Kx = P.Kx * (TT(1) - B[1]^2/NB)
                Ky = P.Ky * (TT(1) - B[2]^2/NB)
                Kxy = -P.Kx * B[1]*B[2]/NB
            else
                Kx = P.Kx
                Ky = P.Ky
                Kxy = 0.0
            end
            K[1][i] = Kx * G.J[i] * (G.qx[i]^2 + G.qy[i]^2)
            K[2][i] = Ky * G.J[i] * (G.rx[i]^2 + G.ry[i]^2)

            K[3][i] = Kxy * G.J[i] * (G.qx[i]*G.rx[i] + G.qy[i]*G.ry[i])
        end
    end

    return K
end
function _BuildDiffusionMatrix(G::LocalGridType{TT,2,CartesianMetric},P,Para::PT) where {TT,PT<:ParallelData}

    existb = true
    try
        existb = !(typeof(Para.MagneticField).parameters[1] == Nothing)
    catch
        existb = false
    end

    function MagField(X)
        if existb
            return Para.MagneticField.B(X,TT(0))
        else
            return zeros(TT,3)
        end
    end

    if existb
        K = [zeros(TT,size(G)), zeros(TT,size(G)), zeros(TT,size(G))]
    else
        K = [zeros(TT,size(G)), zeros(TT,size(G))]
    end
    # K = [zeros(TT,size(G)), zeros(TT,size(G)), zeros(TT,size(G))]

    if typeof(P.Kx) <: Real
        for i in eachindex(G)
            B = MagField(G[i])
            NB = norm(B,2)^2
            if NB == TT(0) # Ensure there are no divisions by zero
                NB = TT(1)
                B[1] = B[2] = TT(0)
                if existb
                    # @warn "norm(B)=0 at $(G[i]), assuming isotropic diffusion at this point. Check _BuildDiffusionMatrix in DataStorage.jl for more info."
                end
            end

            K[1][i] = P.Kx * (TT(1) - B[1]^2/NB)
            K[2][i] = P.Ky * (TT(1) - B[2]^2/NB)
            if existb
                K[3][i] = -P.Kx * B[1]*B[2]/NB
            end

        end
    elseif typeof(P.Kx) <: Function
        error("Needs fixing")
    end


    return K
end
function _BuildDiffusionMatrix(G::LocalGridType{TT,2,CartesianMetric},P,Para::Nothing) where {TT}
    K = [zeros(TT,size(G)), zeros(TT,size(G))]

    if typeof(P.Kx) <: Real
        for i in eachindex(G)
            K[1][i] = P.Kx
            K[2][i] = P.Ky
        end
    elseif typeof(P.Kx) <: Function
        error("Needs fixing")
    end
    return K
end
function _BuildDiffusionMatrix(G::LocalGridType{TT,2,CurvilinearMetric},P,Para::Nothing) where {TT}
    K = [zeros(TT,size(G)), zeros(TT,size(G)),zeros(TT,size(G))]

    if typeof(P.Kx) <: Real
        for i in eachindex(G)
            K[1][i] = P.Kx * G.J[i] * (G.qx[i]^2 + G.qy[i]^2)
            K[2][i] = P.Ky * G.J[i] * (G.rx[i]^2 + G.ry[i]^2)
            K[3][i] = P.Kx * G.J[i] * (G.qx[i]*G.rx[i] + G.qy[i]*G.ry[i])
        end
    elseif typeof(P.Kx) <: Function
        error("Needs fixing")
    end
    return K
end



"""
    _LocalDataBlockBlocks build all necessary blocks for a local data block
"""
function _LocalDataBlockBlocks(G::LocalGridType{TT,DIM,MET},Para::PT) where {TT,DIM,MET,PT}
    u       = zeros(TT,size(G))
    uₙ₊₁    = zeros(TT,size(G))
    if DIM == 1
        K       = zeros(TT,size(G))
        # K = (zeros(TT,G),)
    elseif DIM == 2
        existb = true
        try
            existb = !(typeof(Para.MagneticField).parameters[1] == Nothing)
        catch
            existb = false
        end

        if (MET == CurvilinearMetric) || existb
            K       = [zeros(TT,size(G)), zeros(TT,size(G)), zeros(TT,size(G))]
        else
            K       = [zeros(TT,size(G)), zeros(TT,size(G))]
        end
        # K = (zeros(TT,size(G)),zeros(TT,size(G)))
    end
    cache   = zeros(TT,size(G))
    rₖ      = zeros(TT,size(G))
    dₖ      = zeros(TT,size(G))
    b       = zeros(TT,size(G))

    return u, uₙ₊₁, K , cache, rₖ, dₖ, b
end
"""
    LocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT}
Contains the data for a local block
"""
mutable struct LocalDataBlock{TT<:Real,
        DIM,    # Dimension
        COORD,  # Coordinate type
        AT  <: AbstractArray{TT},
        KT,     # Coefficient storage
        DCT,    # Diffusion coefficient type
        GT <: GridType, # Grid type
        BT,     # Boundary conditions
        DT  <: DerivativeOperatorType,  # Derivative operator
        ST  <: SourceTerm,          # Source term
        PT, # Parallel map or Nothing
        } <: LocalDataBlockType{TT,DIM,AT}
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
    LocalDataBlock(P::PDEProblem{TT,1},G::LocalGridType) where {TT}
Initialise a data block for a 1D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function LocalDataBlock(P::PDEProblem{TT,1},G::LocalGridType,SC::StepConfig) where {TT}

    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(G)

    K = zeros(TT,size(G))
    PK = P.K
    if typeof(PK) <: Real
        K .= P.K
    elseif typeof(PK) <: Function
        for i in eachindex(G)
            K[i] = P.K(G[i])
        end
    end

    BS = GenerateBoundaries(P,G)

    IP = innerH(G.Δx,G.n,P.order)
    D = DiffusionOperator(G.n,G.Δx,P.order,false,:Constant)

    PMap = P.Parallel
    source = P.source

    return LocalDataBlock{TT,1,:Constant,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end
"""
    LocalDataBlock(P::PDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
Initialise a data block for a 1D multiblock problem
"""
function LocalDataBlock(P::PDEProblem{TT,1},G::GridMultiBlock,I::Integer,SC::StepConfig) where {TT}
    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(G.Grids[I])

    K = zeros(TT,size(G.Grids[I]))

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

    BS = GenerateBoundaries(P,G,I,K)
    IP = innerH(G.Grids[I].Δx,G.Grids[I].n,P.order)
    # D = DerivativeOperator{TT,1,typeof(P.order),:Constant}(P.order,G.Grids[I].n,0,G.Grids[I].Δx,TT(0),false,false)
    D = DiffusionOperator(G.Grids[I].n,G.Grids[I].Δx,P.order,false,:Constant)
    source = P.source

    return LocalDataBlock{TT,1,:Constant,typeof(u),typeof(K),typeof(P.K),typeof(G.Grids[I]),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K, P.K, G.Grids[I], BS, D,source,PMap, IP, cache,rₖ,dₖ,b,SC)
end
"""
    LocalDataBlock(P::PDEProblem{TT,2},G::LocalGridType) where {TT}
Initialise a data block for a 2D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function LocalDataBlock(P::PDEProblem{TT,2},G::LocalGridType,SC::StepConfig) where {TT}

    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(G)

    K = _BuildDiffusionMatrix(G,P,P.Parallel)

    PK = (P.Kx,P.Ky)

    BS = GenerateBoundaries(P,G,K)
    IP = innerH(G.Δx,G.Δy,G.nx,G.ny,P.order)

    # if length(K) == 3
    if GetMetricType(G) == CurvilinearMetric
        difftype = :Variable
    else
        difftype = :Constant
    end

    typeof(P.BoundaryConditions[1]).parameters[2] == :Cartesian ? sattype = :Constant : sattype = :Variable
    # sattype = :Variable
    # D = DerivativeOperator{TT,2,typeof(P.order),:Constant}(P.order,G.nx,G.ny,G.Δx,G.Δy,false,false)
    Dx = DiffusionOperator(G.nx,G.Δx,P.order,false,difftype)
    Dy = DiffusionOperator(G.ny,G.Δy,P.order,false,difftype)
    D = DiffusionOperatorND(Dx,Dy)
    PMap = P.Parallel
    source = P.source

    return LocalDataBlock{TT,2,sattype,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end
"""
    LocalDataBlock(P::PDEProblem{TT,2},G::GridMultiBlock,I::Integer)
Initialise a data block for a 2D multiblock problem
"""
function LocalDataBlock(P::PDEProblem{TT,2},G::GridMultiBlock{TT,2,MET},I::Integer,SC::StepConfig) where {TT,MET}
    LG = G.Grids[I]

    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(LG)

    PK = (P.Kx,P.Ky) #need to fix this for material properties

    if typeof(P.Parallel) <: ParallelMultiBlock
        PMap = P.Parallel.PData[I]
    else
        PMap = P.Parallel
    end

    K = _BuildDiffusionMatrix(LG,P,PMap)

    BS = GenerateBoundaries(P,G,I,K)
    IP = innerH(LG.Δx,LG.Δy,LG.nx,LG.ny,P.order)

    if GetMetricType(G) == CurvilinearMetric
        difftype = :Variable
    else
        difftype = :Constant
    end

    typeof(BS[Left].BoundaryOperator).parameters[2] == :Cartesian ? sattype = :Constant : sattype = :Variable
    # sattype = :Constant
    # @show sattype, I, typeof(BS[Left].BoundaryOperator)
    # @show typeof(BS[Left].BoundaryOperator).parameters
    for (side,boundary) in BS
        # @show typeof(boundary)
        # @show :Variable ∈ typeof(BC.Boundary).parameters
        if :Variable ∈ typeof(boundary.BoundaryOperator).parameters
            sattype = :Variable
        end
    end

    # @show typeof(BS[1].Boundary)
    # @show sattype, I
    Dx = DiffusionOperator(LG.nx,LG.Δx,P.order,false,difftype)
    Dy = DiffusionOperator(LG.ny,LG.Δy,P.order,false,difftype)
    D = DiffusionOperatorND(Dx,Dy)

    if typeof(P.source.source) <: Vector
        source = SourceTerm{typeof(P.source.source[I])}(P.source.source[I]) #TODO: Hacked in, fix
    else
        source = P.source
    end

    return LocalDataBlock{TT,2,sattype,typeof(u),typeof(K),typeof(PK),typeof(LG),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,LG,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end


"""
    getproperty(D::LocalDataBlock,s::Symbol)
"""
@inline function Base.getproperty(D::LocalDataBlock,s::Symbol)
    rt = getfield(D,s)
    return rt :: typeof(rt)
end
"""
    getarray(D::LocalDataBlock,s::Symbol)
Type stable getfield for arrays
"""
@inline function getarray(D::LocalDataBlockType{TT,DIM,AT},s::Symbol) where {TT,DIM,AT}
    rt = getfield(D,s)
    return rt :: AT
end

"""
    DataMultiBlock
Data structure for multiblock problems
"""
struct DataMultiBlock{TT<:Real,
        DIM,
        NB,
        TDBLOCK <: NTuple{NB,LocalDataBlockType},
        PTBLOCK <: Union{Nothing,ParallelData,ParallelMultiBlock}
            } <: DataBlockType{TT,DIM}

    Block   :: TDBLOCK
    ParallelData    :: PTBLOCK
    SC      :: StepConfig{TT}
    nblock  :: Int64
    parallel:: Bool

    function DataMultiBlock(P::PDEProblem{TT,DIM},G::LocalGridType{TT},Δt::TT,t::TT;θ=TT(1)) where {TT,DIM}
        SC = StepConfig{TT}(t,Δt,θ)
        DTA = (LocalDataBlock(P,G,SC),)
        new{TT,DIM,1,typeof(DTA),typeof(P.Parallel)}(DTA,P.Parallel,SC,length(DTA),false)
    end
    function DataMultiBlock(P::PDEProblem{TT,DIM},G::GridMultiBlock{TT,DIM},Δt::TT,t::TT;θ=TT(1)) where {TT,DIM}
        SC = StepConfig{TT}(t,Δt,θ)
        DTA = map((x)->LocalDataBlock(P,G,x,SC),eachgrid(G))
        DTA = tuple(DTA...)
        new{TT,DIM,length(DTA),typeof(DTA),typeof(P.Parallel)}(DTA,P.Parallel,SC,length(DTA),false)
    end
    function DataMultiBlock(D::LocalDataBlock{TT,DIM}) where {TT,DIM} #== TESTING METHOD ==#
        DTA = (D,)
        new{TT,DIM,1,typeof(DTA),typeof(P.Parallel)}(DTA,P.Parallel,D.SC,length(DTA),false)
    end
end




function Base.iterate(DMB::DataMultiBlock,state=1)
    state >= length(DMB) && return
    return DMB.Block[state], state+1
end

Base.getindex(DB::DataMultiBlock,i::Integer) = Base.getindex(DB.Block,i)
Base.length(DB::DataMultiBlock) = DB.nblock
@inline eachblock(DB::DataMultiBlock) = Base.OneTo(DB.nblock)

Base.ndims(DB::DataMultiBlock{TT,DIM}) where {TT,DIM} = DIM

# @inline eachjoint(BB::BoundaryData1D) = Base.OneTo(2)



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
function collectΔu(D::LocalDataBlock) end
