
abstract type DataBlockType{dtype,DIM} end
abstract type LocalDataBlockType{dtype,DIM,atype} <: DataBlockType{dtype,DIM} end



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

struct BoundaryNull <: BoundaryStorage{Float64,0,Vector{Float64}} end

"""
    BoundaryData
Container for Dirichlet, Neumann and Robin boundary conditions.
"""
struct BoundaryData{
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

    function BoundaryData(G::Grid1D{TT},BC,Joint,order::Int) where {TT}

        BufferRHS = zeros(TT,1)

        new{TT,1,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,[TT(0)],1,1)
    end
    function BoundaryData(G::Grid2D{TT},BC,Joint,order::Int) where {TT}

        if BC.side ∈ [Left,Right]
            n = G.ny
            BufferRHS = zeros(TT,(1,n))
            # @show size(BufferRHS), size(G.gridx)
            if BC.side == Left
                X = G.gridy[1,:]
            else
                X = G.gridy[G.nx,:]
            end
        elseif BC.side ∈ [Up,Down]
            n = G.nx
            BufferRHS = zeros(TT,(n,1))
            if BC.side == Up
                X = G.gridx[:,1]
            else
                X = G.gridx[:,G.ny]
            end
        end

        new{TT,2,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,X,n,2)
    end
end


"""
    InterfaceBoundaryData
Container for Periodic and Interface boundary conditions.
"""
struct InterfaceBoundaryData{
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

    function InterfaceBoundaryData{TT}(G1::Grid1D,G2::Grid1D,BC,Joint,order::Int) where {TT}

        BufferIn    = zeros(TT,order)
        BufferOut   = zeros(TT,order)
        if BC.side == Right
            I = CartesianIndices((1:order,))
        elseif BC.side == Left
            I = CartesianIndices((length(G2)-order+1:length(G2),))
        end

        new{TT,1,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint,1,I)
    end
    function InterfaceBoundaryData{TT}(G1::Grid2D,G2::Grid2D,BC,Joint,order::Int) where {TT}

        if BC.side ∈ [Left,Right]
            n = G2.ny
            BufferIn    = zeros(TT,(order,n))
            BufferOut   = zeros(TT,(order,n))

            if BC.side == Left
                I = CartesianIndices((1:order,1:n))
            else
                I = CartesianIndices((G2.nx-order+1:G2.nx,1:n))
            end

        elseif BC.side ∈ [Up,Down]
            n = G2.nx
            BufferIn    = zeros(TT,(n,order))
            BufferOut   = zeros(TT,(n,order))

            if BC.side == Up
                I = CartesianIndices((1:n,1:order))
            else
                I = CartesianIndices((1:n,G2.ny-order+1:G2.ny))
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
        if !(typeof(P.BoundaryConditions.BoundaryLeft) <: SimultanousApproximationTerm{:Periodic})
            BCL = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryLeft,1,P.order)
            BCR = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryRight,1,P.order)
        else
            BCL = _newBoundaryCondition(G,G,P.BoundaryConditions.BoundaryLeft,1,P.order)
            BCR = _newBoundaryCondition(G,G,P.BoundaryConditions.BoundaryRight,1,P.order)
        end
        if DIM == 1
            # NB = 2
            BCU = nothing
            BCD = nothing
            # BC = (BCL,BCR)
        elseif DIM ==2
            # NB = 4
            if !(typeof(P.BoundaryConditions.BoundaryUp) <: SimultanousApproximationTerm{:Periodic})
                BCU = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryUp,1,P.order)
                BCD = _newBoundaryCondition(G,P.BoundaryConditions.BoundaryDown,1,P.order)
            else
                BCU = _newBoundaryCondition(G,G,P.BoundaryConditions.BoundaryUp,1,P.order)
                BCD = _newBoundaryCondition(G,G,P.BoundaryConditions.BoundaryDown,1,P.order)
            end
            # BC = (BCL,BCR,BCU,BCD)
        end
        new{1,typeof(BCL),typeof(BCR),typeof(BCU),typeof(BCD)}(BCL,BCR,BCU,BCD)
    end
    """
        newBoundaryConditions(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Int64)
    """
    function newBoundaryConditions(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Int64) where {TT}

        J = G.Joint[I]
        if typeof(J) <: Joint
            if J.side == Left
                BCLt = SAT_Interface(G.Grids[J.index].Δx,G.Grids[I].Δx,J.side,1,GetOrder(P.order))
                BCL = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCLt,J.index,P.order)

                if P.BoundaryConditions.BoundaryRight.type != Periodic
                    BCR = _newBoundaryCondition(G.Grids[I],P.BoundaryConditions.BoundaryRight,J.index,P.order)
                else
                    BCR = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],P.BoundaryConditions.BoundaryRight,J.index,P.order)
                end

            elseif J.side == Right
                BCRt = SAT_Interface(G.Grids[I].Δx,G.Grids[J.index].Δx,J.side,1,GetOrder(P.order))
                BCR = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCRt,J.index,P.order)

                if P.BoundaryConditions.BoundaryLeft.type != Periodic
                    BCL = _newBoundaryCondition(G.Grids[I],P.BoundaryConditions.BoundaryLeft,J.index,P.order)
                else
                    BCL = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],P.BoundaryConditions.BoundaryLeft,J.index,P.order)
                end
            end
        elseif length(J) == 2 #TODO

            JL = G.Joint[J[1].index]
            JR = G.Joint[J[2].index][1]
            length(JL) == 1 ? JL = JL[1] : JL = JL[2] # Correct for second node where `JL = ((1,Right),)`

            # println(J)
            # println(JL)
            # println(JR)

            BCLt = SAT_Interface(G.Grids[J[1].index].Δx,G.Grids[JL.index].Δx,Left,1,GetOrder(P.order))
            BCRt = SAT_Interface(G.Grids[J[2].index].Δx,G.Grids[JR.index].Δx,Right,1,GetOrder(P.order))

            BCL = _newBoundaryCondition(G.Grids[I],G.Grids[J[1].index],BCLt,J[1].index,P.order)
            BCR = _newBoundaryCondition(G.Grids[I],G.Grids[J[2].index],BCRt,J[2].index,P.order)
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

        # Loop over joints
        if typeof(J) <: Joint
            # If there is a single joint
            for i in 1:length(faces)
                if J.side == faces[i]
                    BCt = SAT_Interface(G.Grids[J.index].Δx,G.Grids[I].Δx,J.side,GetAxis(J.side),GetOrder(P.order))
                    BC = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCt,J.index,P.order)
                else

                    BCt = getfield(P.BoundaryConditions,Pfaces[i])
                    if BCt.type != Periodic
                        BC = _newBoundaryCondition(G.Grids[I],BCt,J.index,P.order)
                    else
                        BC = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCt,J.index,P.order)
                    end

                    # BC = _newBoundaryCondition(G.Grids[I],BCt,J.index,P.order)
                end
                push!(BCS,BC)
            end
        else
            for j in 1:length(J)
                # If there are multiple joints
                JC = J[j]
                for i in 1:length(faces)
                    if JC.side == faces[i]
                        BCt = SAT_Interface(G.Grids[J.index].Δx,G.Grids[I].Δx,J.side,GetAxis(J.side),GetOrder(P.order))
                        BC = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCt,J.index,P.order)
                        push!(BCS,BC)
                    else
                        BCt = getfield(P.BoundaryConditions,Pfaces[i])
                        if BCt.type != Periodic
                            BC = _newBoundaryCondition(G.Grids[I],BCt,J.index,P.order)
                        else
                            BC = _newBoundaryCondition(G.Grids[I],G.Grids[J.index],BCt,J.index,P.order)
                        end
                        # BC = _newBoundaryCondition(G.Grids[I],BCt,J.index,P.order)
                        push!(BCS,BC)
                    end
                end
            end
        end
        new{2,typeof(BCS[1]),typeof(BCS[2]),typeof(BCS[3]),typeof(BCS[4])}(BCS[1],BCS[2],BCS[3],BCS[4])
    end
end


function getjoint(BC::BoundaryData) end
getjoint(BC::InterfaceBoundaryData) = BC.Joint



function Base.iterate(B::newBoundaryConditions{DIM,BCLT,BCRT,BCUT,BCDT},state=0) where {DIM,BCLT,BCRT,BCUT,BCDT}
    state >= 2DIM && return
    return Base.getfield(B,state+1), state+1
end
Base.getindex(B::newBoundaryConditions,N::NodeType{:Left})  = B.BC_Left
Base.getindex(B::newBoundaryConditions,N::NodeType{:Right}) = B.BC_Right
Base.getindex(B::newBoundaryConditions,N::NodeType{:Up})    = B.BC_Up
Base.getindex(B::newBoundaryConditions,N::NodeType{:Down})  = B.BC_Down

_newBoundaryCondition(G1::GridType{TT},G2::GridType{TT},BC::SimultanousApproximationTerm{:Periodic},J,order) where TT  = InterfaceBoundaryData{TT}(G1,G2,BC,J,order)
_newBoundaryCondition(G1::GridType{TT},G2::GridType{TT},BC::SimultanousApproximationTerm{:Interface},J,order) where TT = InterfaceBoundaryData{TT}(G1,G2,BC,J,order)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Dirichlet},J,order) where TT = BoundaryData(G,BC,J,order)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Neumann},J,order) where TT   = BoundaryData(G,BC,J,order)





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
                Kxy = -P.Kxy * B[1]*B[2]/NB
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

            x,y = G[i]
            if (x == 0.5 && y == 0.5) | (x==-0.5 && y==-0.5) #NIMROD TEST
                K[1][i] = P.Kx * 0.5
                K[2][i] = P.Ky * 0.5
                if existb
                    K[3][i] = P.Kx*0.5
                end
            elseif (x==0.5 && y==-0.5) | (x==-0.5 && y==0.5)
                K[1][i] = P.Kx * 0.5
                K[2][i] = P.Ky * 0.5
                if existb
                    K[3][i] = -P.Kx*0.5
                end
            elseif x==0 && y==0
                K[1][i] = P.Kx*0.5
                K[2][i] = P.Ky*0.5
                if existb
                    K[3][i] = P.Kx*0.5
                end
            end

            # if (x==0) && (y==0) #SINGLE ISLAND TEST
            #     K[1][i] = P.Kx
            #     K[2][i] = P.Ky
            #     if !(PT<:Nothing)
            #         K[3][i] = 0.0
            #     end
            # end

            # if G[i] == (0.5,0.5)
            #     println("Kx = ",K[1][i])
            #     println("Ky = ",K[2][i])
            #     println("Kxy = ",K[3][i])
            # end
        end
        # for I in eachrow(K[1])
        #     println(I)
        # end
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
    _newLocalDataBlockBlocks build all necessary blocks for a local data block
"""
function _newLocalDataBlockBlocks(G::LocalGridType{TT,DIM,MET},Para::PT) where {TT,DIM,MET,PT}
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
    newLocalDataBlock{TT,DIM,AT,KT,DCT,GT,BT,DT}
Contains the data for a local block
"""
mutable struct newLocalDataBlock{TT<:Real,
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
    newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType) where {TT}
Initialise a data block for a 1D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function newLocalDataBlock(P::newPDEProblem{TT,1},G::LocalGridType,SC::StepConfig) where {TT}

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

    BStor = newBoundaryConditions(P,G)
    BS = (BStor.BC_Left,BStor.BC_Right)

    IP = innerH(G.Δx,G.n,P.order)
    D = DiffusionOperator(G.n,G.Δx,P.order,false,:Constant)
    PMap = nothing
    source = P.source

    return newLocalDataBlock{TT,1,:Constant,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end
"""

    newLocalDataBlock(P::newPDEProblem{TT,2},G::LocalGridType) where {TT}
Initialise a data block for a 2D problem with only 1 grid.

    *THIS METHOD IS PRIMARILY FOR TESTING*
"""
function newLocalDataBlock(P::newPDEProblem{TT,2},G::LocalGridType,SC::StepConfig) where {TT}

    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(G)
    
    K = _BuildDiffusionMatrix(G,P,P.Parallel)

    PK = (P.Kx,P.Ky)

    BStor = newBoundaryConditions(P,G)
    BS = (BStor.BC_Left,BStor.BC_Right,BStor.BC_Up,BStor.BC_Down)
    # BS = (BStor.BC_Left,BStor.BC_Right,nothing,nothing); @warn "Using periodic stencil"
    IP = innerH(G.Δx,G.Δy,G.nx,G.ny,P.order)

    # @show GetMetricType(G)

    # if length(K) == 3
    if GetMetricType(G) == CurvilinearMetric
        difftype = :Variable
    else
        difftype = :Constant
    end

    typeof(P.BoundaryConditions.BoundaryLeft).parameters[2] == :Cartesian ? sattype = :Constant : sattype = :Variable
    sattype = :Variable
    # D = DerivativeOperator{TT,2,typeof(P.order),:Constant}(P.order,G.nx,G.ny,G.Δx,G.Δy,false,false)
    Dx = DiffusionOperator(G.nx,G.Δx,P.order,false,difftype)
    Dy = DiffusionOperator(G.ny,G.Δy,P.order,false,difftype)
    D = DiffusionOperatorND(Dx,Dy)
    PMap = P.Parallel
    source = P.source

    return newLocalDataBlock{TT,2,sattype,typeof(u),typeof(K),typeof(PK),typeof(G),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,G,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
end


"""
    newLocalDataBlock(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
Initialise a data block for a 1D multiblock problem
"""
function newLocalDataBlock(P::newPDEProblem{TT,1},G::GridMultiBlock,I::Integer) where {TT}
    u, uₙ₊₁, cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(G.Grids[I])

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
    IP = innerH(G.Grids[I].Δx,G.Grids[I].n,GetOrder(P.order))
    D = DerivativeOperator{TT,1,typeof(P.order),:Constant}(P.order,G.Grids[I].n,0,G.Grids[I].Δx,TT(0),false,false)
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,1,:Cartesian,typeof(u),typeof(K),typeof(P.K),typeof(G.Grids[I]),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K, P.K, G.Grids[I], BS, D,source,PMap, IP, cache,rₖ,dₖ,b,SC)
end
"""
    newLocalDataBlock(P::newPDEProblem{TT,2},G::GridMultiBlock,I::Integer)
Initialise a data block for a 2D multiblock problem
"""
function newLocalDataBlock(P::newPDEProblem{TT,2},G::GridMultiBlock{TT,2,MET},I::Integer) where {TT,MET}
    LG = G.Grids[I]

    u, uₙ₊₁, K , cache, rₖ, dₖ, b = _BuildGenericLocalDataBlocks(LG)
    
    _setKoefficient!(K,P,LG)
    PK = (P.Kx,P.Ky)

    if typeof(P.Parallel) <: Vector
        PMap = P.Parallel[I]
    else
        PMap = P.Parallel
    end

    BStor = newBoundaryConditions(P,G,I)
    BS = (BStor.BC_Left,BStor.BC_Right,BStor.BC_Up,BStor.BC_Down)
    IP = innerH(LG.Δx,LG.Δy,LG.nx,LG.ny,GetOrder(P.order))
    D = DerivativeOperator{TT,2,typeof(P.order),:Constant}(P.order,LG.nx,LG.ny,LG.Δx,LG.Δy,false,false)
    source = SourceTerm{Nothing}(nothing)
    SC = StepConfig{TT}()

    return newLocalDataBlock{TT,2,:Cartesian,typeof(u),typeof(K),typeof(PK),typeof(LG),typeof(BS),typeof(D),typeof(source),typeof(PMap)}(u,uₙ₊₁,K,PK,LG,BS,D,source,PMap,IP,cache,rₖ,dₖ,b,SC)
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
        TDBLOCK <: NTuple{NB,LocalDataBlockType}
            } <: DataBlockType{TT,DIM}
    
    Block   :: TDBLOCK
    SC      :: StepConfig{TT}
    nblock  :: Int64
    parallel:: Bool

    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::LocalGridType{TT},Δt::TT,t::TT;θ=TT(1)) where {TT,DIM}
        SC = StepConfig{TT}(t,Δt,θ)
        DTA = (newLocalDataBlock(P,G,SC),)
        new{TT,DIM,1,typeof(DTA)}(DTA,SC,length(DTA),false)
    end
    function DataMultiBlock(P::newPDEProblem{TT,DIM},G::GridMultiBlock{TT,DIM},Δt::TT,t::TT) where {TT,DIM}
        DTA = map((x)->newLocalDataBlock(P,G,x),eachgrid(G))
        DTA = tuple(DTA...)
        SC = StepConfig{TT}(t,Δt)
        new{TT,DIM,length(DTA),typeof(DTA)}(DTA,SC,length(DTA),false)
    end
    function DataMultiBlock(D::newLocalDataBlock{TT,DIM}) where {TT,DIM} #== TESTING METHOD ==#
        DTA = (D,)
        new{TT,DIM,1,typeof(DTA)}(DTA,D.SC,length(DTA),false)
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
function collectΔu(D::newLocalDataBlock) end
