
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

    BoundaryOperator:: BCT
    RHS             :: F1
    BufferRHS       :: AT
    X               :: Union{Vector{TT},Vector{NTuple{2,TT}}}   # Grid points along boundary
    n               :: Int64        # Length of boundary
    DIM             :: Int64
end
function BoundaryData(G::Grid1D{TT},BC,order::Int64) where {TT}

    BufferRHS = zeros(TT,1)

    BoundaryData{TT,1,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,[TT(0)],1,1)
end
function BoundaryData(G::Grid2D{TT,COORD},BC,order::Int64) where {TT,COORD}

    # Store the boundary points and generate buffers for storing RHS
    if BC.side ∈ [Left,Right]
        n = G.ny
        BufferRHS = zeros(TT,(1,n))
    elseif BC.side ∈ [Up,Down]
        n = G.nx
        BufferRHS = zeros(TT,(n,1))
    end

    X = GetBoundaryCoordinates(G,BC.side)
    BoundaryData{TT,2,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,X,n,2)
end
# function BoundaryData(G::Grid2D{TT,CurvilinearMetric},BC,order::Int64) where {TT}

#     # Store the boundary points and generate buffers for storing RHS
#     if BC.side ∈ [Left,Right]
#         n = G.ny
#         BufferRHS = zeros(TT,(1,n))
#     elseif BC.side ∈ [Up,Down]
#         n = G.nx
#         BufferRHS = zeros(TT,(n,1))
#     end

#     X = GetBoundaryCoordinates(G,BC.side)
#     BoundaryData{TT,2,typeof(BC.RHS),typeof(BC),typeof(BufferRHS)}(BC,BC.RHS,BufferRHS,X,n,2)
# end


"""
    InterfaceBoundaryData
Container for Periodic and Interface boundary conditions.
`OutgoingJoint` is the joint FROM this block TO the other block.
`IncomingJoint` is the joint TO this block FROM the other block.
"""
struct InterfaceBoundaryData{
        TT<:Real,
        DIM,
        BCT,
        AT} <: BoundaryStorage{TT,0,AT}

    BoundaryOperator:: BCT
    BufferOut       :: AT
    BufferIn        :: AT
    OutgoingJoint   :: Joint
    IncomingJoint   :: Joint
end
function InterfaceBoundaryData{TT}(G1::Grid1D,G2::Grid1D,BC,Joint1::Joint,Joint2::Joint) where {TT}

    BufferIn    = zeros(TT,2)
    BufferOut   = zeros(TT,2)

    InterfaceBoundaryData{TT,1,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint1,Joint2)
end
function InterfaceBoundaryData{TT}(G1::Grid2D,G2::Grid2D,BC,Joint1::Joint,Joint2::Joint) where {TT}
    if BC.side ∈ [Left,Right]
        n = G1.ny
        BufferIn    = zeros(TT,(2,n))
        BufferOut   = zeros(TT,(2,n))

    elseif BC.side ∈ [Up,Down]
        n = G1.nx
        BufferIn    = zeros(TT,(n,2))
        BufferOut   = zeros(TT,(n,2))
    end

    InterfaceBoundaryData{TT,2,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint1,Joint2)
end
function InterfaceBoundaryData{TT}(G1::Grid2D,BC) where {TT}
    if BC.side ∈ [Left,Right]
        n = G1.ny
        BufferIn    = zeros(TT,(BC.order,n))
        BufferOut   = zeros(TT,(BC.order,n))

    elseif BC.side ∈ [Up,Down]
        n = G1.nx
        BufferIn    = zeros(TT,(n,BC.order))
        BufferOut   = zeros(TT,(n,BC.order))
    end

    InterfaceBoundaryData{TT,2,typeof(BC),typeof(BufferOut)}(BC,BufferOut,BufferIn,Joint(0,Left),Joint(0,Left))
end


"""
    GenerateBoundaries
"""
function GenerateBoundaries end
function GenerateBoundaries(P::Problem1D,G::LocalGridType{TT,1}) where TT
    tmpDict = Dict()
    for BC in P.BoundaryConditions
        if typeof(BC) <: SAT_Periodic
            tmpDict[BC.side] = _newBoundaryCondition(G,G,BC,BC.side,_flipside(BC.side),P.order)
        else
            tmpDict[BC.side] = BoundaryData(G,BC,P.order)
        end
    end

    return tmpDict
end
function GenerateBoundaries(P::Problem1D,G::GridMultiBlock{TT,1},I::Int64) where TT
    jts = G.Joint[I]

    tmpDict = Dict()

    for Joint in jts
        BC = SAT_Interface(G.Grids[Joint.index].Δx,G.Grids[I].Δx,Joint.side,1,P.order)
        tmpDict[Joint.side] = _newBoundaryCondition(G.Grids[I],G.Grids[Joint.index],BC,Joint.index,P.order)
    end

    if haskey(P.BoundaryConditions,I)
        for BC in P.BoundaryConditions[I]
            if typeof(BC) <: SAT_Periodic
                # In 1D the boundary must be on the other side
                I - 1 == 0 ? j = G.ngrids : j = 1 # Tell the periodic BC where to find the source data
                tmpDict[BC.side] = _newBoundaryCondition(G.Grids[I],BC)
            else
                # If non-periodic we can just assign the BC
                tmpDict[BC.side] = BoundaryData(G.Grids[I],BC,P.order)
            end
        end
    end

    return tmpDict
end
function GenerateBoundaries(P::Problem2D,G::LocalGridType{TT,2},K) where TT
    tmpDict = Dict()
    for BC in P.BoundaryConditions
        if typeof(BC) <: SAT_Periodic
            tmpDict[BC.side] = _newBoundaryCondition(G,BC)
        else
            if (GetMetricType(G) == CurvilinearMetric) & (typeof(BC) <: SAT_Dirichlet)
                if typeof(BC.side).parameters[2] == 1
                    typeof(BC.side).parameters[1] == :Left ? τ = minimum(K[1][1,:]) : τ = minimum(K[1][end,:])
                else
                    typeof(BC.side).parameters[1] == :Left ? τ = minimum(K[2][:,1]) : τ = minimum(K[2][:,end])
                end
                BC.τ₀[1] = -(10 + 1/τ)
            end
            tmpDict[BC.side] = BoundaryData(G,BC,P.order)
        end
    end

    return tmpDict
end
"""
2D multiblock version of GenerateBoundaries.
"""
function GenerateBoundaries(P::Problem2D,G::GridMultiBlock{TT,2,COORD},I::Int64,K) where {TT,COORD}
    Joints = G.Joint[I]

    tmpDict = Dict{NodeType,BoundaryStorage}()

    COORD == CurvilinearMetric ? sattype = :Curvilinear : sattype = :Cartesian

    for Joint in Joints
        normal = TT(1)
        NeighbouringJoint = nothing
        for tmpJoint in G.Joint[Joint.index] # Scan all joints to work out how the other block sees this one
            if (tmpJoint.index == I)
                NeighbouringJoint = tmpJoint
                if (_flipside(tmpJoint.side) != Joint.side)# & (Joint.index < I) # flip the normal vector if required
                    if typeof(Joint.side).parameters[1] == typeof(tmpJoint.side).parameters[1]
                        # Flipping normal
                        normal = -TT(1)
                    end
                end
            end
        end

        G1 = G.Grids[I]
        G2 = G.Grids[Joint.index] #Neighbouring grid

        MySide = Joint.side
        NeighbourSide = NeighbouringJoint.side

        # Select the correct Δx in the neighbouring domain
        if NeighbourSide ∈ [Left, Right]
            Δx₂ = G2.Δx
        elseif NeighbourSide ∈ [Up, Down]
            Δx₂ = G2.Δy
        end

        if (Joint.side == Left) || (Joint.side == Right)
            Δx₁ = G1.Δx# * G.Grids[I].qx[1]
            Δy₁ = G1.Δy

            buffer = zeros(TT,(1,G1.ny))

            @. buffer[1,:] = P.Kx * max(G1.J[MySide] * (G1.qx[MySide]^2 + G1.qy[MySide]^2),
                G2.J[NeighbourSide] * (G2.qx[NeighbourSide]^2 + G2.qy[NeighbourSide]^2))

        elseif (Joint.side == Down) | (Joint.side == Up)
            Δx₁ = G1.Δy
            Δy₁ = G1.Δx

            buffer = zeros(TT,(G1.nx,1))

            @. buffer[:,1] = P.Ky * max(G1.J[MySide] * (G1.qx[MySide]^2 + G1.qy[MySide]^2),
                G2.J[NeighbourSide]*(G2.qx[NeighbourSide]^2 + G2.qy[NeighbourSide]^2))
        end
        τ₀ = maximum(buffer)

        BC = SAT_Interface(Δx₁,Δx₂,τ₀,Joint.side,P.order,Δy=Δy₁,coordinates=sattype,normal=normal)
        tmpDict[Joint.side] = _newBoundaryCondition(G.Grids[I],G.Grids[Joint.index],BC,Joint,NeighbouringJoint)
    end

    if haskey(P.BoundaryConditions,I)
        for BC in P.BoundaryConditions[I]
            if typeof(BC) <: SAT_Periodic
                # In 2D we need to find which boundary matches the periodic BC
                @warn "Periodic boundary condition for 2D multiblock problem assumes the periodic boundary is contained to one block, if you want to use periodic conditions across blocks you should connect them via a Joint."
                tmpDict[BC.side] = _newBoundaryCondition(G.Grids[I],BC)
            else
                # If non-periodic we can just assign the BC
                if (GetMetricType(G) == CurvilinearMetric) | (typeof(BC) <: SAT_Dirichlet)
                    if typeof(BC.side).parameters[2] == 1
                        typeof(BC.side).parameters[1] == :Left ? τ = maximum(K[1][1,:]) : τ = maximum(K[1][end,:])
                    else
                        typeof(BC.side).parameters[1] == :Left ? τ = maximum(K[2][:,1]) : τ = maximum(K[2][:,end])
                    end
                    BC.τ₀[1] = -(10 + 1/τ)
                end
                tmpDict[BC.side] = BoundaryData(G.Grids[I],BC,P.order)

                # @show typeof(BoundaryData(G.Grids[I],BC,P.order))

            end
        end
    end

    return tmpDict
    # return (tmpDict[Left],tmpDict[Right],tmpDict[Up],tmpDict[Down])
end

_newBoundaryCondition(G1::GridType{TT},BC::SimultanousApproximationTerm{:Periodic}) where TT  = InterfaceBoundaryData{TT}(G1,BC)
_newBoundaryCondition(G1::GridType{TT},G2::GridType{TT},BC::SimultanousApproximationTerm{:Interface},Joint1,Joint2) where TT = InterfaceBoundaryData{TT}(G1,G2,BC,Joint1,Joint2)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Dirichlet},order) where TT = BoundaryData(G,BC,order)
_newBoundaryCondition(G::GridType{TT},BC::SimultanousApproximationTerm{:Neumann},order) where TT   = BoundaryData(G,BC,order)
