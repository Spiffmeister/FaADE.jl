

"""
    addSource!
Adds PDE forcing term
"""
function addSource! end
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
function addSource!(S::SourceTerm{F},dest::Symbol,D::DataMultiBlock,G::GridType) where {TT,F<:Function}
    for I in eachblock(D)
        write = getproperty(D[I],dest)
        addSource!(S.F, write, G, D.SC.t, D.SC.Δt)
    end
end
function addSource!(S::SourceTerm{Nothing},dest::Symbol,D::DataMultiBlock,G::GridType) end


"""
BoundaryConditions
Sets the value of the boundary.
"""
function setBoundaryConditions! end
function setBoundaryConditions!(RHS::Function,Bound::AT,t::T,Δt::T) where {AT,T}
    Bound[1] = Δt*RHS(t)
end
function setBoundaryConditions!(RHS::Function,Bound::AT,grid::Vector{T},n::Int,t::T,Δt::T) where {AT,T}
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end
"""
    BoundaryConditions
Sets the value of the boundary.
"""
function setBoundaryCondition! end
function setBoundaryCondition!(B::newBoundaryData{TT,1,Fn},Δt::TT,args...) where {TT,Fn<:Real}
    @. B.BufferRHS = Δt*B.RHS
end
function setBoundaryCondition!(B::newBoundaryData{TT,1,Fn},Δt::TT,t::TT) where {TT,Fn<:Function}
    B.BufferRHS[1] = Δt*B.RHS(t)
end
function setBoundaryCondition!(B::newBoundaryData{TT,2,Fn},Δt::TT,t::TT) where {TT,Fn<:Function}
    for i = 1:B.n
        B.BufferRHS[i] = Δt*B.RHS(B.X[i],t)
    end
end
"""
If no boundary condition set or the boundary condition is periodic or an interface, ignore this step
"""
function setBoundaryCondition!(BC::Nothing,args...) end
function setBoundaryCondition!(BC::newInterfaceBoundaryData,args...) where TT end
"""
Calling boundaries for data blocks
"""
function setBoundaryConditions!(D::newLocalDataBlockType{TT,DIM}) where {TT,DIM}
    # for B in D.boundary
    setBoundaryCondition!(D.boundary.BC_Left,   D.SC.Δt, D.SC.t)
    setBoundaryCondition!(D.boundary.BC_Right,  D.SC.Δt, D.SC.t)
    if DIM == 2
        setBoundaryCondition!(D.boundary.BC_Up,     D.SC.Δt, D.SC.t)
        setBoundaryCondition!(D.boundary.BC_Down,   D.SC.Δt, D.SC.t)
    end
    # end
end
"""
Calling boundaries from multiblocks
"""
function setBoundaryConditions!(D::DataMultiBlock) where TT
    for I in eachblock(D)
        setBoundaryConditions!(D[I])
    end
end


"""
    fillBuffers
"""
function fillBuffers end
function fillBuffers(B::newBoundaryData,args...) end
function fillBuffers(source::Symbol,DB::newLocalDataBlock)
    S = getproperty(DB,source)
    copyto!(DB.boundary.RHS_Left,   S)
    copyto!(DB.boundary.RHS_Right,  S)
end
function fillBuffer!(source,B::newBoundaryData,args...) end
function fillBuffer!(source::AT,B::newInterfaceBoundaryData,K::AT) where AT
    B.BufferIn .= source
end
function fillBuffer(source::Symbol,B::newBoundaryData,DB::DataMultiBlock) end
function fillBuffer(source::Symbol,B::newInterfaceBoundaryData,DB::DataMultiBlock)
    cache = getproperty(DB[B.Joint],source)
    # println(B.I)
    # println("----------")
    copyto!(B.BufferIn,CartesianIndices(B.BufferIn),cache,B.I)
end

function fillBuffers(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(DB)
        # println(I)
        BC = DB[I].boundary.BC_Left
        # S = getproperty(DB[BC.Joint],source)
        fillBuffer(source,BC,DB)

        BC = DB[I].boundary.BC_Right
        # S = getproperty(DB[BC.Joint],source)
        fillBuffer(source,BC,DB)
    end
    # println("buffer ",DB[1].boundary.BC_Right.BufferIn)
end


"""
    applySAT!
"""
function applySAT! end
"""
Apply the SAT
"""
@inline applySAT!(SAT::SimultanousApproximationTerm,dest::AT,K::AT,u::AT,mode::SATMode) where {AT} = SAT(dest,K,u,mode)
@inline applySAT!(SAT::SimultanousApproximationTerm,dest::AT,K::AT,u::AT,buffer::AT,mode::SATMode) where {AT} = SAT(dest,K,u,buffer,mode)
"""
Applying SATs in DataMode ignoring interface terms
"""
@inline function applySAT!(BC::newInterfaceBoundaryData,dest,K,mode::SATMode{:DataMode}) end
@inline function applySAT!(BC::newBoundaryData,dest::AT,K::AT,mode::SATMode{:DataMode}) where {TT,AT}
    applySAT!(BC.Boundary,dest,K,BC.BufferRHS,mode)
    # BC.Boundary(dest,K,BC.BufferRHS,mode)
end
"""
Applying SATs in SolutionMode
"""
@inline function applySAT!(BC::newInterfaceBoundaryData,dest::AT,K::AT,source::AT,mode::SATMode{:SolutionMode}) where {AT}
    applySAT!(BC.Boundary,dest,K,source,BC.BufferIn,mode)
end
@inline function applySAT!(BC::newBoundaryData,dest::AT,K::AT,source::AT,mode::SATMode{:SolutionMode}) where {AT}
    applySAT!(BC.Boundary,dest,K,source,mode)
end
"""
    applySATs
"""
function applySATs end
function applySATs(dest::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT}
        applySAT!(D.boundary.BC_Left,   dest, D.K, mode)
        applySAT!(D.boundary.BC_Right,  dest, D.K, mode)
end
function applySATs(dest::VT,source::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT}
    applySAT!(D.boundary.BC_Left,   dest, D.K, source, mode)
    applySAT!(D.boundary.BC_Right,  dest, D.K, source, mode)
end
function applySATs(dest::AT,D::newLocalDataBlock{TT,2,AT},mode) where {TT,AT}
    applySAT!(D.B.BC_Left,   dest, D.K[1], mode)
    applySAT!(D.B.BC_Right,  dest, D.K[1], mode)
    applySAT!(D.B.BC_Up,     dest, D.K[2], mode)
    applySAT!(D.B.BC_Down,   dest, D.K[2], mode)
end
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,AT},mode) where {TT,AT}
    applySAT!(D.B.BC_Left,   dest, source, D.K[1], mode)
    applySAT!(D.B.BC_Right,  dest, source, D.K[1], mode)
    applySAT!(D.B.BC_Up,     dest, source, D.K[2], mode)
    applySAT!(D.B.BC_Down,   dest, source, D.K[2], mode)
end
"""
Multiblock version
"""
function applySATs(dest::Symbol,D::DataMultiBlock,mode::SATMode)
    for I in eachblock(D)
        write = getproperty(D[I],dest)
        applySATs(write,D[I],mode)
    end
end

"""
    setDiffusionCoefficient!
Sets the diffusion coefficient
"""
function setDiffusionCoefficient! end
function setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid1D)
    for i = 1:grid.n
        K[i] = κ(grid[i])
    end
end
function setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid2D)
    for i = 1:grid.nx
        for j = 1:grid.ny
            K[i,j] = κ(grid[i,j]...)
        end
    end
end
@inline function setDiffusionCoefficient!(κ::TT,K::AbstractArray{TT},grid::GridType) where {TT<:Real} end
@inline setDiffusionCoefficient!(D::newLocalDataBlock,grid::GridType) = setDiffusionCoefficient!(D.κ,D.K,grid)
function setDiffusionCoefficient!(D::DataMultiBlock,grid::GridType)
    for I in eachblock(D)
        setDiffusionCoefficient!(D[I],grid)
    end
end



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
function copyUtoSAT!(D::DataMultiBlock{TT,DIM,DT},order::DerivativeOrder) where {TT,DIM,DT<:newLocalDataBlockType}
    copyUtoSAT!(D.Data.boundary,D.Data.u,GetOrder(order))
end



"""
    muladd!
multiply-add for multiblock problems
"""
function muladd!(dest::AT,source::AT,α::TT,β::TT) where {TT,AT<:AbstractArray{TT}}
    @. dest = α*dest + β*source
end
function muladd!(dest::Symbol,source::Symbol,D::newLocalDataBlock,α::TT,β::TT) where {TT}
    W = getproperty(D,dest)
    R = getproperty(D,source)
    muladd!(W,R,α,β)
end
function muladd!(dest::Symbol,source::Symbol,D::DataMultiBlock{TT};α=TT(1),β=TT(1)) where TT
    for I in eachblock(D)
        muladd!(dest,source,D[I],α,β)
    end
end
"""
    setValue
"""
@inline function setValue(dest::AT,source::AT,α::TT) where {TT<:Real,AT<:AbstractArray{TT}}
    @. dest = α*source
end
function setValue(dest::Symbol,source::Symbol,D::newLocalDataBlock{TT},α::TT) where TT
    W = getproperty(D,dest)
    R = getproperty(D,source)
    setValue(W,R,α)
end
function setValue(dest::Symbol,source::Symbol,D::DataMultiBlock{TT},α=TT(1)) where TT
    for I in eachblock(D)
        setValue(dest,source,D[I],α)
    end
end

"""
    innerprod
"""
innerprod(u::AT,v::AT,IP::innerH{TT}) where {TT,AT<:AbstractArray{TT}} = IP(u,v)
@inline function innerprod(u::Symbol,v::Symbol,D::newLocalDataBlock{TT,DIM}) :: TT where {TT,DIM} 
    U = getproperty(D,u)
    V = getproperty(D,v)
    IP = getproperty(D,:innerprod)
    return IP(U,V)
end
function innerprod(u::Symbol,v::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    local r::TT
    r = TT(0)
    for I in eachblock(DB)
        r += innerprod(u,v,DB[I])
    end
    return r
end


"""
    copyto!
"""
@inline function Base.copyto!(dest,source,D::DataMultiBlock)
    for I in eachblock(D)
        d = getproperty(D[I],dest)
        # d .= getproperty(D[I],source)
        s = getproperty(D[I],source)
        copyto!(d,s)
    end
end