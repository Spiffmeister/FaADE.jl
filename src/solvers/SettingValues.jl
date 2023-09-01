

"""
    addSource!
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
function setBoundaryConditions!(D::newLocalDataBlockType,G::Grid1D,t::TT,Δt::TT) where TT
    setBoundaryConditions!(D.boundary.Left.RHS,    D.boundary.RHS_Left, t,Δt)
    setBoundaryConditions!(D.boundary.Right.RHS,   D.boundary.RHS_Right,t,Δt)
end
function setBoundaryConditions!(D::newLocalDataBlockType,G::Grid2D,t::TT,Δt::TT) where TT
    setBoundaryConditions!(B.BoundaryLeft.RHS, D.Data.RHS_Left,G.nx,G.gridx,t,Δt)
    setBoundaryConditions!(B.BoundaryRight.RHS,D.Data.RHS_Right,G.nx,G.gridx,t,Δt)
    setBoundaryConditions!(B.BoundaryUp.RHS,   D.Data.RHS_Up,  G.ny,G.gridy,t,Δt)
    setBoundaryConditions!(B.BoundaryDown.RHS, D.Data.RHS_Down,G.ny,G.gridy,t,Δt)
end
function setBoundaryConditions!(D::DataMultiBlock,G::GridType)
    for I in eachblock(D)
        setBoundaryConditions!(D[I],   G,  D.SC.t, D.SC.Δt)
    end
end


"""
    fillBuffers
"""
function fillBuffers end
function fillBuffers(source::Symbol,DB::newLocalDataBlock)
    S = getproperty(DB,source)
    copyto!(DB.boundary.RHS_Left, getproperty(DB,source))
    copyto!(DB.boundary.RHS_Right, getproperty(DB,source))
end

function fillBuffers(source,DB::newLocalDataBlock)
    for J in eachjoint(DB)
        SB = getproperty(source,DB)

    end
end

function fillBuffers(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(DB)
        BB = DB[I].boundary
        cache =  getproperty(DB[BB.JointLeft],source)
        copyto!(BB.u_Left,BB.LRStorageIndex,cache,BB.LeftIndexFrom)
        cache =  getproperty(DB[BB.JointRight],source)
        copyto!(BB.u_Right,BB.LRStorageIndex,cache,BB.RightIndexFrom)

        if DIM == 2
            copyto!(BB.u_Up,1,getproperty(DB[BB.JointUp],source),BB.UpIndex)
            copyto!(BB.u_Down,1,getproperty(DB[BB.JointDown],source),BB.DownIndex)
        end

    end
end



"""
    applySATs
"""
function applySAT!(Boundary::Nothing,Block,Prob,mode) end
function applySAT!(SAT::SimultanousApproximationTerm,cache::AT,u::AT,k::AT,mode::SATMode) where AT
    SAT(cache,k,u,mode)
end
function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,1},mode::SATMode{:DataMode}) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.RHS_Left,    D.K, mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.RHS_Right,   D.K, mode)
end
function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,2},mode::SATMode{:DataMode}) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.RHS_Left,    D.K[1], mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.RHS_Right,   D.K[1], mode)
    applySAT!(D.boundary.Up,     dest, D.boundary.RHS_Left,    D.K[2], mode)
    applySAT!(D.boundary.Down,   dest, D.boundary.RHS_Right,   D.K[2], mode)
end
function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,1},mode::SATMode{:SolutionMode}) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.u_Left,    D.K, mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.u_Right,   D.K, mode)
end
function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,2},mode::SATMode{:SolutionMode}) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.u_Left,   D.K[1], mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.u_Right,  D.K[1], mode)
    applySAT!(D.boundary.Up,     dest, D.boundary.u_Up,     D.K[2], mode)
    applySAT!(D.boundary.Down,   dest, D.boundary.u_Down,   D.K[2], mode)
end
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
function setDiffusionCoefficient!(κ::TT,K::AbstractArray{TT},grid::GridType) where {TT<:Real}
    K .= κ
end
function setDiffusionCoefficient!(D::newLocalDataBlock,grid::GridType)
    setDiffusionCoefficient!(D.κ,D.K,grid)
end
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
@inline function muladd!(dest::Array{TT},source::Array{TT},α::TT,β::TT) where {TT}
    @. dest = α*dest + β*source
end
@inline function muladd!(dest::Symbol,source::Symbol,D::DataMultiBlock{TT};α=TT(1),β=TT(1)) where TT

    for I in eachblock(D)
        W = getproperty(D[I],dest)
        A = getproperty(D[I],source)
        muladd!(W,A,α,β)
        # @. W = α*W + β*A
    end
end
"""
    setValue
"""
@inline function setValue(dest::AT,source::AT,α::TT) where {AT,TT}
    @. dest = α*source
end
@inline function setValue(dest::Symbol,source::Symbol,D::newDataBlockType{TT},α=TT(1)) where TT
    
    for I in eachblock(D)
        W = getproperty(D[I],dest)
        R = getproperty(D[I],source)

        setValue(W,R,α)
    end
end

"""
    innerprod
"""
innerprod(u::Array{TT},v::Array{TT},IP::innerH) where TT = IP(u,v)
function innerprod(u::Symbol,v::Symbol,D::DataMultiBlock{TT,DIM}) where {TT,DIM}
    local r::TT
    r = TT(0)
    for I in eachblock(D)
        # r += D[I].innerprod(getproperty(D[I],u),getproperty(D[I],v))
        U = getproperty(D[I],u)
        V = getproperty(D[I],v)
        r += innerprod(U,V,D[I].innerprod)
    end
    return r
end


"""
    copyto!
"""
@inline function Base.copyto!(dest,source,D::DataMultiBlock)
    for I in eachblock(D)
        d = getproperty(D[I],dest)
        d .= getproperty(D[I],source)
    end
end