

"""
    addSource!
Adds PDE forcing term
"""
function addSource! end
function addSource!(S::Function,u::AbstractArray{TT},grid::Grid2D{TT},t::TT,Δt::TT) where TT
    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] += Δt*S(grid.gridx[i],grid.gridy[j],t)
        end
    end
    u
end
function addSource!(S::SourceTerm{Function},u::AbstractArray{TT},grid::Grid1D{TT},t::TT,Δt::TT) where TT
    u .+= Δt*S.source.(grid.grid,t)
end
function addSource!(S::SourceTerm{Function},u::AbstractArray{TT},grid::Grid2D{TT},t::TT,Δt::TT) where TT
    for j in 1:grid.ny
        for i in 1:grid.nx
            u[i,j] += Δt*S.source(grid.gridx[i],grid.gridy[j],t)
        end
    end
end
function addSource!(dest::Symbol,D::DataMultiBlock)
    for I in eachblock(D)
        write = getproperty(D[I],dest)
        addSource!(D[I].source, write, D[I].grid, D[I].SC.t, D[I].SC.Δt)
    end
end
function addSource!(S::SourceTerm{Nothing},tmp...) end


"""
    setBoundaryConditions!
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
function setBoundaryCondition!(BC::newInterfaceBoundaryData,args...) end
"""
Calling boundaries for data blocks
"""
function setBoundaryConditions!(D::newLocalDataBlockType{TT,DIM}) where {TT,DIM}
    # for I in eachboundary(D)
        setBoundaryCondition!(D.boundary[1], D.SC.Δt, D.SC.t)
        setBoundaryCondition!(D.boundary[2], D.SC.Δt, D.SC.t)
        setBoundaryCondition!(D.boundary[3], D.SC.Δt, D.SC.t)
        setBoundaryCondition!(D.boundary[4], D.SC.Δt, D.SC.t)
    # end
end
"""
Calling boundaries from multiblocks
"""
function setBoundaryConditions!(D::DataMultiBlock)
    for I in eachblock(D)
        setBoundaryConditions!(D[I])
    end
end



"""
    fillBuffer!
"""
function fillBuffer! end
function fillBuffer!(source,B::newBoundaryData,args...) end
function fillBuffer!(source::AT,B::newInterfaceBoundaryData,K::AT) where AT
    B.BufferIn .= source
end
"""
    fillBuffer
"""
function fillBuffer end
function fillBuffer(source::Symbol,B::newBoundaryData,DB::DataMultiBlock) end
function fillBuffer(source::Symbol,B::newInterfaceBoundaryData,DB::DataMultiBlock)
    cache = getproperty(DB[B.Joint],source)
    copyto!(B.BufferIn,CartesianIndices(B.BufferIn),cache,B.I)
end
function fillBuffer(source::Symbol,DB::DataMultiBlock,I::Int64,side::NodeType)
    BC = DB[I].boundary[side]
    fillBuffer(source,BC,DB)
end
"""
    fillBuffers
"""
function fillBuffers end
function fillBuffers(B::newBoundaryData,args...) end
function fillBuffers(source::Symbol,DB::newLocalDataBlock{TT,DIM}) where {TT,DIM}
    S = getproperty(DB,source)
    copyto!(DB.boundary.RHS_Left,   S)
    copyto!(DB.boundary.RHS_Right,  S)
    if DIM == 2
        copyto!(DB.boundary.RHS_Up,     S)
        copyto!(DB.boundary.RHS_Down,   S)
    end
end
function fillBuffers(source::Symbol,BC::NTuple{NB,BoundaryStorage},DB::DataMultiBlock{TT}) where {TT,NB}
    fillBuffer(source,BC[1],DB)
    fillBuffer(source,BC[2],DB)
    if NB == 4
        fillBuffer(source,BC[3], DB)
        fillBuffer(source,BC[4], DB)
    end
end
function fillBuffer(source::Symbol,DB::DataMultiBlock{TT,DIM},I::Int64) where {TT,DIM}
    D = DB[I]
    fillBuffer(source,D.boundary[1],DB)
    fillBuffer(source,D.boundary[2],DB)
    if DIM == 2
        fillBuffer(source,D.boundary[3],DB)
        fillBuffer(source,D.boundary[4],DB)
    end
end
function fillBuffers(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(DB)
        fillBuffer(source,DB,I)
    end
end

# @generated function fillBuffers(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
#     quote
#         for I in eachblock(DB)
#             fillBuffers(source,DB[I].boundary,DB)
#         end
#     end
# end




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
@inline function applySAT!(BC::newBoundaryData,dest::AT,K::AT,mode::SATMode{:DataMode}) where {AT}
    # if (BC.Boundary.side == Left)
    if typeof(BC.Boundary) <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_implicit_data!(dest,BC.BufferRHS,K,BC.Boundary,mode)
    elseif typeof(BC.Boundary) <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_implicit_data!(dest,BC.BufferRHS,K,BC.Boundary,mode)
    else
        applySAT!(BC.Boundary,dest,K,mode)
    end
    # else
        # applySAT!(BC.Boundary,dest,K,mode)
        # applySAT!(BC.Boundary,dest,K,BC.BufferRHS,mode)
    # end
    # BC.Boundary(dest,K,BC.BufferRHS,mode)
end
"""
Applying SATs in SolutionMode
"""
@inline function applySAT!(BC::newInterfaceBoundaryData,dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {AT}
    # applySAT!(BC.Boundary,dest,K,source,BC.BufferIn,mode)
    if typeof(BC.Boundary) <: SimultanousApproximationTerm{:Interface}
        SAT_Interface!(dest,source,K,BC.BufferIn,BC.Boundary,mode)
    elseif typeof(BC.Boundary) <: SimultanousApproximationTerm{:Periodic}
    end
end
function applySAT!(BC::newBoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {TT,AT,DIM,FT,BCT}
    # if (BC.Boundary.side ∈ [Left,Up]) & (typeof(BC.Boundary) <: SimultanousApproximationTerm{:Dirichlet})
    # println(BCT)
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_implicit!(dest,source,K,BC.Boundary,mode)
    else
        error("Not implemented")
    end
end
"""
    applySATs
"""
function applySATs end
function applySATs(dest::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT}
    for I in eachboundary(D)
        applySAT!(D.boundary[I],  dest, D.K, mode)
    end
        # applySAT!(D.boundary[1],  dest, D.K, mode)
        # applySAT!(D.boundary[2],  dest, D.K, mode)
end
function applySATs(dest::VT,source::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT}
    # for I in eachboundary(D)
    #     applySAT!(D.boundary[I],   dest, D.K, source, mode)
    # end
    
    applySAT!(D.boundary[1],  dest, source, D.K, mode)
    applySAT!(D.boundary[2],  dest, source, D.K, mode)
    
    # applySAT!(D.boundary.BC_Left,   dest, D.K, source, mode)
    # applySAT!(D.boundary.BC_Right,  dest, D.K, source, mode)
end
function applySATs(dest::AT,D::newLocalDataBlock{TT,2,AT},mode) where {TT,AT}
    applySAT!(D.boundary[1],    dest, D.K[1], mode)
    applySAT!(D.boundary[2],    dest, D.K[1], mode)
    applySAT!(D.boundary[3],    dest, D.K[2], mode)
    # println(D.boundary[4])
    applySAT!(D.boundary[4],    dest, D.K[2], mode)
    # println(dest)
    # println()

    # dest[1,:] .= 0.0
    # dest[end,:] .= 0.0
    # dest[:,1] .= 0.0
    # dest[:,end] .= 0.0
end
"""
    applySATs
applySATs for 2D local block
"""
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,AT},mode) where {TT,AT}
    applySAT!(D.boundary[1],    dest, source, D.K[1], mode)
    applySAT!(D.boundary[2],    dest, source, D.K[1], mode)
    applySAT!(D.boundary[3],    dest, source, D.K[2], mode)
    applySAT!(D.boundary[4],    dest, source, D.K[2], mode)
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
Sets the diffusion coefficient for real or functional diffusion coefficients
"""
function setDiffusionCoefficient! end
"""
    setDiffusionCoefficient!
Sets the diffusion coefficient for scalar diffusion coefficient
"""
@inline function setDiffusionCoefficient!(k::TT,K::AbstractArray{TT},grid::GridType) where {TT<:Real}
    @. K = k
end
"""
    setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid1D)
1D functional diffusion coefficient
"""
@inline function setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid1D)
    for i = 1:grid.n
        K[i] = κ(grid[i])
    end
end
"""
    setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid2D)
2D functional diffusion coefficient
"""
@inline function setDiffusionCoefficient!(κ::Function,K::AbstractArray,grid::Grid2D)
    for i = 1:grid.nx
        for j = 1:grid.ny
            K[i,j] = κ(grid[i,j]...)
        end
    end
end
"""
    setDiffusionCoefficient!(D::newLocalDataBlock{TT,1})
Set diffusion coefficient for 1D grids
"""
setDiffusionCoefficient!(D::newLocalDataBlock{TT,1}) where TT = setDiffusionCoefficient!(D.κ,D.K,D.grid)
"""
    setDiffusionCoefficient!(D::newLocalDataBlock{TT,2})
Set diffusion coefficient for 2D grids
"""
function setDiffusionCoefficient!(D::newLocalDataBlock{TT,2}) where TT
    setDiffusionCoefficient!(D.κ[1],D.K[1],D.grid)
    setDiffusionCoefficient!(D.κ[2],D.K[2],D.grid)
    # for i in eachdim(DIM)
    #     setDiffusionCoefficient!(D.κ[i],D.K[i],D.grid)
    # end
end
"""
    setDiffusionCoefficient!(D::DataMultiBlock)
Set diffusion coefficient for each block
"""
function setDiffusionCoefficient!(D::DataMultiBlock)
    for I in eachblock(D)
        setDiffusionCoefficient!(D[I])
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
@inline function muladd!(dest::AT,source::AT,α::TT,β::TT) where {TT,AT<:AbstractArray{TT}}
    @. dest = α*dest + β*source
end
function muladd!(dest::Symbol,source::Symbol,D::newLocalDataBlock,α,β)
    W = getarray(D,dest)
    R = getarray(D,source)
    @. W = α*W + β*R
    W
end
function muladd!(dest::Symbol,source::Symbol,D::DataMultiBlock{TT};α=TT(1),β=TT(1)) where TT
    for I in eachblock(D)
        muladd!(dest,source,D[I],α,β)
    end
end

"""
    setValue
"""
function setValue(dest::AT,source::AT,α::TT) where {TT<:Real,AT<:AbstractArray{TT}}
    @. dest = α*source
end
function setValue(dest::Symbol,source::Symbol,D::newLocalDataBlock,α)
    W = getarray(D,dest)
    R = getarray(D,source)
    @. W = α*R
    # setValue(W,R,α)
    W
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
function innerprod(u::Symbol,v::Symbol,D::newLocalDataBlock{TT}) :: TT where {TT}
    V = getarray(D,v)
    U = getarray(D,u)
    IP = getproperty(D,:innerprod)
    ret = innerprod(U,V,IP)
    return ret
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
@inline function Base.copyto!(dest::Symbol,source::Symbol,D::newLocalDataBlock)
    d = getarray(D,dest)
    s = getarray(D,source)
    copyto!(d,s)
end
@inline function Base.copyto!(dest::Symbol,source::Symbol,D::DataMultiBlock)
    for I in eachblock(D)
        d = getproperty(D[I],dest)
        # d .= getproperty(D[I],source)
        s = getproperty(D[I],source)
        copyto!(d,s)
    end
end




"""
"""
function relerr(D::newLocalDataBlock{TT}) where {TT}
    u = getarray(D,:u)
    v = getarray(D,:uₙ₊₁)
    cache = getarray(D,:cache)
    cache = u - v
    D.SC.Δu = innerprod(cache,cache,D.innerprod)/innerprod(u,u,D.innerprod)
end


"""
    applyParallelPenalty
"""
applyParallelPenalty(D::newLocalDataBlock) = applyParallelPenalty!(D.uₙ₊₁,D.u,D.SC.Δt,D.PGrid)
"""
    applyParallelPenalties
"""
function applyParallelPenalties(DB::DataMultiBlock)
    for I in eachblock(DB)
        applyParallelPenalty(DB[I])
    end
end

