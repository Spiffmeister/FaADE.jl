

"""
    addSource!
Adds PDE forcing term
"""
function addSource! end
function addSource!(S::SourceTerm{Nothing},tmp...) end
function addSource!(S::SourceTerm{FT},u::AT,grid::LocalGridType,t::TT,Δt::TT,θ::TT) where {TT,AT<:AbstractArray{TT},FT<:Function}
    f = S.source
    for I in eachindex(grid)
        # u[I] = Δt*θ*f(grid[I]...,t+Δt)
        u[I] = u[I] + Δt*(1-θ)*f(grid[I],t) + Δt*θ*f(grid[I],t+Δt)
    end
    u
end
function addSource!(dest::Symbol,D::DataMultiBlock{TT},θ=TT(1)) where {TT}
    for I in eachblock(D)
        write = getproperty(D[I],dest)
        addSource!(D[I].source, write, D[I].grid, D[I].SC.t, D[I].SC.Δt,D[I].SC.θ)
    end
end


"""
    BoundaryConditions
Sets the value of the boundary.
"""
function setBoundaryCondition! end
function setBoundaryCondition!(B::BoundaryData{TT,1,Fn},Δt::TT,args...) where {TT,Fn<:Real}
    @. B.BufferRHS = Δt*B.RHS
end
function setBoundaryCondition!(B::BoundaryData{TT,1,Fn},Δt::TT,t::TT,θ::TT) where {TT,Fn<:Function}
    B.BufferRHS[1] = Δt*(1-θ)*B.RHS(t) + Δt*θ*B.RHS(t+Δt)
end
function setBoundaryCondition!(B::BoundaryData{TT,2,Fn},Δt::TT,t::TT,θ::TT) where {TT,Fn<:Function}
    for i = 1:B.n
        B.BufferRHS[i] = Δt*(1-θ)*B.RHS(B.X[i],t) + Δt*θ*B.RHS(B.X[i],t+Δt)
    end
end
"""
If no boundary condition set or the boundary condition is periodic or an interface, ignore this step
"""
function setBoundaryCondition!(BC::Nothing,args...) end
function setBoundaryCondition!(BC::InterfaceBoundaryData,args...) end
"""
Calling boundaries for data blocks
"""
function setBoundaryConditions!(D::LocalDataBlockType{TT,DIM}) where {TT,DIM}
    setBoundaryCondition!(D.boundary[1], D.SC.Δt, D.SC.t, D.SC.θ)
    setBoundaryCondition!(D.boundary[2], D.SC.Δt, D.SC.t, D.SC.θ)
    if DIM == 2
        setBoundaryCondition!(D.boundary[3], D.SC.Δt, D.SC.t, D.SC.θ)
        setBoundaryCondition!(D.boundary[4], D.SC.Δt, D.SC.t, D.SC.θ)
    end
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
function fillBuffer!(source,B::BoundaryData,args...) end
function fillBuffer!(source::AT,B::InterfaceBoundaryData,K::AT) where AT
    B.BufferIn .= source
end
"""
fillBuffer
"""
function fillBuffer end
function fillBuffer(source,B::Nothing,args...) end # no boundary data
function fillBuffer(source::Symbol,B::BoundaryData,DB::DataMultiBlock) end
function fillBuffer(source::Symbol,B::InterfaceBoundaryData,DB::DataMultiBlock)
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
function fillBuffers(B::BoundaryData,args...) end
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
Applies SATs to the given data block.
"""
function applySAT! end
"""
Apply the SAT
"""
@inline applySAT!(SAT::SimultanousApproximationTerm,dest::AT,K::KT,u::AT,mode::SATMode) where {AT,KT} = SAT(dest,K,u,mode)
@inline applySAT!(SAT::SimultanousApproximationTerm,dest::AT,K::KT,u::AT,buffer::AT,mode::SATMode) where {AT,KT} = SAT(dest,K,u,buffer,mode)
"""
Applying SATs in DataMode ignoring interface terms
"""
@inline function applySAT!(BC::Nothing,tmp...) end # Do not apply and boundary conditions
@inline function applySAT!(BC::InterfaceBoundaryData,dest,K,mode::SATMode{:DataMode}) end # Do nothing
@inline function applyCartesianSAT!(BC::InterfaceBoundaryData,dest,K,mode::SATMode{:DataMode}) end # Do nothing
@inline function applyCurvilinearSAT!(BC::InterfaceBoundaryData,dest,K,mode::SATMode{:DataMode}) end # Do nothing

"""
Decide which SAT to apply
"""
@inline function applyCartesianSAT!(BC::BoundaryData,dest::AT,K::AT,mode::SATMode{:DataMode}) where {AT}
    if typeof(BC.Boundary) <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_data!(dest,BC.BufferRHS,K,BC.Boundary)
    elseif typeof(BC.Boundary) <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_data!(dest,BC.BufferRHS,BC.Boundary)
    else
        applySAT!(BC.Boundary,dest,K,mode)
    end
end

"""
Decide which Cartesian SAT to apply
"""
# function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Cartesian,TT,VT,FT1,PT,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,PT,LAT} where {SIDE,AX}
#     SAT_Dirichlet_data!(dest,BC.BufferRHS,K,BC.Boundary)
# end
"""
Decide which curvilinear SAT to apply
"""
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,PT,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,PT,LAT} where {SIDE,AX}
    SAT_Dirichlet_data!(dest,BC.BufferRHS,K[AX],K[3],BC.Boundary)
end


"""
Applying SATs in SolutionMode
"""
@inline function applyCartesianSAT!(BC::InterfaceBoundaryData,dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {AT}
    # applySAT!(BC.Boundary,dest,K,source,BC.BufferIn,mode)
    SAT = BC.Boundary
    if typeof(SAT) <: SimultanousApproximationTerm{:Interface}
        SAT_Interface!(dest,source,K,BC.BufferIn,SAT,mode)
    elseif typeof(SAT) <: SimultanousApproximationTerm{:Periodic}
        SAT_Periodic!(dest,source,K,SAT)
    end
end
function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {TT,AT,DIM,FT,BCT}
    SAT = BC.Boundary
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_solution!(dest,source,K,SAT)
    elseif BCT <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_solution!(dest,source,K,SAT)
    else
        error("Not implemented")
    end
end



# function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Cartesian,TT,VT,FT1,PT,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,PT,LAT} where {SIDE,AXIS}#,BCT<:SAT_Dirichlet}
#     SAT_Dirichlet_solution!(dest,source,K[AXIS],BC.Boundary)
# end
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,PT,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,PT,LAT} where {SIDE,AXIS}#,BCT<:SAT_Dirichlet}
    SAT_Dirichlet_solution!(dest,source,K[AXIS],K[3],BC.Boundary)
end

# function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Periodic{TN,:Cartesian,TT,VT,FT1,FT2},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,FT2} where {SIDE,AXIS}
#     SAT_Periodic!(dest,source,K[AXIS],BC.Boundary)
# end
function applyCurvilinearSAT!(BC::InterfaceBoundaryData{TT,DIM,SAT_Periodic{TN,:Curvilinear,TT,VT,FT1,FT2},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,TN<:NodeType{SIDE,AXIS},VT,FT1,FT2} where {SIDE,AXIS}
    SAT_Periodic!(dest,source,K[AXIS],K[3],BC.Boundary)
end



function applySAT!(BC::BoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:ExplicitMode}) where {TT,AT,DIM,FT,BCT}
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_explicit!(dest,source,BC.BufferRHS,K,BC.Boundary,mode)
    end
end

"""
    applySATs
Solution and Data modes
"""
function applySATs end
function applySATs(dest::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT} #data mode
    applySAT!(D.boundary[1],  dest, D.K, mode)
    applySAT!(D.boundary[2],  dest, D.K, mode)
end
function applySATs(dest::VT,source::VT,D::newLocalDataBlock{TT,1,VT},mode) where {TT,VT} #solution mode
    applySAT!(D.boundary[1],  dest, source, D.K, mode)
    applySAT!(D.boundary[2],  dest, source, D.K, mode)
end
"""
    applySATs
applySATs for 2D local block
"""
function applySATs(dest::AT,D::newLocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT}
    applyCurvilinearSAT!(D.boundary[1],    dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[2],    dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[3],    dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[4],    dest, D.K, mode)
end
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT}
    applyCurvilinearSAT!(D.boundary[1],    dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[2],    dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[3],    dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[4],    dest, source, D.K, mode)
end

function applySATs(dest::AT,D::newLocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT}
    applyCartesianSAT!(D.boundary[1],    dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[2],    dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[3],    dest, D.K[2], mode)
    applyCartesianSAT!(D.boundary[4],    dest, D.K[2], mode)
end
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT}
    applyCartesianSAT!(D.boundary[1],    dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[2],    dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[3],    dest, source, D.K[2], mode)
    applyCartesianSAT!(D.boundary[4],    dest, source, D.K[2], mode)
end




function applyCurvilinearSATs(dest::AT,D::newLocalDataBlock{TT,2,COORD,AT},mode) where {TT,COORD,AT}
    applySAT!(D.boundary[1],    dest,   source, D.K[1], D.K[3], mode)
    applySAT!(D.boundary[2],    dest,   source, D.K[1], D.K[3], mode)
    applySAT!(D.boundary[3],    dest,   source, D.K[2], D.K[3], mode)
    applySAT!(D.boundary[4],    dest,   source, D.K[2], D.K[3], mode)
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
@inline function setDiffusionCoefficient!(k::TT,K::AbstractArray{TT},grid::GridType) where {TT<:Real} end
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
    for I in eachindex(grid)
        K[I] = κ(grid[I]...)
    end
    # for i = 1:grid.nx
    #     for j = 1:grid.ny
    #         K[i,j] = κ(grid[i,j]...)
    #     end
    # end
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
function copyUtoSAT!(D::DataMultiBlock{TT,DIM,DT},order::Int) where {TT,DIM,DT<:LocalDataBlockType}
    copyUtoSAT!(D.Data.boundary,D.Data.u,order)
end



"""
    muladd!(dest,source,α,β)
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
innerprod(u::AT,J::AT,v::AT,IP::innerH{TT}) where {TT,AT<:AbstractArray{TT}} = IP(u,J,v)
function innerprod(u::Symbol,v::Symbol,D::newLocalDataBlock{TT}) :: TT where {TT}
    V = getarray(D,v)
    U = getarray(D,u)
    IP = getproperty(D,:innerprod)
    J = D.grid.J
    ret = innerprod(U,J,V,IP)
    # ret = innerprod(U,V,IP)
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
@inline applyParallelPenalty(D::newLocalDataBlock) = applyParallelPenalty!(D.uₙ₊₁,D.u,D.SC.Δt,D.PGrid)
"""
    applyParallelPenalties
"""
function applyParallelPenalties(DB::DataMultiBlock)
    for I in eachblock(DB)
        applyParallelPenalty(DB[I])
    end
end





#=
"""
    kIminusB!
Compute ``I-BBᵀ/||B||²`` for each point in the grid
"""
function IminusB!(MagField::Function,t::TT) where {TT,AT::AbstractArray{TT}}
    for I in eachindex(G)
        B = MagField(grid[I],t)
        tmp = I - B*B'/norm(B)^2
    end
end
=#


# K[I] = κ(grid[I],t)
# K[I] += IminusB(MagField,t,grid[I])


