

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
    setBoundaryCondition!(D.boundary[Left], D.SC.Δt, D.SC.t, D.SC.θ)
    setBoundaryCondition!(D.boundary[Right], D.SC.Δt, D.SC.t, D.SC.θ)
    if DIM == 2
        setBoundaryCondition!(D.boundary[Up], D.SC.Δt, D.SC.t, D.SC.θ)
        setBoundaryCondition!(D.boundary[Down], D.SC.Δt, D.SC.t, D.SC.θ)
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
function fillBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},DB::DataMultiBlock) where {TT,DIM,BCT<:SAT_Periodic}
    cache = getproperty(DB[B.Joint],source)
    copyto!(B.BufferIn,CartesianIndices(B.BufferIn),cache,B.I)
end
function fillBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},DB::DataMultiBlock) where {TT,DIM,BCT<:SAT_Interface}
    cache = getproperty(DB[B.Joint],source)
    K = getproperty(DB[B.Joint],:K)
    SAT = B.BoundaryOperator
    # Kq = K[SAT.axis]
    # Kqr = K[3]

    if DIM == 1
        SAT_Interface_cache!(B.BufferIn,cache,K,SAT)
    elseif DIM == 2
        if SAT.coordinates == :Cartesian
            SAT_Interface_cache!(B.BufferIn,cache,K[SAT.axis],SAT)
        elseif SAT.coordinates == :Curvilinear
            SAT_Interface_cache!(B.BufferIn,cache,K[SAT.axis],K[3],SAT)
        end
    end
end
function fillBuffer(source::Symbol,DB::DataMultiBlock,I::Int64,side::NodeType)
    BC = DB[I].boundary[side]
    fillBuffer(source,BC,DB)
end


function _filllocalBuffer(source::Symbol,B::BoundaryData,D::newLocalDataBlock) end
function _filllocalBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},D::newLocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD,BCT<:SAT_Periodic} end
function _filllocalBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},D::newLocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD,BCT<:SAT_Interface}
    cache = getproperty(D,source)
    K = getproperty(D,:K)
    SAT = B.BoundaryOperator

    if DIM == 1
        SAT_Interface_cache!(B.BufferIn,cache,D.K,B.BoundaryOperator)
    elseif DIM == 2
        if B.BoundaryOperator.coordinates == :Cartesian
            SAT_Interface_cache!(B.BufferOut,cache,D.K[SAT.axis],B.BoundaryOperator)
        elseif B.BoundaryOperator.coordinates == :Curvilinear
            SAT_Interface_cache!(B.BufferOut,cache,D.K[SAT.axis],D.K[3],B.BoundaryOperator)
        end
    end
end
function _fillLocalBuffers(source::Symbol,D::newLocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD}
    _filllocalBuffer(source,D.boundary[Left],D)
    _filllocalBuffer(source,D.boundary[Right],D)
    if DIM == 2
        _filllocalBuffer(source,D.boundary[Up],D)
        _filllocalBuffer(source,D.boundary[Down],D)
    end
end
function _fillLocalBuffers(source::Symbol,D::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(D)
        _fillLocalBuffers(source,D[I])
    end
end


function _tradeBuffer!(B::BoundaryData,DB) end
function _tradeBuffer!(B::InterfaceBoundaryData,DB)
    side = B.BoundaryOperator.side
    # @show B.Joint

    Myjoint = B.OutgoingJoint
    Tojoint = B.IncomingJoint
    
    MyBuffer = B.BufferOut # Get this buffer
    ToBuffer = DB[Myjoint.index].boundary[Tojoint.side].BufferIn # Get the buffer we are writing to

    # @show Myjoint
    # @show Tojoint
    # @show DB[Myjoint.index].boundary[Tojoint.side].BufferIn

    if (Myjoint.side == _flipside(Tojoint.side)) 
        # if they are the same dimension and they match we can just write to array
        ToBuffer .= MyBuffer
    elseif Myjoint.side == Tojoint.side
        # if they are the same dimension but they don't match we must reverse
        # @show Myjoint, Tojoint
        # @show Tojoint.side
        ToBuffer .= MyBuffer
        reverse!(ToBuffer,dims=mod1(typeof(Tojoint.side).parameters[2]+1,2))
        # reverse!(ToBuffer,dims=2)
    else
        # if they are not the the same dimension we need to rearrange things
        for (BI,BO) in zip(eachrow(ToBuffer),eachcol(MyBuffer))
            BI .= BO
        end
        if typeof(Myjoint.side).parameters[1] != typeof(Tojoint.side).parameters[1]
            # If they are not the same axis we need to reverse
            # @show Myjoint, Tojoint
            # @show ToBuffer
            reverse!(ToBuffer,dims=mod1(typeof(Tojoint.side).parameters[2]+1,2))
        end
    end
end
function _tradeBuffers!(D,DB)
    for J in eachindex(D.boundary)
        _tradeBuffer!(D.boundary[J],DB)
    end
end
function _tradeBuffers!(DB::DataMultiBlock{TT,DIM,COORD}) where {TT,DIM,COORD}
    for I in eachblock(DB)
        _tradeBuffers!(DB[I],DB)
    end
end

"""
    fillBuffers
"""
function fillBuffers end
function fillBuffers(B::BoundaryData,args...) end
function fillBuffer(source::Symbol,DB::DataMultiBlock{TT,DIM},I::Int64) where {TT,DIM}
    D = DB[I]
    fillBuffer(source,D.boundary[Left],DB)
    fillBuffer(source,D.boundary[Right],DB)
    if DIM == 2
        fillBuffer(source,D.boundary[Up],DB)
        fillBuffer(source,D.boundary[Down],DB)
    end
end
function fillBuffers(source::Symbol,DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(DB)
        fillBuffer(source,DB,I)
    end
end



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
    if typeof(BC.BoundaryOperator) <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_data!(dest,BC.BufferRHS,K,BC.BoundaryOperator)
    elseif typeof(BC.BoundaryOperator) <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_data!(dest,BC.BufferRHS,BC.BoundaryOperator)
    else
        applySAT!(BC.BoundaryOperator,dest,K,mode)
    end
end

"""
Decide which Cartesian SAT to apply
"""
# function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Cartesian,TT,VT,FT1,PT,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,PT,LAT} where {SIDE,AX}
#     SAT_Dirichlet_data!(dest,BC.BufferRHS,K,BC.BoundaryOperator)
# end
"""
Decide which curvilinear SAT to apply
"""
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,LAT} where {SIDE,AX}
    SAT_Dirichlet_data!(dest,BC.BufferRHS,K[AX],K[3],BC.BoundaryOperator)
end


"""
Applying SATs in SolutionMode
"""
function applyCartesianSAT!(BC::InterfaceBoundaryData,dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {AT}
    # applySAT!(BC.BoundaryOperator,dest,K,source,BC.BufferIn,mode)
    SAT = BC.BoundaryOperator
    if typeof(SAT) <: SimultanousApproximationTerm{:Interface}
        SAT_Interface!(dest,source,K,BC.BufferIn,SAT,mode)
    elseif typeof(SAT) <: SimultanousApproximationTerm{:Periodic}
        SAT_Periodic!(dest,source,K,SAT)
    end
end
function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {TT,AT,DIM,FT,BCT}
    SAT = BC.BoundaryOperator
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_solution!(dest,source,K,SAT)
    elseif BCT <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_solution!(dest,source,K,SAT)
    else
        error("Not implemented")
    end
end


"""
    applyCurvilinearSAT!
Apply SATs for curvilinear grids
"""
function applyCurvilinearSAT! end
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,LAT} where {SIDE,AXIS}#,BCT<:SAT_Dirichlet}
    SAT_Dirichlet_solution!(dest,source,K[AXIS],K[3],BC.BoundaryOperator)
end
function applyCurvilinearSAT!(BC::InterfaceBoundaryData{TT,DIM,SAT_Periodic{TN,:Curvilinear,TT,VT,FT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,TN<:NodeType{SIDE,AXIS},VT,FT} where {SIDE,AXIS}
    SAT_Periodic!(dest,source,K[AXIS],K[3],BC.BoundaryOperator)
end
function applyCurvilinearSAT!(BC::InterfaceBoundaryData{TT,DIM,SAT_Interface{TN,:Curvilinear,TT,VT,FT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,TN<:NodeType{SIDE,AXIS},VT,FT} where {SIDE,AXIS}
    SAT_Interface!(dest,source,K[AXIS],K[3],BC.BufferIn,BC.BoundaryOperator,mode)
end


function applySAT!(BC::BoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:ExplicitMode}) where {TT,AT,DIM,FT,BCT}
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_explicit!(dest,source,BC.BufferRHS,K,BC.BoundaryOperator,mode)
    end
end

"""
    applySATs
Solution and Data modes
"""
function applySATs end
function applySATs(dest::VT,D::newLocalDataBlock{TT,1,MET,VT},mode) where {TT,VT<:Vector{TT},MET} #data mode
    applyCartesianSAT!(D.boundary[Left],  dest, D.K, mode)
    applyCartesianSAT!(D.boundary[Right],  dest, D.K, mode)
end
function applySATs(dest::VT,source::VT,D::newLocalDataBlock{TT,1,MET,VT},mode) where {TT,VT<:Vector{TT},MET} #solution mode
    applyCartesianSAT!(D.boundary[Left],  dest, source, D.K, mode)
    applyCartesianSAT!(D.boundary[Right],  dest, source, D.K, mode)
end
"""
    applySATs
applySATs for 2D local block
"""
function applySATs(dest::AT,D::newLocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT} #data mode
    applyCurvilinearSAT!(D.boundary[Left],  dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Right], dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Up],    dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Down],  dest, D.K, mode)
end
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT} #solution mode
    applyCurvilinearSAT!(D.boundary[Left],  dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Right], dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Up],    dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Down],  dest, source, D.K, mode)
end

function applySATs(dest::AT,D::newLocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT} #data mode
    applyCartesianSAT!(D.boundary[Left],    dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Right],   dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Up],      dest, D.K[2], mode)
    applyCartesianSAT!(D.boundary[Down],    dest, D.K[2], mode)
end
function applySATs(dest::AT,source::AT,D::newLocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT} #solution mode
    applyCartesianSAT!(D.boundary[Left],    dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Right],   dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Up],      dest, source, D.K[2], mode)
    applyCartesianSAT!(D.boundary[Down],    dest, source, D.K[2], mode)
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
function innerprod(u::Symbol,v::Symbol,D::newLocalDataBlock{TT,DIM,:Constant}) :: TT where {TT,DIM}
    V = getarray(D,v)
    U = getarray(D,u)
    IP = getproperty(D,:innerprod)
    ret = innerprod(U,V,IP)
    return ret
end
function innerprod(u::Symbol,v::Symbol,D::newLocalDataBlock{TT,DIM,:Variable}) :: TT where {TT,DIM}
    V = getarray(D,v)
    U = getarray(D,u)
    IP = getproperty(D,:innerprod)
    J = D.grid.J
    ret = innerprod(U,J,V,IP)
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


