

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
function setBoundaryCondition!(BC::Nothing,args...) end
function setBoundaryCondition!(BC::InterfaceBoundaryData,args...) end
"""
    setBoundaryConditions!(D::LocalDataBlockType{TT,DIM}) where {TT,DIM}

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
    setBoundaryConditions!(D::DataMultiBlock)

Calling boundaries from multiblocks
"""
function setBoundaryConditions!(D::DataMultiBlock)
    for I in eachblock(D)
        setBoundaryConditions!(D[I])
    end
end



"""
Move interface data from `getproperty(DataBlock,source)` into internal buffers `DataBlock.InterfaceBoundaryData.BufferIn`
"""
function _filllocalBuffer(source::Symbol,B::BoundaryData,D::LocalDataBlock) end
function _filllocalBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},D::LocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD,BCT<:SAT_Periodic}
    cache = getarray(D,source)
    BufferIn = B.BufferIn

    side = B.BoundaryOperator.side
    order = B.BoundaryOperator.order

    if side == Left
        BufferIn .= cache[end-order+1:end,:]
    elseif side == Right
        BufferIn .= cache[1:order,:]
    elseif side == Up
        BufferIn .= cache[:,end-order+1:end]
    elseif side == Down
        BufferIn .= cache[:,1:order]
    end
end
"""
Fill `InterfaceBoundaryData.BufferOut` with data from this block for sending to `InterfaceBoundaryData.BufferIn`
"""
function _filllocalBuffer(source::Symbol,B::InterfaceBoundaryData{TT,DIM,BCT},D::LocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD,BCT<:SAT_Interface}
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
"""
Loops over all `DataBlock.boundary` objects in a single `DataBlock`
"""
function _fillLocalBuffers(source::Symbol,D::LocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD}
    _filllocalBuffer(source,D.boundary[Left],D)
    _filllocalBuffer(source,D.boundary[Right],D)
    if DIM == 2
        _filllocalBuffer(source,D.boundary[Up],D)
        _filllocalBuffer(source,D.boundary[Down],D)
    end
end
"""
Loops over all `DataBlock`'s in a `DataMultiBlock`
"""
function _fillLocalBuffers(source::Symbol,D::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(D)
        _fillLocalBuffers(source,D[I])
    end
end


"""
Moves a `InterfaceBoundaryData.BufferOut` to the relevant `InterfaceBoundaryData.BufferOut` neighbouring `DataBlock.BufferIn`
"""
function _tradeBuffer!(B::BoundaryData,DB) end
function _tradeBuffer!(B::InterfaceBoundaryData{TT,DIM,BCT,AT},DB) where {TT,DIM,AT,BCT<:SAT_Periodic} end
function _tradeBuffer!(B::InterfaceBoundaryData{TT,DIM,BCT,AT},DB) where {TT,DIM,AT,BCT<:SAT_Interface}
    Myjoint = B.OutgoingJoint
    Tojoint = B.IncomingJoint

    MyBuffer = B.BufferOut :: AT # Get this buffer
    ToBuffer = DB[Myjoint.index].boundary[Tojoint.side].BufferIn :: AT # Get the buffer we are writing to

    if (Myjoint.side == _flipside(Tojoint.side)) 
        # if they are the same dimension and they match we can just write to array
        @. ToBuffer = MyBuffer
    elseif Myjoint.side == Tojoint.side
        # if they are the same dimension but they don't match we must reverse
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
            reverse!(ToBuffer,dims=mod1(typeof(Tojoint.side).parameters[2]+1,2))
        end
    end
end
"""
Loop over each `DataBlock.boundary`
"""
function _tradeBuffers!(D::LocalDataBlock{TT,DIM,COORD,AT,KT,DCT,GT,BT},DB::DataMultiBlock{TT,DIM}) where {TT,DIM,COORD,AT,KT,DCT,GT,BT}
    boundary = D.boundary :: BT
    # @show BT
    # boundaryLeft = boundary[Left]
    _tradeBuffer!(boundary[Left],DB)
    _tradeBuffer!(boundary[Right],DB)
    if DIM == 2
        _tradeBuffer!(boundary[Down],DB)
        _tradeBuffer!(boundary[Up],DB)
    end
    # @show D.boundary
    # for J in eachindex(D.boundary)
    #     _tradeBuffer!(D.boundary[J],DB)
    # end
    # map(B->_tradeBuffer!(B,DB),D.boundary)
end
"""
Loop over each `DataMultiBlock.DataBlock`
"""
function _tradeBuffers!(DB::DataMultiBlock{TT,DIM,COORD}) where {TT,DIM,COORD}
    for I in eachblock(DB)
        _tradeBuffers!(DB[I],DB)
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
    elseif typeof(BC.BoundaryOperator) <: SimultanousApproximationTerm{:Robin}
        SAT_Robin_data!(dest,BC.BufferRHS,BC.BoundaryOperator)
    else
        applySAT!(BC.BoundaryOperator,dest,K,mode)
    end
end

"""
Decide which curvilinear SAT to apply
"""
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,LAT} where {SIDE,AX}
    SAT_Dirichlet_data!(dest,BC.BufferRHS,K[AX],K[3],BC.BoundaryOperator)
end
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Neumann{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},VT,FT1,LAT} where {SIDE,AX}
    SAT_Neumann_data!(dest,BC.BufferRHS,BC.BoundaryOperator)
end
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Robin{TN,:Curvilinear,TT,FT1,LAT},AT},dest::AT,K::KT,mode::SATMode{:DataMode}) where {AT,KT,TT,DIM,FT,TN<:NodeType{SIDE,AX},FT1,LAT} where {SIDE,AX}
    SAT_Robin_data!(dest,BC.BufferRHS,BC.BoundaryOperator)
end

"""
Applying periodic or interface SATs in `SolutionMode`
"""
function applyCartesianSAT!(BC::InterfaceBoundaryData,dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {AT}
    # applySAT!(BC.BoundaryOperator,dest,K,source,BC.BufferIn,mode)
    SAT = BC.BoundaryOperator
    if typeof(SAT) <: SimultanousApproximationTerm{:Interface}
        SAT_Interface!(dest,source,K,BC.BufferIn,SAT,mode)
    elseif typeof(SAT) <: SimultanousApproximationTerm{:Periodic}
        SAT_Periodic!(dest,source,K,SAT)
    elseif typeof(SAT) <: SimultanousApproximationTerm{:Robin}
        SAT_Robin_solution!(dest,source,K,SAT)
    else
        error("Not implemented")
    end
end
"""
Apply Dirichlet, Neumann, Robin SATs in `SolutionMode`
"""
function applyCartesianSAT!(BC::BoundaryData{TT,DIM,FT,BCT},dest::AT,source::AT,K::AT,mode::SATMode{:SolutionMode}) where {TT,AT,DIM,FT,BCT}
    SAT = BC.BoundaryOperator
    if BCT <: SimultanousApproximationTerm{:Dirichlet}
        SAT_Dirichlet_solution!(dest,source,K,SAT)
    elseif BCT <: SimultanousApproximationTerm{:Neumann}
        SAT_Neumann_solution!(dest,source,K,SAT)
    elseif BCT <: SimultanousApproximationTerm{:Robin}
        SAT_Robin_solution!(dest,source,K,SAT)
    else
        error("Not implemented")
    end
end


"""
    applyCurvilinearSAT!
Apply SATs for curvilinear grids
"""
function applyCurvilinearSAT! end
"""
Dirichlet curvilinear SAT
"""
function applyCurvilinearSAT!(BC::BoundaryData{TT,DIM,FT,SAT_Dirichlet{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,LAT} where {SIDE,AXIS}#,BCT<:SAT_Dirichlet}
    SAT_Dirichlet_solution!(dest,source,K[AXIS],K[3],BC.BoundaryOperator)
end
"""
Neumann curvilinear SAT
"""
function applyCurvilinearSAT!(BD::BoundaryData{TT,DIM,FT,SAT_Neumann{TN,:Curvilinear,TT,VT,FT1,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},VT,FT1,LAT} where {SIDE,AXIS}
    SAT_Neumann_solution!(dest,source,K[AXIS],K[3],BD.BoundaryOperator)
end
"""
Neumann curvilinear SAT
"""
function applyCurvilinearSAT!(BD::BoundaryData{TT,DIM,FT,SAT_Robin{TN,:Curvilinear,TT,FT1,LAT},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,FT,TN<:NodeType{SIDE,AXIS},FT1,LAT} where {SIDE,AXIS}
    SAT_Robin_solution!(dest,source,K[AXIS],K[3],BD.BoundaryOperator)
end
"""
Interface curvilinear SAT
"""
function applyCurvilinearSAT!(BC::InterfaceBoundaryData{TT,DIM,SAT_Periodic{TN,:Curvilinear,TT,VT,FT1,FT2},AT},dest::AT,source::AT,K::KT,mode::SATMode{:SolutionMode}) where {TT,AT,KT<:Vector{AT},DIM,TN<:NodeType{SIDE,AXIS},VT,FT1,FT2} where {SIDE,AXIS}
    SAT_Periodic!(dest,source,K[AXIS],K[3],BC.BoundaryOperator)
end
"""
Interface curvilinear SAT
"""
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
function applySATs(dest::VT,D::LocalDataBlock{TT,1,MET,VT},mode) where {TT,VT<:Vector{TT},MET} #data mode
    applyCartesianSAT!(D.boundary[Left],  dest, D.K, mode)
    applyCartesianSAT!(D.boundary[Right],  dest, D.K, mode)
end
function applySATs(dest::VT,source::VT,D::LocalDataBlock{TT,1,MET,VT},mode) where {TT,VT<:Vector{TT},MET} #solution mode
    applyCartesianSAT!(D.boundary[Left],  dest, source, D.K, mode)
    applyCartesianSAT!(D.boundary[Right],  dest, source, D.K, mode)
end
"""
    applySATs
applySATs for 2D local block
"""
function applySATs(dest::AT,D::LocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT} #data mode
    applyCurvilinearSAT!(D.boundary[Left],  dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Right], dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Up],    dest, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Down],  dest, D.K, mode)
end
function applySATs(dest::AT,source::AT,D::LocalDataBlock{TT,2,:Variable,AT},mode) where {TT,AT} #solution mode
    applyCurvilinearSAT!(D.boundary[Left],  dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Right], dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Up],    dest, source, D.K, mode)
    applyCurvilinearSAT!(D.boundary[Down],  dest, source, D.K, mode)
end

function applySATs(dest::AT,D::LocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT} #data mode
    applyCartesianSAT!(D.boundary[Left],    dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Right],   dest, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Up],      dest, D.K[2], mode)
    applyCartesianSAT!(D.boundary[Down],    dest, D.K[2], mode)
end
function applySATs(dest::AT,source::AT,D::LocalDataBlock{TT,2,:Constant,AT},mode) where {TT,AT} #solution mode
    applyCartesianSAT!(D.boundary[Left],    dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Right],   dest, source, D.K[1], mode)
    applyCartesianSAT!(D.boundary[Up],      dest, source, D.K[2], mode)
    applyCartesianSAT!(D.boundary[Down],    dest, source, D.K[2], mode)
end




function applyCurvilinearSATs(dest::AT,D::LocalDataBlock{TT,2,COORD,AT},mode) where {TT,COORD,AT}
    applySAT!(D.boundary[Left], dest,   source, D.K[1], D.K[3], mode)
    applySAT!(D.boundary[Right],dest,   source, D.K[1], D.K[3], mode)
    applySAT!(D.boundary[Up],   dest,   source, D.K[2], D.K[3], mode)
    applySAT!(D.boundary[Down], dest,   source, D.K[2], D.K[3], mode)
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
    setDiffusionCoefficient!(D::LocalDataBlock{TT,1})
Set diffusion coefficient for 1D grids
"""
setDiffusionCoefficient!(D::LocalDataBlock{TT,1}) where TT = setDiffusionCoefficient!(D.κ,D.K,D.grid)
"""
    setDiffusionCoefficient!(D::LocalDataBlock{TT,2})
Set diffusion coefficient for 2D grids
"""
function setDiffusionCoefficient!(D::LocalDataBlock{TT,2}) where TT
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
function muladd!(dest::Symbol,source::Symbol,D::LocalDataBlock,α,β)
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
function setValue(dest::Symbol,source::Symbol,D::LocalDataBlock,α)
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
function innerprod(u::Symbol,v::Symbol,D::LocalDataBlock{TT,DIM,:Constant}) :: TT where {TT,DIM}
    V = getarray(D,v)
    U = getarray(D,u)
    IP = getproperty(D,:innerprod)
    ret = innerprod(U,V,IP)
    return ret
end
function innerprod(u::Symbol,v::Symbol,D::LocalDataBlock{TT,DIM,:Variable}) :: TT where {TT,DIM}
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
@inline function Base.copyto!(dest::Symbol,source::Symbol,D::LocalDataBlock)
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
Compute relative error
"""
function relerr(D::LocalDataBlock{TT}) where {TT}
    u = getarray(D,:u)
    v = getarray(D,:uₙ₊₁)
    cache = getarray(D,:cache)
    cache = u - v
    D.SC.Δu = innerprod(cache,cache,D.innerprod)/innerprod(u,u,D.innerprod)
end


"""
    applyParallelPenalty
"""
@inline applyParallelPenalty(D::LocalDataBlock) = applyParallelPenalty!(D.uₙ₊₁,D.u,D.SC.Δt,D.PGrid)
"""
    applyParallelPenalties
"""
function applyParallelPenalties(DB::DataMultiBlock)
    for I in eachblock(DB)
        applyParallelPenalty(DB[I])
    end
end


"""
    Updates the `BivariateCHSInterpolation` object from `CubicHermiteSpline.jl`
"""
function _updateCHSinterp(D::LocalDataBlock{TT,2,COORD,AT}) where {TT,COORD,AT}
    
    rₖ = D.rₖ   :: AT
    dₖ = D.dₖ   :: AT

    u = D.uₙ₊₁

    Dx = D.Derivative.DO[1]
    Dy = D.Derivative.DO[2]
    Interp = D.Parallel.Interpolant

    grid = D.grid

    nx = Dx.n
    ny = Dy.n
    #overwrite rₖ
    QX = grid.qx    :: AT
    RX = grid.rx    :: AT
    QY = grid.qy    :: AT
    RY = grid.ry    :: AT
    
    # ∂u/∂x = ∂u/∂q ∂q/∂x + ∂u/∂r ∂r/∂x ≈ q_x D_q u + r_x D_r u
    for (cache,qx,U) in zip(eachcol(rₖ),eachcol(QX),eachcol(u))
        D₁!(cache, qx, U, nx, Dx.Δx, Dx.order, TT(0))  ## Need q_x
    end
    for (cache,rx,U) in zip(eachrow(rₖ),eachrow(RX),eachrow(u))
        D₁!(cache, rx, U, ny, Dy.Δx, Dx.order, TT(1))  ## Need r_x
    end
    #overwrite dₖ
    # ∂u/∂y = ∂u/∂y ∂y/∂q + ∂u/∂y ∂y/∂r ≈ q_y D_q u + r_y D_r u
    for (cache,qy,U) in zip(eachcol(dₖ),eachcol(QY),eachcol(u))
        D₁!(cache, qy, U, nx, Dx.Δx, Dy.order, TT(0))  ## Need r_y
    end
    for (cache,ry,U) in zip(eachrow(dₖ),eachrow(RY),eachrow(u))
        D₁!(cache, ry, U, ny, Dy.Δx, Dy.order, TT(1))  ## Need r_x
    end
    

    # If the domain is periodic i.e. for a hollow torus then we have to remove
    #   some elements
    if nx != size(u)[1]
        nx = Dx.n-1
    end
    if ny != size(u)[2]
        ny = Dy-1
    end

    u_resize = view(u, 1:nx, 1:ny)
    rₖ_resize = view(rₖ, 1:nx, 1:ny)
    dₖ_resize = view(dₖ, 1:nx, 1:ny)

    if COORD == :Constant # Cartesian coordinates do not require Jacobian
        for I in eachindex(Interp.z)
            Interp.z[I] = u_resize[I]
            Interp.dzdx[I] = rₖ_resize[I]
            Interp.dzdy[I] = dₖ_resize[I]
        end
    else
        grid = D.grid #TODO: Function barrier probably required for speed/allocations
        J = grid.J
        # J_resize = view(J, 1:nx, 1:ny)

        for I in eachindex(Interp.z)
            Interp.z[I] = u_resize[I]
            Interp.dzdx[I] = rₖ_resize[I] #/ J_resize[I]
            Interp.dzdy[I] = dₖ_resize[I] #/ J_resize[I]
        end
    end

end
function _updateCHSinterp(DB::DataMultiBlock)
    _fillLocalBuffers(:uₙ₊₁,DB)
    _tradeBuffers!(DB)
    _average_boundaries(DB)

    for I in eachblock(DB)
        _updateCHSinterp(DB[I])
    end
end

"""
    _setglobalu!
"""
function _setglobalu! end
function _setglobalu!(DB::DataMultiBlock,uglobal::Vector{AT}) where AT
    #TODO: CURRENTLY HACKED - FIX IF CHS INTERP WORKS
    for I in eachblock(DB)
        DB[I].uₙ₊₁ .= uglobal[I]
    end
end

"""
    setglobalu!
"""
function _setglobalu!(uglobal::Vector{AT},DB::DataMultiBlock) where AT

    # To perform the averaging at boundaries we'll trade the buffers then average them in-block
    _fillLocalBuffers(:uₙ₊₁,DB)
    _tradeBuffers!(DB)
    _average_boundaries(DB)

    # Now we can call them into the global space
    for I in eachblock(DB)
        uglobal[I] .= getarray(DB[I],:uₙ₊₁)
    end

end

function _average_boundaries(DB::DataMultiBlock{TT,DIM}) where {TT,DIM}
    for I in eachblock(DB)
        _average_local_buffer(DB[I])
    end
end
function _average_local_buffer(D::LocalDataBlock{TT,DIM,COORD}) where {TT,DIM,COORD}
    _average_boundary_data(D.boundary[Left],D)
    _average_boundary_data(D.boundary[Right],D)
    if DIM == 2
        _average_boundary_data(D.boundary[Up],D)
        _average_boundary_data(D.boundary[Down],D)
    end
end
function _average_boundary_data(B::BoundaryData,D::LocalDataBlock) end
function _average_boundary_data(B::InterfaceBoundaryData{TT,DIM,BCT,AT},D::LocalDataBlock) where {TT,DIM,AT,BCT<:SAT_Periodic} end
function _average_boundary_data(B::InterfaceBoundaryData{TT,DIM,BCT,AT},D::LocalDataBlock) where {TT,DIM,AT,BCT<:SAT_Interface}
    MyBuffer = B.BufferOut
    IncomingBuffer = B.BufferIn

    Myjoint = B.OutgoingJoint #What side is it on for this block

    Myside = Myjoint.side

    if Myside == Left
        uₙ₊₁ = @view D.uₙ₊₁[1,:]
        @. uₙ₊₁ = (MyBuffer[1,:] + IncomingBuffer[1,:])/2
    elseif Myside == Right
        uₙ₊₁ = @view D.uₙ₊₁[end,:]
        @. uₙ₊₁ = (MyBuffer[1,:] + IncomingBuffer[1,:])/2
    elseif Myside == Up
        uₙ₊₁ = @view D.uₙ₊₁[:,end]
        @. uₙ₊₁ = (MyBuffer[:,1] + IncomingBuffer[:,1])/2
    else
        uₙ₊₁ = @view D.uₙ₊₁[:,1]
        @. uₙ₊₁ = (MyBuffer[:,1] + IncomingBuffer[:,1])/2
    end

end
