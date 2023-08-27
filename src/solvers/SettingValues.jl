

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
function BoundaryConditions! end
function BoundaryConditions!(RHS::Function,Bound::AT,t::T,Δt::T) where {AT,T}
    Bound[1] = Δt*RHS(t)
end
function BoundaryConditions!(RHS::Function,Bound::AT,grid::Vector{T},n::Int,t::T,Δt::T) where {AT,T}
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end
function BoundaryConditions(D::newLocalDataBlockType,G::Grid1D,t::TT,Δt::TT) where TT
    BoundaryConditions!(D.boundary.Left.RHS,    D.boundary.RHS_Left, t,Δt)
    BoundaryConditions!(D.boundary.Right.RHS,   D.boundary.RHS_Right,t,Δt)
end
function BoundaryConditions(D::newLocalDataBlockType,G::Grid2D,t::TT,Δt::TT) where TT
    BoundaryConditions!(B.BoundaryLeft.RHS, D.Data.RHS_Left,G.nx,G.gridx,t,Δt)
    BoundaryConditions!(B.BoundaryRight.RHS,D.Data.RHS_Right,G.nx,G.gridx,t,Δt)
    BoundaryConditions!(B.BoundaryUp.RHS,   D.Data.RHS_Up,  G.ny,G.gridy,t,Δt)
    BoundaryConditions!(B.BoundaryDown.RHS, D.Data.RHS_Down,G.ny,G.gridy,t,Δt)
end
function BoundaryConditions(D::DataMultiBlock,G::GridType)
    for I in eachblock(D)
        BoundaryConditions(D[I],   G,  D.SC.t, D.SC.Δt)
    end
end

"""
    fillBuffers()
"""
function fillBuffers()
end



"""
    applySATs
"""
function applySATs(dest::Array{TT},D::LocalDataBlock{TT,DIM,NBLOCK},mode::SATMode) where {TT,DIM,NBLOCK}
    applySAT!(D.boundary.Left,   dest, D.boundary.RHS_Left,    D.K, mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.RHS_Right,   D.K, mode)
    if DIM == 2
        applySAT!(D.boundary.Up,     dest, D.boundary.RHS_Left,    D.K, mode)
        applySAT!(D.boundary.Down,   dest, D.boundary.RHS_Right,   D.K, mode)
    end
end
function applySATs(CG,D::DataMultiBlock{TT,DIM,NBLOCK},dest::Symbol,mode::SATMode) where {TT,DIM,NBLOCK}
    for i = 1:NBLOCK
        tmp = getproperty(CG,dest)
        applySATs(CG,DB[i],mode)

        # applySAT!(D[i].boundary.Left,   CG[i].b, D[i].boundary.RHS_Left,    D[i].K, mode)
        # applySAT!(D[i].boundary.Right,  CG[i].b, D[i].boundary.RHS_Right,   D[i].K, mode)
        # if DIM == 2
        #     applySAT!(D[i].boundary.Up,     CG[i].b, D[i].boundary.RHS_Left,    D[i].K, mode)
        #     applySAT!(D[i].boundary.Down,   CG[i].b, D[i].boundary.RHS_Right,   D[i].K, mode)
        # end
    end
end

function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,1},mode::SATMode) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.RHS_Left,    D.K, mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.RHS_Right,   D.K, mode)
end
function applySATs(dest::Array{TT},D::newLocalDataBlock{TT,2},mode::SATMode) where {TT}
    applySAT!(D.boundary.Left,   dest, D.boundary.RHS_Left,    D.K[1], mode)
    applySAT!(D.boundary.Right,  dest, D.boundary.RHS_Right,   D.K[1], mode)
    applySAT!(D.boundary.Up,     dest, D.boundary.RHS_Left,    D.K[2], mode)
    applySAT!(D.boundary.Down,   dest, D.boundary.RHS_Right,   D.K[2], mode)
end
function applySATs(dest::Symbol,D::DataMultiBlock,mode::SATMode)
    for I in eachblock(D)
        write = getproperty(D[I],dest)
        applySATs(write,D[I],mode)
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
@inline function muladd!(dest::Symbol,source::Symbol,D::DataMultiBlock{TT};α=TT(1),β=TT(1)) where TT

    for I in eachblock(D)
        W = getproperty(D[I],dest)
        A = getproperty(D[I],source)

        @. W = α*W + β*A
    end
end
"""
    setValue
"""
@inline function setValue(dest::Symbol,source::Symbol,D::newDataBlockType{TT},α=TT(1)) where TT

    for I in eachblock(D)
        W = getproperty(D[I],dest)
        R = getproperty(D[I],source)

        @. W = α*R
    end
end

"""
    innerprod
"""
function innerprod(u::Symbol,v::Symbol,D::DataMultiBlock{TT}) where TT
    local r::TT
    r = TT(0)
    for I in eachblock(D)
        r += D[I].innerprod(getproperty(D[I],u),getproperty(D[I],v))
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