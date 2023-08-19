
#=
"""
    setBoundary!
Sets the value of the boundary.
"""
function setBoundary! end
function setBoundary!(RHS::F,Bound::AT,grid::Vector{T},n::Int,t::T,Δt::T) where {F,AT,T}
    for i = 1:n
        Bound[i] = Δt*RHS(grid[i],t)
    end
end
function setBoundary!(RHS::F,Bound::AT,t::T,Δt::T) where {F,AT,T}
    Bound[1] = Δt*RHS(t)
end


function setBoundaries(B::SATBoundaries,G::GridType,D::newDataBlockType,t::TT,Δt::TT) where TT
    # setBoundary!(B.Boundary1,)
end
=#

#=
"""
    setCoefficient!
Sets the diffusion coefficient
"""
function setCoefficient! end
function setCoefficient!(K::Function,κ::AbstractArray,grid::Grid1D)
    for i = 1:grid.n
        κ[i] = K(grid.grid[i])
    end
end
function setCoefficient!(K::Function,κ::AbstractArray,grid::Grid2D)
    for i = 1:grid.nx
        for j = 1:grid.ny
            κ[i,j] = K(grid.gridx[i],grid.gridy[j])
        end
    end
end
function setCoefficient!(DC::DiffusionCoefficient{F},κ::AbstractArray,grid::Grid2D) where {F<:Function}
    setCoefficient!(DC.coeff,κ,grid)
end
function setCoefficient!(DC::DiffusionCoefficient{TT},κ::AbstractArray{TT},grid::GridType) where {TT}
    κ .= DC.coeff
end
=#


"""
    addSource!
"""
function addSource! end
function addSource!(S::SourceTerm{TT},u::AbstractArray{TT},grid::Grid1D{TT},t::TT,Δt::TT) where TT
    u .+= Δt*F.(grid.grid,t)
end
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
function addSource!(S::SourceTerm{F},D::newDataBlock,G::GridType,t::TT,Δt::TT) where {TT,F<:Function}
    for (LDB,LG) in zip(D.DataBlock,G.Grid)
        addSource(S.F,D.DataBlock,G.Grid,t,Δt)
    end
end
function addSource!(S::SourceTerm{Nothing},D::newDataBlock,grid::GridType,t::TT,Δt::TT) where TT end




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
function copyUtoSAT!(D::newDataBlock{TT,DIM,DT},order::DerivativeOrder) where {TT,DIM,DT<:newLocalDataBlockType}
    copyUtoSAT!(D.Data.boundary,D.Data.u,GetOrder(order))
end