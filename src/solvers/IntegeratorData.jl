
abstract type IntegratorDataMultiBlock{dtype,DIM,atype} <: DataBlockType{dtype,DIM,atype} end

abstract type IntegratorDataBlock{dtype,DIM,atype} <: DataBlockType{dtype,DIM,atype} end

mutable struct ConjGradScalarBlock{TT}
    rnorm       :: TT
    unorm       :: TT
    dₖAdₖ       :: TT
    βₖ          :: TT
    αₖ          :: TT
    converged   :: Bool
    function ConjGradScalarBlock{TT}() where TT
        new{TT}(TT(0),TT(0),TT(0),TT(0),TT(0),true)
    end
end


"""
    ConjGradBlock
Data storage for the conjugate gradient method. See [`conj_grad!`](@ref)
"""
struct ConjGradBlock{TT,DIM,AT} <: IntegratorDataBlock{TT,DIM,AT}
    b   :: AT # b = uⁿ⁺¹ + F
    rₖ  :: AT # (uⁿ⁺¹ - Δt*uₓₓⁿ⁺¹) - b
    Adₖ :: AT # Adₖ = dₖ - Δt*D(dₖ)
    Drₖ :: AT # Drₖ = rₖ - Δt*D(rₖ)
    dₖ  :: AT # dₖ = -rₖ, -rₖ .+ βₖ*dₖ

    scalar :: ConjGradScalarBlock

    innerprod   :: innerH{TT,DIM,Vector{TT}}

    function ConjGradBlock{TT}(grid::GridType{TT,DIM},order::Union{Int,DerivativeOrder}) where {TT,DIM}

        if typeof(grid) <: Grid1D
            Δ = grid.Δx
            n = grid.n
        elseif typeof(grid) <: Grid2D
            Δ = min(grid.Δx,grid.Δy)
            n = (grid.nx,grid.ny)
        end

        b   = zeros(TT,n)
        rₖ  = zeros(TT,n)
        Adₖ = zeros(TT,n)
        Arₖ = zeros(TT,n)
        dₖ  = zeros(TT,n)

        innerprod = innerH(grid,order)

        scalars = ConjGradScalarBlock{TT}()

        new{TT,DIM, typeof(b)}(b, rₖ, Adₖ, Arₖ, dₖ, scalars, innerprod)
    end
end

struct ConjGradMultiBlock{TT,DIM,NBLOCK,AT} <: IntegratorDataBlock{TT,DIM,AT}
    Block :: Vector{ConjGradBlock}
    function ConjGradMultiBlock(G::GridMultiBlock{TT,DIM},DO::DerivativeOrder) where {TT,DIM}
        CGB = [ConjGradBlock{TT}(grid,GetOrder(DO)) for grid in G.Grids]
        new{TT,DIM,length(CBG),AT}(CGB)
    end
    function ConjGradMultiBlock(G::LocalGridType{TT,DIM},DO::DerivativeOrder) where {TT,DIM}
        CGB = [ConjGradBlock{TT}(G,GetOrder(DO))]
        new{TT,DIM,1,typeof(CGB[1].b)}(CGB)
    end
end

Base.getindex(CG::ConjGradMultiBlock,i::Integer) = CG.Block[i]

"""
    ExplicitBlock{T,N,O}
Data storage block for explicit integrators.
O<:Integer is the order of the integrator

Currently only Forward Euler and Runge-Kutta 4 are implemented for this solver.

TODO: Improve construction
"""
mutable struct ExplicitBlock{T,N,AT, O} <: DataBlockType{T,N, AT}
    k1  :: AT
    k2  :: AT
    k3  :: AT
    k4  :: AT
    Δt  :: T

    function ExplicitBlock{T}(grid::GridType,Δt::T,integrator::Symbol=:RK4) where T
        if typeof(grid) <: Grid1D
            n = grid.n
        elseif typeof(grid) <: Grid2D
            n = (grid.nx,grid.ny)
        end

        dims = length(n)

        k1 = zeros(T,n)
        if integrator == :RK4
            k2 = zeros(T,n)
            k3 = zeros(T,n)
            k4 = zeros(T,n)
            O = 4
        elseif integrator == :euler
            k2 = k3 = k4 = zeros(T,0)
            O = 1
        end

        new{T,dims,typeof(k1),O}(k1,k2,k3,k4,Δt)
    end
end

