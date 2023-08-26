
abstract type IntegratorDataMultiBlock{dtype,DIM} <: newDataBlockType{dtype,DIM} end
abstract type IntegratorDataBlock{dtype,DIM} <: newDataBlockType{dtype,DIM} end
abstract type IntegratorConfig end

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
struct ConjGradBlock{TT,DIM,AT} <: IntegratorDataBlock{TT,DIM}
    b       :: AT # b = uⁿ⁺¹ + F
    rₖ      :: AT # (uⁿ⁺¹ - Δt*uₓₓⁿ⁺¹) - b
    Adₖ     :: AT # Adₖ = dₖ - Δt*D(dₖ)
    Drₖ     :: AT # Drₖ = rₖ - Δt*D(rₖ)
    dₖ      :: AT # dₖ = -rₖ, -rₖ .+ βₖ*dₖ
    cache   :: AT   # Adₖ and Drₖ

    # buffer_left     :: AT
    # buffer_right    :: AT
    # buffer_up       :: AT
    # buffer_down     :: AT

    scalar :: ConjGradScalarBlock{TT}

    innerprod   :: innerH{TT,DIM,Vector{TT}}

    function ConjGradBlock{TT}(grid::GridType{TT,DIM},order::Union{Int,DerivativeOrder}) where {TT,DIM}

        n = size(grid)

        b   = zeros(TT,n)
        rₖ  = zeros(TT,n)
        Adₖ = zeros(TT,n)
        Arₖ = zeros(TT,n)
        dₖ  = zeros(TT,n)
        cache = zeros(TT,n)

        nnodes = SATNodeOutput(order)
        if DIM == 1
            buffer_left = zeros(TT,nnodes)
            buffer_up   = zeros(TT,1)
        elseif DIM == 2
            buffer_left = zeros(TT,(nnodes,grid.ny))
            buffer_up   = zeros(TT,(grid.nx,nnodes))
        end
        buffer_right    = similar(buffer_left)
        buffer_down     = similar(buffer_up)

        innerprod = innerH(grid,order)

        scalars = ConjGradScalarBlock{TT}()

        new{TT,DIM, typeof(b)}(b, rₖ, Adₖ, Arₖ, dₖ, cache, scalars, innerprod)
    end
end


mutable struct ConjugateGradientConfig <: IntegratorConfig
    converged   :: Bool
end

"""
    ConjGradMultiBlock
"""
struct ConjGradMultiBlock{TT<:Real,
        DIM,
        TID <: IntegratorConfig} <: IntegratorDataBlock{TT,DIM}

    Block   :: Vector{ConjGradBlock}
    nblock  :: Int64    # Number of data blocks
    SC      :: TID

    function ConjGradMultiBlock(G::GridMultiBlock{TT,DIM},DO::DerivativeOrder) where {TT,DIM}
        CGB = [ConjGradBlock{TT}(grid,GetOrder(DO)) for grid in G.Grids]

        new{TT,DIM,ConjugateGradientConfig}(CGB,length(CGB),ConjugateGradientConfig(true))
    end
    function ConjGradMultiBlock(G::LocalGridType{TT,DIM},DO::DerivativeOrder) where {TT,DIM}
        CGB = [ConjGradBlock{TT}(G,GetOrder(DO))]
        new{TT,DIM,ConjugateGradientConfig}(CGB,length(CGB),ConjugateGradientConfig(true))
    end
end

Base.getindex(CG::ConjGradMultiBlock,i::Integer) = CG.Block[i]
eachblock(CG::ConjGradMultiBlock) = Base.OneTo(CG.nblock)

# CG[i].Block.scalar.converged

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

