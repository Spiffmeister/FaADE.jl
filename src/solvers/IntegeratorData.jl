


"""
    ConjGradBlock
Data storage for the conjugate gradient method. See [`conj_grad!`](@ref)
"""
mutable struct ConjGradBlock{T,N, AT} <: DataBlockType{T,N,AT}
    b   :: AT # b = uⁿ⁺¹ + F
    rₖ  :: AT # (uⁿ⁺¹ - Δt*uₓₓⁿ⁺¹) - b
    Adₖ :: AT # Adₖ = dₖ - Δt*D(dₖ)
    Drₖ :: AT # Drₖ = rₖ - Δt*D(rₖ)
    dₖ  :: AT # dₖ = -rₖ, -rₖ .+ βₖ*dₖ

    converged   :: Bool

    innerprod   :: innerH{T,N,Vector{T}}

    function ConjGradBlock{T}(grid::GridType,order::Int) where T

        if typeof(grid) <: Grid1D
            Δ = grid.Δx
            n = grid.n
        elseif typeof(grid) <: Grid2D
            Δ = min(grid.Δx,grid.Δy)
            n = (grid.nx,grid.ny)
        end

        b   = zeros(T,n)
        rₖ  = zeros(T,n)
        Adₖ = zeros(T,n)
        Arₖ = zeros(T,n)
        dₖ  = zeros(T,n)

        dims = length(n)

        innerprod = innerH(grid,order)
        

        new{T,dims, typeof(b)}(b, rₖ, Adₖ, Arₖ, dₖ, true, innerprod)
    end
end




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

