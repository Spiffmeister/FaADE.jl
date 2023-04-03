


"""
    ConjGradBlock
"""
mutable struct ConjGradBlock{T,N} <: DataBlockType{T,N}
    b   :: AbstractArray{T,N} # b = uⁿ⁺¹ + F
    rₖ  :: AbstractArray{T,N} # (uⁿ⁺¹ - Δt*uₓₓⁿ⁺¹) - b
    Adₖ :: AbstractArray{T,N} # Adₖ = dₖ - Δt*D(dₖ)
    Drₖ :: AbstractArray{T,N} # Drₖ = rₖ - Δt*D(rₖ)
    dₖ  :: AbstractArray{T,N} # dₖ = -rₖ, -rₖ .+ βₖ*dₖ

    converged   :: Bool

    innerprod   :: innerH{T,N}

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
        

        new{T,dims}(b, rₖ, Adₖ, Arₖ, dₖ, true, innerprod)
    end
end




"""
    ExplicitBlock{T,N,O}
Data storage block for explicit integrators.
O<:Integer is the order of the integrator
"""
mutable struct ExplicitBlock{T,N,O} <: DataBlockType{T,N}
    k1  :: AbstractArray{T,N}
    k2  :: AbstractArray{T,N}
    k3  :: AbstractArray{T,N}
    k4  :: AbstractArray{T,N}
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

        new{T,dims,O}(k1,k2,k3,k4,Δt)
    end
end

