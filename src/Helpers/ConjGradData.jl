



"""
innerH
Method for a 1 or 2 dimensional H inner product
"""
struct innerH{T,N}
    Hx  :: AbstractArray{T}
    Hy  :: AbstractArray{T}
    nx  :: Int
    ny  :: Int
    Δ   :: T

    function innerH(grid::GridType{T,N},order::Int) where {T,N}
        if typeof(grid) <: Grid1D
            H = build_H(order,grid.n)
            new{T,N}(H,[1.0], grid.n,1, grid.Δx)
            
        elseif typeof(grid) <: Grid2D
            Hx = build_H(order,grid.nx)
            Hy = build_H(order,grid.ny)
            new{T,N}(Hx,Hy, grid.nx,grid.ny, grid.Δx*grid.Δy)
        end
    end
end

function build_H(order,n)
    H = ones(Float64,n)
    if order == 2
        H[1] = H[end] = 0.5
    elseif order == 4
        H[1] = H[end]   = 17.0/48.0
        H[2] = H[end-1] = 59.0/48.0
        H[3] = H[end-2] = 43.0/48.0
        H[4] = H[end-3] = 49.0/48.0
    else
        error("Order must be 2 or 4")
    end
    return H
end



function (H::innerH{T,1})(u::AbstractArray{T},v::AbstractArray{T}) where T
    tmp = 0.0
    for i = 1:H.nx
        tmp += u[i] * H.Hx[i] * v[i] * H.Δ
    end
    return tmp
end

function (H::innerH{T,2})(u::AbstractArray{T},v::AbstractArray{T}) where T
    tmp = 0.0
    for j = 1:H.ny
        for i = 1:H.nx
            tmp += u[i,j] * H.Hx[i] * H.Hy[j] * v[i,j] * H.Δ
        end
    end
    return tmp
end



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



mutable struct ExplicitBlock{T,N} <: DataBlockType{T,N}
    k1  :: AbstractArray{T,N}
    k2  :: AbstractArray{T,N}
    k3  :: AbstractArray{T,N}
    k4  :: AbstractArray{T,N}

    function ExplicitBlock{T}(grid::GridType,order::Int;ntegrator::Symbol=:RK4) where T
        if typeof(grid) <: Grid1D
            Δ = mind(grid.Δx,grid.Δy)
            n = grid.n
        elseif typeof(grid) <: Grid2D
            n = (grid.nx,grid.ny)
        end

        dims = length(n)

        if integrator == :RK4
            k1 = zeros(T,n)
            k2 = zeros(T,n)
            k3 = zeros(T,n)
            k4 = zeros(T,n)
        end

        new{T,dims}(k1,k2,k3,k4)
    end
end

