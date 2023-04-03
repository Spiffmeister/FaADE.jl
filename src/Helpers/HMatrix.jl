



"""
    innerH
Method for a 1 or 2 dimensional H inner product.
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
    build_H
Build the H matrix in the SBP operator Dₓ = H⁻¹Q
"""
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



