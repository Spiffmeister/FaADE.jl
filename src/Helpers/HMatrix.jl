



"""
    innerH{T,N,VT<:AbstractVector{T}}
Method for a 1 or 2 dimensional H inner product.
        
Inputs:
- [`GridType`](@ref)
- Order of method

Returns:
- Method-like struct for H-inner products
"""
struct innerH{T<:Real,
        N,
        VT<:AbstractVector{T}}
    Hx  :: VT
    Hy  :: VT
    Δ   :: T
end
function innerH(Δx::T,n,order::Int) where {T}
    H = build_H(order,n)
    new{T,1,typeof(H)}(H,[1.0], Δx)
end
function innerH(Δx::T,Δy::T,nx,ny,order::Int) where {T} 
    Hx = build_H(order,nx)
    Hy = build_H(order,ny)
    return innerH{T,2,typeof(Hx)}(Hx,Hy, Δx*Δy)
end

"""
    (H::InnerH)
1D or 2D H-inner product constructed from [`innerH`](@ref)
"""
function (H::innerH{T,1})(u::AbstractArray{T},v::AbstractArray{T}) :: T where T
    local tmp = T(0)
    for i in eachindex(H.Hx)
        tmp += u[i] * H.Hx[i] * v[i] * H.Δ
    end
    return tmp
end
function (H::innerH{T,2})(u::AbstractArray{T,2},v::AbstractArray{T,2}) :: T where T
    local tmp::T
    tmp = T(0)
    for j in eachindex(H.Hx)
        for i in eachindex(H.Hy)
            @inbounds tmp += u[i,j] * H.Hx[i] * H.Hy[j] * v[i,j] * H.Δ
        end
    end
    return tmp
end



"""
    build_H
Build the H matrix in the SBP operator D₁ = H⁻¹Q
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



