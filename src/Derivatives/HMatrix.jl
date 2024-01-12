

abstract type MassMatrix{TYPE,TT,VT} end

struct DiagonalH{TT<:Real,VT<:AbstractVector{TT}} <: MassMatrix{:Diagonal,TT,VT}
    Boundary    :: VT
    Interior    :: TT
    Δ           :: TT
    n           :: Int64
    Bwidth      :: Int64

    function DiagonalH(order::Int,Δx::TT,n::Int) where {TT}
        if order == 2
            H = TT.([1//2])
        elseif order == 4
            H = TT.([17//48, 59//48, 43//48, 49//48])
        end

        interior = one(TT)

        new{TT,Vector{TT}}(H,interior,Δx,n,length(H))
    end
end


struct CompositeH{DIM,TT<:Real,VT<:AbstractArray,HTYPE}
    H :: NTuple{DIM,HTYPE}
    sz :: Vector{Int}
    function CompositeH(H...)
        sz = [H[i].n for i in 1:length(H)]
        new{length(H),eltype(H[1].Boundary),typeof(H[1].Boundary),typeof(H[1])}(H,sz)
    end
end




function mul!(u::VT,H::DiagonalH{TT},v::VT) where {TT,VT<:AbstractVector{TT}}
    tmp = TT(0)
    for i = 1:B.n
        tmp += u[i] * H[i] * v[i]
    end
    return tmp
end
function mul!(u::AT,H::CompositeH{DIM,TT},v::AT) where {DIM,TT,AT<:AbstractArray{TT}}
    tmp = TT(0)
    for j = 1:H.sz[2]
        for i = 1:H.sz[1]
            tmp += u[i,j] * H[1][i] * H[2][j] * v[i,j]
        end
    end
    return tmp
end



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
function innerH(n)
    H = build_H(2,n)
    H[1] = 1.0
    H[end] = 1.0
    return innerH{Float64,1,typeof(H)}(H,[1.0], 1.0)
end
function innerH(Δx::T,n,order) where T
    H = build_H(order,n)
    return innerH{T,1,typeof(H)}(H,[1.0], Δx)
end
function innerH(Δx::T,Δy::T,nx,ny,order) where T
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
            @inbounds tmp += u[j,i] * H.Hx[j] * H.Hy[i] * v[j,i] * H.Δ
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





Base.size(H::DiagonalH) = (H.n,H.n)
Base.size(H::CompositeH) = (H.sz[1],H.sz[2])

Base.length(H::DiagonalH) = H.n
Base.length(H::CompositeH) = prod(H.sz)

Base.lastindex(H::CompositeH) = length(H)

function getH(H::CompositeH{DIM,TT,VT,HTYPE},i::Int) where {DIM,TT,VT<:AbstractVector{TT},HTYPE}
    
end

function Base.getindex(H::DiagonalH,i::Int)
    if i < 1 || i > H.n
        error("Index out of bounds")
    elseif i <= H.Bwidth
        return H.Boundary[i] * H.Δ
    elseif i > H.n - H.Bwidth
        return H.Boundary[H.n-i+1] * H.Δ
    else
        return H.Interior * H.Δ
    end
end
function Base.getindex(H::DiagonalH{TT},i::Int,j::Int) where {TT}
    if i != j
        return TT(0)
    else
        return H[i]*H.Δ
    end
end
function Base.getindex(H::CompositeH{2,TT,VT,HTYPE},i::Int,j::Int) where {TT,VT<:Vector{TT},HTYPE}
    Hx = H.H[1] :: HTYPE
    Hy = H.H[2] :: HTYPE
    return Hx[i]*Hy[j]
end

