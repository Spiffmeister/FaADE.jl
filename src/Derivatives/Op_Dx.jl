#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru

"""
First derivative SBP operator. 
    
Also available as and internally uses in place operator [`D₁!`](@ref).
"""
function D₁ end
"""
    D₁(u::AbstractVector{T},Δx::T;order::Integer=2) where T
1D implementation of `D₁` operator.

```
julia> n = 101
julia> x = collect(LinRange(0.0,1.0,n))

julia> u = sin.(x)
julia> Δx = 1/(n-1)

julia> D₁(u,Δx,order=2)
```
"""
function D₁(u::AbstractVector{T},Δx::T;order::Integer=2) where T
    uₓ = zeros(T,length(u))
    DO = Val(order)
    D₁!(uₓ,u,length(u),Δx,DO,0.0)
    return uₓ
end
"""
    D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;order::Integer=2) where T
2D implementation of `D₁` operator.
"""
function D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;order::Integer=2) where T
    uₓ = zeros(T,size(u))
    D₁!(uₓ,u,nx,ny,Δx,Δy,order,order,0.0)
    return uₓ
end


"""
1D and 2D in place first derivative operator.

See also [`FirstDerivativeBoundary!`](@ref) and [`FirstDerivativeInternal!`](@ref).
"""
function D₁! end
"""
    D₁!(uₓ::AbstractVector{T},u::AbstractVector{T},n::Integer,Δx::T,order::Integer)
1D [`D₁!`](@ref).
"""
function D₁!(uₓ::AT,u::AT,n::Integer,Δx::T,order::Val,α::T) where {T,AT<:AbstractVector{T}}
    FirstDerivativeBoundary!(uₓ,u,Δx,Left,order,α)
    FirstDerivativeInternal!(uₓ,u,Δx,n,order,α)
    FirstDerivativeBoundary!(uₓ,u,Δx,Right,order,α)
end
function D₁!(uₓ::AT,c::AT,u::AT,n::Integer,Δx::T,order::Int,α::T) where {T,AT<:AbstractVector{T}}
    O = Val(order)
    FirstDerivativeBoundary!(uₓ,c,u,Δx,Left,O,α)
    FirstDerivativeInternal!(uₓ,c,u,Δx,n,O,α)
    FirstDerivativeBoundary!(uₓ,c,u,Δx,Right,O,α)
end
"""
    D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,order::Integer,α::T,dim::Integer) where T
1D implementation for 2D problems for [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,order::Integer,α::T,dim::Integer) where T
    loopdir = _SelectLoopDirection(dim)
    ORD = Val(order)
    foreach(zip(eachslice(uₓ,dims=dim),eachslice(u,dims=dim))) do (cache,U)
    # for (cache,U) in zip(loopdir(uₓ),loopdir(u))
    # foreach(zip(loopdir(uₓ),loopdir(u))) do (cache,U)
        D₁!(cache,U,n,Δ,ORD,α)
    end
    # uₓ
end
"""
    D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,order::Int,ordery::Int,α::T) where T
2D [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,order::Int,ordery::Int,α::T) where T
    
    order == 2 ? ret = 1 : ret = order
    order == 2 ? nin = 2 : nin = order + halforder(order)
    order == 2 ? ein = 1 : ein = halforder(order)

    order_x = Val(order)
    order_y = Val(ordery)
    
    for (A,B) in zip(eachcol(uₓ),eachcol(u))
        D₁!(A,B,nx,Δx,order_x,α)
    end
    for (A,B) in zip(eachrow(uₓ),eachrow(u))
        D₁!(A,B,ny,Δy,order_y,1.0)
    end

    uₓ
end

"""
    D₁ᵀ!(dest::VT,u::VT,n::Int,Δx::TT,order::Int,α::TT) where {TT,VT<:AbstractVector{TT}}
1D and 2D in place first derivative transpose operator.
"""
function D₁ᵀ!(dest::VT,u::VT,n::Int,Δx::TT,order::Int,α::TT) where {TT,VT<:AbstractVector{TT}}
    order == 2 ? m = 2 : m = 7

    for i = m:n-m+1
        @inbounds dest[i] = α*dest[i] + FirstDerivativeInternal(u,Δx,Val(order),i,TT(1))
    end
    
    FirstDerivativeBoundaryTranspose!(dest,u,Δx,Left,   Val(order),α)
    FirstDerivativeBoundaryTranspose!(dest,u,Δx,Right,  Val(order),α)
end
