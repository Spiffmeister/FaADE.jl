#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru

"""
    D₁
1D and 2D first derivative operator. 
    
Also available as and internally uses in place operator [`D₁!`](@ref).
"""
function D₁ end
"""
    D₁function D₁(u::AbstractVector{T},Δx::T;order::Integer=2)
1D implementation of ``D_x`` operator.
"""
function D₁(u::AbstractVector{T},Δx::T;
        order::Integer=2) where T
    uₓ = zeros(T,length(u))
    DO = Val(order)
    D₁!(uₓ,u,length(u),Δx,DO,0.0)
    return uₓ
end
"""
    D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;order::Integer=2)
2D implementation of ``D_x`` operator
"""
function D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;
        order::Integer=2) where T
    uₓ = zeros(T,size(u))
    D₁!(uₓ,u,nx,ny,Δx,Δy,order,order,0.0)
    return uₓ
end


"""
    D₁!
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
function D₁!(uₓ::AT,c::AT,u::AT,n::Integer,Δx::T,order::Val,α::T) where {T,AT<:AbstractVector{T}}
    FirstDerivativeBoundary!(uₓ,c,u,Δx,Left,order,α)
    FirstDerivativeInternal!(uₓ,c,u,Δx,n,order,α)
    FirstDerivativeBoundary!(uₓ,c,u,Δx,Right,order,α)
end
"""
    function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,order::Integer,dim::Integer)
1D implementation for 2D problems for [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,order::Integer,α::T,dim::Integer) where T
    loopdir = SelectLoopDirection(dim)
    for (cache,U) in zip(loopdir(uₓ),loopdir(u))
        D₁!(cache,U,n,Δ,order,α)
    end
    uₓ
end
"""
    function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,order::Integer)
2D [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,
        order::Int,ordery::Int,α::T) where T
    
    order == 2 ? ret = 1 : ret = order
    order == 2 ? nin = 2 : nin = order + halforder(order)
    order == 2 ? ein = 1 : ein = halforder(order)

    order_x = Val(order)
    order_y = Val(ordery)
    
    for (A,B,C) in zip(eachcol(uₓ),eachcol(u),eachcol(cx))
        D₁!(A,B,C,nx,Δx,order_x,α)
    end
    for (A,B,C) in zip(eachrow(uₓ),eachrow(u),eachrow(cy))
        D₁!(A,B,C,ny,Δy,order_y,1.0)
    end

    
    uₓ
end

"""
    D₁ᵀ!
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
