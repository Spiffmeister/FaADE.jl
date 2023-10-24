#======================================#
#====== Second Derivative Methods =====#
#======================================#
# Author: Dean Muir, Kenneth Duru



"""
    D₂
Not in place second derivative operators
"""
function D₂ end
function D₂(u::TA,c::TA,n::Integer,Δx::TT,order::DerivativeOrder) where {TT,TA <: AbstractVector{TT}}
    uₓₓ = zeros(eltype(u),length(u))
    D₂!(uₓₓ,u,c,n,Δx,order,0.0)
    return uₓₓ
end
function D₂(u::TA,cx::TA,cy::TA,nx::Integer,ny::Integer,Δx::TT,Δy::TT,
    orderx::DO,ordery::DO) where {TT,TA <: AbstractVector{TT},DO<:DerivativeOrder}
    uₓₓ = zeros(eltype(u),size(u))
    D₂!(uₓₓ,u,cx,cy,nx,ny,Δx,Δy,orderx,ordery)
    return uₓₓ
end


"""
    D₂!
1D and 2D in place second derivative operator.

Internally uses [`SecondDerivativeInternal`](@ref) and [`SecondDerivativeInternal`](@ref).
"""
function D₂! end
### 1D second derviative function
function D₂!(uₓₓ::AbstractVector{T},u::AbstractVector{T},c::AbstractVector{T},n::Integer,Δx::T,order::DerivativeOrder,α::T) where T
    SecondDerivativeInternal!(uₓₓ,u,c,Δx,n,order,α)
    SecondDerivativeBoundary!(uₓₓ,u,c,Δx,Left,order,α)
    SecondDerivativeBoundary!(uₓₓ,u,c,Δx,Right,order,α)
    uₓₓ
end
### Multidimensional second derivative SBP operator, select the axis to differentiate across by dim
function D₂!(uₓₓ::AT,u::AT,c::AT,n::Integer,Δ::T,order::DerivativeOrder,α::T,dim::Integer) where {T,AT<:AbstractArray{T}}
    loopdir = SelectLoopDirection(dim)    
    for (cache,U,C) in zip(loopdir(uₓₓ),loopdir(u),loopdir(c))
        D₂!(cache,U,C,n,Δ,order,α)
    end
    uₓₓ
end
### 2D second derviative function
function D₂!(uₓₓ::AT,u::AT,cx::AT,cy::AT,
        nx::Integer,ny::Integer,Δx::T,Δy::T,
        # order_x::DerivativeOrder,order_y::DerivativeOrder) where {T,AT}
        order_x::DerivativeOrder,order_y::DerivativeOrder,α::T) where {T,AT<:AbstractMatrix{T}}
    
    for (A,B,C) in zip(eachcol(uₓₓ),eachcol(u),eachcol(cx))
        D₂!(A,B,C,nx,Δx,order_x,α)
    end
    for (A,B,C) in zip(eachrow(uₓₓ),eachrow(u),eachrow(cy))
        D₂!(A,B,C,ny,Δy,order_y,T(1))
    end
    

    uₓₓ
end

