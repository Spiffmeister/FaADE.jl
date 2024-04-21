#======================================#
#====== Second Derivative Methods =====#
#======================================#
# Author: Dean Muir, Kenneth Duru



"""
    D₂
Not in place second derivative operators
"""
function D₂ end
"""
Vector version of [`D₂`](@ref).
"""
function D₂(u::VT,c::VT,Δx::TT,order::Int=2) where {TT,VT<:AbstractVector{TT}}
    uₓₓ = zeros(eltype(u),length(u))
    DO = Val(order)
    D₂!(uₓₓ,u,c,length(u),Δx,DO,0.0)
    return uₓₓ
end
"""
Matrix version of [`D₂`](@ref).
"""
function D₂(u::AT,cx::AT,cy::AT,Δx::TT,Δy::TT,order::Int=2) where {TT,AT <: AbstractMatrix{TT}}
    n = size(u)
    uₓₓ = zeros(eltype(u),size(u))
    DO = Val(order)
    D₂!(uₓₓ,u,cx,cy,n[1],n[2],Δx,Δy,DO,DO,TT(0))
    return uₓₓ
end


"""
    D₂!
1D and 2D in place second derivative operator.

Internally uses [`SecondDerivativeInternal`](@ref) and [`SecondDerivativeInternal`](@ref).
"""
function D₂! end
### 1D second derviative function
function D₂!(uₓₓ::AbstractVector{T},u::AbstractVector{T},c::AbstractVector{T},n::Integer,Δx::T,order::Val,α::T) where T
    SecondDerivativeInternal!(uₓₓ,u,c,Δx,n,order,α)
    SecondDerivativeBoundary!(uₓₓ,u,c,Δx,Left,order,α)
    SecondDerivativeBoundary!(uₓₓ,u,c,Δx,Right,order,α)
    uₓₓ
end
### Multidimensional second derivative SBP operator, select the axis to differentiate across by dim
function D₂!(uₓₓ::AT,u::AT,c::AT,n::Integer,Δ::T,order::Val,α::T,dim::Integer) where {T,AT<:AbstractArray{T}}
    loopdir = SelectLoopDirection(dim)    
    for (cache,U,C) in zip(loopdir(uₓₓ),loopdir(u),loopdir(c))
        D₂!(cache,U,C,n,Δ,order,α)
    end
    uₓₓ
end
### 2D second derviative function
function D₂!(uₓₓ::AT,u::AT,cx::AT,cy::AT,
        nx::Integer,ny::Integer,Δx::T,Δy::T,
        order_x::Val,order_y::Val,α::T) where {T,AT<:AbstractMatrix{T}}
    
    for (A,B,C) in zip(eachcol(uₓₓ),eachcol(u),eachcol(cx))
        D₂!(A,B,C,nx,Δx,order_x,α)
    end
    for (A,B,C) in zip(eachrow(uₓₓ),eachrow(u),eachrow(cy))
        D₂!(A,B,C,ny,Δy,order_y,T(1))
    end

    uₓₓ
end


