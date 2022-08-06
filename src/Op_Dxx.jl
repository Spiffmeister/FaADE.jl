
"""
Dₓₓ(u::Vector{Float64},c::Vector{Float64},n::Int64,Δx::Float64;order::Int64=2)
or
Dₓₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::Matrix{Float64};dim::Int64=1,order::Int64=2)

- If `typeof(u) <: Vector{Float64}` then uses a 1D second derviative SBP operator
- Returns a `Vector{Float64}` of length `n`.
- If `typeof(u) <: Matrix{Float64}` then uses a 2D second derivative SBP operator
- dim ∈ [1,2]
    - If `dim==1` then takes derivative along rows (`u[:,i]`)
    - If `dim==2` then takes derivative along columns (`u[i,:]`)
- Δ should be the grid spacing along the dimension requested
Returns a `Matrix{Float64}` of size `nx × ny`.

This is also available is an iterator [`Dₓₓ!`](@ref)

Internally uses [`SecondDerivative`](@ref)
"""
function Dₓₓ end
### 1D second derivative operator
function Dₓₓ(u::AbstractVector{Float64},c::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    # Function call for the 1D 2nd derivative SBP operator

    uₓₓ = zeros(Float64,n)

    Dₓₓ!(uₓₓ,u,c,n,Δx,order=order)

    return uₓₓ
end
### 2D second derivative operator
function Dₓₓ(u::AbstractMatrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix{Float64};dim::Int64=1,order::Int64=2)
# Multidimensional call for 2nd derivative SBP operator

    uₓₓ = zeros(Float64,nx,ny)
    if dim == 1
        for i = 1:ny #column derivatives
            uₓₓ[:,i] = Dₓₓ!(uₓₓ[:,i],u[:,i],c[:,i],nx,Δ,order=order)
        end
    elseif dim == 2
        for i = 1:nx #row derivative
            uₓₓ[i,:] = Dₓₓ!(uₓₓ[i,:],u[i,:],c[i,:],ny,Δ,order=order)
        end
    else
        error("dim must be 1 or 2.")
    end

    return uₓₓ
end



"""
Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
or
Dₓₓ!(uₓₓ::AbstractMatrix{Float64},u::AbstractMatrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix{Float64};dim::Int64=1,order::Int64=2)
or
Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},m::Int64,Δ::Float64;order::Int64=2)

Iterator for 1D and 2D second derviative SBP operator
"""
function Dₓₓ! end
### Vector based second derviative iterator
function Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)

    adj = Int64(order/2)

    for i = order:n-order+1
        uₓₓ[i] = SecondDerivative(u[i-adj:i+adj],c[i-adj:i+adj],Δx,Internal,order=order)
    end

    adj += order

    if order == 2
        uₓₓ[1] = 0.0
        uₓₓ[n] = 0.0
    else
        uₓₓ[1:adj] = SecondDerivative(u[1:2order],c[1:2order],Δx,Left,order=order)
        uₓₓ[n-adj+1:n] = SecondDerivative(u[n-2order+1:n],c[n-2order+1:n],Δx,Right,order=order)
    end

    return uₓₓ
end
### Multidimensional second derivative SBP operator
function Dₓₓ!(uₓₓ::AbstractMatrix,u::AbstractMatrix,nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix;dim::Int64=1,order::Int64=2)

    if dim == 1
        @sync @distributed for i = 1:ny #column derivatives
            uₓₓ[:,i] = Dₓₓ!(uₓₓ[:,i],u[:,i],c[:,i],nx,Δ,order=order)
        end
    elseif dim == 2
        @sync @distributed for i = 1:nx #row derivative
            uₓₓ[i,:] = Dₓₓ!(uₓₓ[i,:],u[i,:],c[i,:],ny,Δ,order=order)
        end
    else
        error("dim must be 1 or 2.")
    end

    return uₓₓ
end



### 2D at nodes

# function Stencil2D(uₓₓ::AbstractMatrix,u::AbstractMatrix,nx,ny,Δx,Δy,cx,cy;order_x=2,order_y=order_x)

#     adjx = Int64(order/2)
#     adjy = Int64(order/2)

#     # Internal nodes
#     @sync @distributed for i = order_x:nx-order_x+1
#         @distributed for j = order_y:ny-order_y+1
#             uₓₓ[i] = SecondDerivative(u[i-adjx:i+adjx,j],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x) + 
#                 SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Internal,order=order_y)
#         end
#     end

# end


"""
    Stencil2D

Computes the 2D finite difference stencil at a given node
"""
function Stencil2D end

### Iternal node stencil
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;order_x=2,order_y=order_x)
    adjx = Int64(order/2)+1 #half way
    adjy = Int64(order/2)+1 #half way

    i = Int64(order/2)+1 #half way
    j = Int64(order/2)+1 #half way

    uₓₓ = zeros(Float64,nx,ny)

    for i = adjx:nx-adjx
        for j = adjy:ny-adjy
            uₓₓ[i] = SecondDerivative(u[:,j],cx[:,j],Δx,Internal,order=order_x) + 
                SecondDerivative(u[i,:],cy[i,:],Δy,Internal,order=order_y)
        end
    end

    return uₓₓ

end
### X boundaries
function Stencil2D(u::AbstractMatrix,::NodeType{:Left},::NodeType{:Internal},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # Left x internal y
    i = Int64(order/2)+1
    j = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[:,j],cx[:,j],Δx,Left,order=order_x)
    for i = 1:order_x
        uₓₓ += SecondDerivative(u[i,:],cy[i,:],Δy,Internal,order=order_y)
    end
    return uₓₓ

end
function Stencil2D(u::AbstractMatrix,::NodeType{:Right},::NodeType{:Internal},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # Right x internal y
    i = Int64(order/2)+1
    j = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[:,j],cx[:,j],Δx,Right,order=order_x) + 
    for i = 1:order_x
        uₓₓ += SecondDerivative(u[i,:],cy[i,:],Δy,Internal,order=order_y)
    end
    return uₓₓ

end
### Y boundaries
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Left},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # internal x left y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Left,order=order_y)
    for j = j:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Right},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # internal x right y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Right,order=order_y)
    for j = 1:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
### Corners
function Stencil2D(u::AbstractMatrix,::NodeType{:Left},::NodeType{Left},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # left x left y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Right,order=order_y)
    for j = 1:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
function Stencil2D(u::AbstractMatrix,::NodeType{:Left},::NodeType{:Right},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # left x right y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Right,order=order_y)
    for j = 1:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
function Stencil2D(u::AbstractMatrix,::NodeType{:Right},::NodeType{:Left},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # right x left y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Right,order=order_y)
    for j = 1:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
function Stencil2D(u::AbstractMatrix,::NodeType{:Right},::NodeType{:Right},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # right x right y
    adjx = Int64(order/2)+1
    adjy = Int64(order/2)+1

    uₓₓ = SecondDerivative(u[i,j-adjy:j+adjy],cy[i,j-adjy:j+adjy],Δy,Right,order=order_y)
    for j = 1:order_y
        uₓₓ += SecondDerivative(u[:,adjx],cx[i-adjx:i+adjx,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end

### Chunked Arrays
# function Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},m::Int64,Δ::Float64;order::Int64=2)
#     adj = Int64(order/2)
#     for i = 1+adj:m-adj # Avoid the "ghost nodes" by doing adj:m-adj
#         uₓₓ[i] = SecondDerivative(u[i-adj:i+adj],c[i-adj:i+adj],Δ,NodeInternal,order=order)
#     end
#     return uₓₓ
# end
