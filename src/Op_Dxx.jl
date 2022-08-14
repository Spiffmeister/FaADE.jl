
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
function Dₓₓ(u::AbstractMatrix{Float64},nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,cx::AbstractMatrix{Float64},cy::AbstractMatrix{Float64};
    order_x::Int64=2,order_y::Int64=order_x)

    # Half order
    halfx = Int64(order_x/2) #half way
    halfy = Int64(order_y/2) #half way

    # Boundary node input selection
    order_x == 2 ? adjx = -1 : adjx = order_x
    order_y == 2 ? adjy = -1 : adjy = order_y

    # Needed for size of return array
    order_x == 2 ? retx = -1 : retx = halfx
    order_y == 2 ? rety = -1 : rety = halfy

    order_x == 2 ? intx = 1 : intx = order_x+1
    order_y == 2 ? inty = 1 : inty = order_y+1

    uₓₓ = SharedArray(zeros(Float64,nx,ny))

    # Internal nodes
    uₓₓ[intx+halfx:nx-(intx+halfx-1),inty+halfy:ny-(inty+halfy-1)] = Stencil2D(u[intx:nx-intx+1,inty:ny-inty+1],Internal,Internal,Δx,Δy,
            cx[intx:nx-intx+1,inty:ny-inty+1],
            cy[intx:nx-intx+1,inty:ny-inty+1],nx-2(order_x+retx),ny-2(order_y+rety),order_x=order_x,order_y=order_y)

    # Boundary nodes - avoiding corner nodes
    for j = order_y+retx+1:ny-(order_y+rety)
        # left and right x boundaries, for a matrix zeros(nx,ny) this is the 'top' and 'bottom' boundaries
        # x left
        uₓₓ[1:order_x+retx,j] = Stencil2D(u[1:order_x+adjx,j-halfy:j+halfy],
                Left,Internal,Δx,Δy,cx[1:order_x+adjx,j-halfy:j+halfy],cy[1:order_x+adjx,j-halfy:j+halfy],order_x=order_x,order_y=order_y)
        # x right
        uₓₓ[nx-(order_x+retx-1):nx,j] = Stencil2D(u[nx-order_x-adjx:nx,j-halfy:j+halfy],
                Right,Internal,Δx,Δy,cx[nx-order_x-adjx:nx,j-halfy:j+halfy],cy[nx-order_x-adjx:nx,j-halfy:j+halfy],order_x=order_x,order_y=order_y)
    end
    for i = order_x+retx+1:nx-(order_x+retx)
        # left and right y boundaries, for a matrix zeros(nx,ny) this is the 'left' and 'right' boundaries
        # y left
        uₓₓ[i,1:order_y+rety] = Stencil2D(u[i-halfx:i+halfx,1:order_y+adjy],
                Internal,Left ,Δx,Δy,cx[i-halfx:i+halfx,1:order_y+adjy],cy[i-halfx:i+halfx,1:order_y+adjy],order_x=order_x,order_y=order_y)
        # y right
        uₓₓ[i,ny-(order_y+rety-1):ny] = Stencil2D(u[i-halfx:i+halfx,ny-order_y-adjy:ny],
                Internal,Right,Δx,Δy,cx[i-halfx:i+halfx,ny-order_y-adjy:ny],cy[i-halfx:i+halfx,ny-order_y-adjy:ny],order_x=order_x,order_y=order_y)
    end

    # Corners - do entire corner as block
    # x left y left
    uₓₓ[1:order_x+retx,1:order_y+rety]          = 
        Stencil2D(u[1:order_x+adjx,1:order_y+adjy],Left,Left,Δx,Δy,
            cx[1:order_x+adjx,1:order_y+adjy],
            cy[1:order_x+adjx,1:order_y+adjy],order_x=order_x,order_y=order_y)
    # x left y right
    uₓₓ[1:order_x+retx,ny-order_y-rety+1:ny]      = 
        Stencil2D(u[1:order_x+adjx,ny-order_y-adjy:ny],Left,Right,Δx,Δy,
            cx[1:order_x+adjx,ny-order_y-adjy:ny],
            cy[1:order_x+adjx,ny-order_y-adjy:ny],order_x=order_x,order_y=order_y)
    # x right y left
    uₓₓ[nx-order_x-retx+1:nx,1:order_y+rety]      = 
        Stencil2D(u[nx-order_x-adjx:nx,1:order_y+adjy],Right,Left,Δx,Δy,
            cx[nx-order_x-adjx:nx,1:order_y+adjy],
            cy[nx-order_x-adjx:nx,1:order_y+adjy],order_x=order_x,order_y=order_y)
    # x right y right
    uₓₓ[nx-order_x-retx+1:nx,ny-order_y-rety+1:ny]  = 
        Stencil2D(u[nx-order_x-adjx:nx,ny-order_y-adjy:ny],Right,Right,Δx,Δy,
            cx[nx-order_x-adjx:nx,ny-order_y-adjy:ny],
            cy[nx-order_x-adjx:nx,ny-order_y-adjy:ny],order_x=order_x,order_y=order_y)

    return uₓₓ
end



"""
    Stencil2D

Computes the 2D finite difference stencil at a given node
"""
function Stencil2D end
### Iternal node stencil
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;order_x=2,order_y=order_x)

    halfx = Int64(order_x/2) #half way
    halfy = Int64(order_y/2) #half way

    uₓₓ = SharedArray(zeros(Float64,nx,ny))
    @sync @distributed for i = 1:nx
        for j = 1:ny
            uₓₓ[i,j] = SecondDerivative(u[i:i+order_x,j+halfx],
                    cx[i:i+order_x,j+halfx],Δx,Internal,order=order_x) + 
                SecondDerivative(u[i+halfy,j:j+order_y],
                    cy[i+halfy,j:j+order_y],Δy,Internal,order=order_y)
        end
    end

    return uₓₓ

end
### X boundaries
function Stencil2D(u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # side x internal y
    order_x == 2 ? adjx = 1 : adjx = order_x + Int64(order_x/2)
    adjy = Int64(order_y/2)+1 # half way
    
    uₓₓ = SecondDerivative(u[:,adjy],cy[:,adjy],Δx,xnode,order=order_x)
    for i = 1:adjx
        uₓₓ[i] += SecondDerivative(u[i,:],cx[i,:],Δy,Internal,order=order_y)
    end
    return uₓₓ

end
### Y boundaries
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # internal x side y
    adjx = Int64(order_x/2)+1 # half way
    order_y == 2 ? adjy = 1 : adjy = order_y + Int64(order_y/2)

    uₓₓ = SecondDerivative(u[adjx,:],cy[adjx,:],Δy,ynode,order=order_y)
    for j = 1:adjy
        uₓₓ[j] += SecondDerivative(u[:,j],cx[:,j],Δx,Internal,order=order_x)
    end
    return uₓₓ

end
### Corners
function Stencil2D(u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy;order_x=2,order_y=order_x)
    # corner x-y
    order_x == 2 ? adjx = -1 : adjx = Int64(order_x/2)
    order_y == 2 ? adjy = -1 : adjy = Int64(order_y/2)

    uₓₓ = zeros(order_x+adjx,order_y+adjy)

    for j = 1:order_y+adjy # x direction - top/bottom
        uₓₓ[:,j] += SecondDerivative(u[:,j],cx[:,j],Δx,xnode,order=order_x)
    end
    for i = 1:order_x+adjx # y direction - left/right
        uₓₓ[i,:] += SecondDerivative(u[i,:],cy[i,:],Δy,ynode,order=order_y)
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
