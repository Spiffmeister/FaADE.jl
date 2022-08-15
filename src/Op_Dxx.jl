
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
    ordx::Int64=2,ordy::Int64=ordx)

    # Half order
    halfx = Int64(ordx/2) #half way
    halfy = Int64(ordy/2) #half way

    # Boundary node input selection
    ordx == 2 ? adjx = -1 : adjx = ordx
    ordy == 2 ? adjy = -1 : adjy = ordy

    # Needed for size of return array
    ordx == 2 ? retx = -1 : retx = halfx
    ordy == 2 ? rety = -1 : rety = halfy

    ordx == 2 ? intx = 1 : intx = ordx+1
    ordy == 2 ? inty = 1 : inty = ordy+1

    uₓₓ = SharedArray(zeros(Float64,nx,ny))

    # Internal nodes
    uₓₓ[intx+halfx:nx-(intx+halfx-1),inty+halfy:ny-(inty+halfy-1)] = Stencil2D!(uₓₓ[intx+halfx:nx-(intx+halfx-1),inty+halfy:ny-(inty+halfy-1)],u[intx:nx-intx+1,inty:ny-inty+1],Internal,Internal,Δx,Δy,
            cx[intx:nx-intx+1,inty:ny-inty+1],
            cy[intx:nx-intx+1,inty:ny-inty+1],nx-2(ordx+retx),ny-2(ordy+rety),ordx=ordx,ordy=ordy)

    # Boundary nodes - avoiding corner nodes
    for j = ordy+retx+1:ny-(ordy+rety)
        # left and right x boundaries, for a matrix zeros(nx,ny) this is the 'top' and 'bottom' boundaries
        # x left
        uₓₓ[1:ordx+retx,j] = Stencil2D(u[1:ordx+adjx,j-halfy:j+halfy],
                Left,Internal,Δx,Δy,cx[1:ordx+adjx,j-halfy:j+halfy],cy[1:ordx+adjx,j-halfy:j+halfy],ordx=ordx,ordy=ordy)
        # x right
        uₓₓ[nx-(ordx+retx-1):nx,j] = Stencil2D(u[nx-ordx-adjx:nx,j-halfy:j+halfy],
                Right,Internal,Δx,Δy,cx[nx-ordx-adjx:nx,j-halfy:j+halfy],cy[nx-ordx-adjx:nx,j-halfy:j+halfy],ordx=ordx,ordy=ordy)
    end
    for i = ordx+retx+1:nx-(ordx+retx)
        # left and right y boundaries, for a matrix zeros(nx,ny) this is the 'left' and 'right' boundaries
        # y left
        uₓₓ[i,1:ordy+rety] = Stencil2D(u[i-halfx:i+halfx,1:ordy+adjy],
                Internal,Left ,Δx,Δy,cx[i-halfx:i+halfx,1:ordy+adjy],cy[i-halfx:i+halfx,1:ordy+adjy],ordx=ordx,ordy=ordy)
        # y right
        uₓₓ[i,ny-(ordy+rety-1):ny] = Stencil2D(u[i-halfx:i+halfx,ny-ordy-adjy:ny],
                Internal,Right,Δx,Δy,cx[i-halfx:i+halfx,ny-ordy-adjy:ny],cy[i-halfx:i+halfx,ny-ordy-adjy:ny],ordx=ordx,ordy=ordy)
    end

    # Corners - do entire corner as block
    # x left y left
    uₓₓ[1:ordx+retx,1:ordy+rety]          = 
        Stencil2D(u[1:ordx+adjx,1:ordy+adjy],Left,Left,Δx,Δy,
            cx[1:ordx+adjx,1:ordy+adjy],
            cy[1:ordx+adjx,1:ordy+adjy],ordx=ordx,ordy=ordy)
    # x left y right
    uₓₓ[1:ordx+retx,ny-ordy-rety+1:ny]      = 
        Stencil2D(u[1:ordx+adjx,ny-ordy-adjy:ny],Left,Right,Δx,Δy,
            cx[1:ordx+adjx,ny-ordy-adjy:ny],
            cy[1:ordx+adjx,ny-ordy-adjy:ny],ordx=ordx,ordy=ordy)
    # x right y left
    uₓₓ[nx-ordx-retx+1:nx,1:ordy+rety]      = 
        Stencil2D(u[nx-ordx-adjx:nx,1:ordy+adjy],Right,Left,Δx,Δy,
            cx[nx-ordx-adjx:nx,1:ordy+adjy],
            cy[nx-ordx-adjx:nx,1:ordy+adjy],ordx=ordx,ordy=ordy)
    # x right y right
    uₓₓ[nx-ordx-retx+1:nx,ny-ordy-rety+1:ny]  = 
        Stencil2D(u[nx-ordx-adjx:nx,ny-ordy-adjy:ny],Right,Right,Δx,Δy,
            cx[nx-ordx-adjx:nx,ny-ordy-adjy:ny],
            cy[nx-ordx-adjx:nx,ny-ordy-adjy:ny],ordx=ordx,ordy=ordy)

    return uₓₓ
end



"""
    Stencil2D
For internal nodes use
    [`Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;ordx=2,ordy=ordx)`](@ref)
For boundary nodes use either
    `Stencil2D(u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy;ordx=2,ordy=ordx)`
    `Stencil2D(u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy;ordx=2,ordy=ordx)`
For corners use
    `Stencil2D(u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy;ordx=2,ordy=ordx)`


Computes the 2D finite difference stencil at a given node
"""
function Stencil2D end
### Iternal node stencil
"""
    Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;ordx=2,ordy=ordx)

See [`Stencil2D`](@ref)
"""
function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;ordx=2,ordy=ordx)

    halfx = Int64(ordx/2) #half way
    halfy = Int64(ordy/2) #half way

    @sync @distributed for j = 1:ny
        for i = 1:nx
            uₓₓ[i,j] = SecondDerivative(u[i:i+ordx,j+halfx],
                    cx[i:i+ordx,j+halfx],Δx,Internal,order=ordx) + 
                SecondDerivative(u[i+halfy,j:j+ordy],
                    cy[i+halfy,j:j+ordy],Δy,Internal,order=ordy)
        end
    end
    return uₓₓ
end
### X boundaries
function Stencil2D(u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy;ordx=2,ordy=ordx)
    # side x internal y
    ordx == 2 ? adjx = 1 : adjx = ordx + Int64(ordx/2)
    adjy = Int64(ordy/2)+1 # half way
    
    uₓₓ = SecondDerivative(u[:,adjy],cy[:,adjy],Δx,xnode,order=ordx)
    for i = 1:adjx
        uₓₓ[i] += SecondDerivative(u[i,:],cx[i,:],Δy,Internal,order=ordy)
    end
    return uₓₓ

end
### Y boundaries
function Stencil2D(u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy;ordx=2,ordy=ordx)
    # internal x side y
    adjx = Int64(ordx/2)+1 # half way
    ordy == 2 ? adjy = 1 : adjy = ordy + Int64(ordy/2)

    uₓₓ = SecondDerivative(u[adjx,:],cy[adjx,:],Δy,ynode,order=ordy)
    for j = 1:adjy
        uₓₓ[j] += SecondDerivative(u[:,j],cx[:,j],Δx,Internal,order=ordx)
    end
    return uₓₓ

end
### Corners
function Stencil2D(u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy;ordx=2,ordy=ordx)
    # corner x-y
    ordx == 2 ? adjx = -1 : adjx = Int64(ordx/2)
    ordy == 2 ? adjy = -1 : adjy = Int64(ordy/2)

    uₓₓ = zeros(ordx+adjx,ordy+adjy)

    for j = 1:ordy+adjy # x direction - top/bottom
        uₓₓ[:,j] += SecondDerivative(u[:,j],cx[:,j],Δx,xnode,order=ordx)
    end
    for i = 1:ordx+adjx # y direction - left/right
        uₓₓ[i,:] += SecondDerivative(u[i,:],cy[i,:],Δy,ynode,order=ordy)
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
