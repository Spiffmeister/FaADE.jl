
"""
    Dₓₓ

The second derivative SBP operator. There are three call options:
    1. Dₓₓ(u::AbstractVector,c::Vector,n::Int64,Δx::Float64;order::Int64=2)
    2. Dₓₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix;dim::Int64=1,order::Int64=2)
    3. Dₓₓ(u::AbstractMatrix,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,cx::AbstractMatrix,cy::AbstractMatrix;
        order_x::Int64=2,order_y::Int64=order_x)

1. 1D second derviative SBP operator
    - Returns a `Vector{Float64}` of length `n`.

2. 1D second derivative SBP operator along the dimension given by `dim`
    - dim ∈ [1,2]
        - If `dim==1` then takes derivative along rows (`u[:,i]`)
        - If `dim==2` then takes derivative along columns (`u[i,:]`)
    - Δ should be the grid spacing along the dimension requested
    - Returns a `Matrix{Float64}` of size `nx × ny`.

3. Computes the second derivative operator using stencils across an enture 2D domain

- order_x and order_y are the order of the derivative approximations



This is also available is an iterator [`Dₓₓ!`](@ref)

Internally uses [`SecondDerivative`](@ref)
"""
function Dₓₓ end
function Dₓₓ(u::AbstractVector,c::AbstractVector,n::Int64,Δx::Float64;order::Int64=2)
    # Function call for the 1D 2nd derivative SBP operator
    uₓₓ = zeros(eltype(u),size(u))
    Dₓₓ!(uₓₓ,u,c,n,Δx,order)
    return uₓₓ
end
### 2D second derivative operator
function Dₓₓ(u::AbstractMatrix,n::Int64,Δ::Float64,c::AbstractMatrix,dim::Int64=1;order::Int64=2)
    # Multidimensional call for 2nd derivative SBP operator
    uₓₓ = zeros(eltype(u),size(u))
    Dₓₓ!(uₓₓ,u,c,n,Δ,dim,order)
    return uₓₓ
end
function Dₓₓ(u::AbstractMatrix,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,cx::AbstractMatrix,cy::AbstractMatrix;
    order_x::Int64=2,order_y::Int64=order_x)
    # 2D 2nd derivative SBP operator
    uₓₓ = SharedArray(zeros(Float64,nx,ny))
    Dₓₓ!(uₓₓ,u,nx,ny,Δx,Δy,cx,cy,order_x,order_y)
    return uₓₓ
end






"""
    generate_Derivative

Returns a function `Diff(uₓₓ,u,c...)` which computes the derivative in 1 or 2 dimensions
"""
function generate_Derivative end
function generate_Derivative(n::Int64,Δx::Float64,order::Int64)
    let n = n, Δx=Δx, order=order
        Diff(uₓₓ,u,c) = Dₓₓ!(uₓₓ,u,c,n,Δx,order)
        return Diff
    end
end
function generate_Derivative(nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,order::Int64)
    let nx=nx, ny=ny,
            Δx=Δx, Δy=Δy,
            order=order
        Diff(uₓₓ,u,cx,cy) = Dₓₓ!(uₓₓ,u,cx,cy, nx,ny,Δx,Δy,order,order)
        return Diff
    end
end


"""
    Dₓₓ!



    Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
or
    Dₓₓ!(uₓₓ::AbstractMatrix{Float64},u::AbstractMatrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix{Float64};dim::Int64=1,order::Int64=2)
or
    Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},m::Int64,Δ::Float64;order::Int64=2)

Mutable function for 1D and 2D second derviative SBP operator
"""
function Dₓₓ! end
### 1D second derviative mutable function
@inline function Dₓₓ!(uₓₓ::AbstractVector,u::AbstractVector,c::AbstractVector,n::Int64,Δx::Float64,order::Int64)
    #TODO HIGHER ORDER METHODS
    # half = halforder(order)
    # adj = BoundaryNodeOutput(order)
    
    SecondDerivativeInternal1D!(@views(uₓₓ[2:n-1]),u,c,Δx,n)

    SecondDerivative!(uₓₓ,u[1:order],c[1:order],Δx,Left,order=order)
    SecondDerivative!(uₓₓ,u[n-order+1:n],c[n-order+1:n],Δx,Right,order=order)
    uₓₓ
end
### Multidimensional second derivative SBP operator
function Dₓₓ!(uₓₓ::AbstractMatrix,u::AbstractMatrix,c::AbstractMatrix,n::Int64,Δ::Float64,dim::Int64=1,order::Int64=2)
    loopdir = SelectLoopDirection(dim)    
    for (T,U,C) in zip(loopdir(uₓₓ),loopdir(u),loopdir(c))
        Dₓₓ!(T,U,C,n,Δ,order)
    end
    return uₓₓ
end
### 2D second  derviative mutable function
function Dₓₓ!(uₓₓ::AbstractMatrix,u::AbstractMatrix,cx::AbstractMatrix,cy::AbstractMatrix,
        nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,
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


    SecondDerivativeInternal2D!(@views(uₓₓ[2:nx-1,2:ny-1]),u,cx,cy,Δx,Δy,nx,ny)

    ### Boundary nodes - avoiding corners
    # left and right x boundaries, for a matrix zeros(nx,ny) this is the 'top' and 'bottom' boundaries
    # x left (top)
    uₓₓ[1:order_x+retx,inty+halfy:ny-(inty+halfy-1)] = 
        Stencil2D!(uₓₓ[1:order_x+retx,inty+halfy:ny-(inty+halfy-1)],
            u[1:order_x+adjx,inty:ny-inty+1],Left,Internal,Δx,Δy,
            cx[1:order_x+adjx,inty:ny-inty+1],
            cy[1:order_x+adjx,inty:ny-inty+1],
            order_x+retx,ny-2(order_y+rety),order_x=order_x,order_y=order_y)
    # x right (bottom)
    uₓₓ[nx-(order_x+retx-1):nx,order_y+rety+1:ny-order_y-rety] =
        Stencil2D!(uₓₓ[nx-(order_x+retx-1):nx,inty+halfy:ny-(inty+halfy-1)],
            u[nx-order_x-adjx:nx,inty:ny-inty+1],Right,Internal,Δx,Δy,
            cx[nx-order_x-adjx:nx,inty:ny-inty+1],
            cy[nx-order_x-adjx:nx,inty:ny-inty+1],
            order_x+retx,ny-2(order_y+rety),order_x=order_x,order_y=order_y)

    # left and right y boundaries, for a matrix zeros(nx,ny) this is the 'left' and 'right' boundaries
    # y left (left)
    uₓₓ[intx+halfx:nx-(intx+halfx-1),1:order_y+retx] = 
        Stencil2D!(uₓₓ[intx+halfx:nx-(intx+halfx-1),1:order_y+retx],
        u[intx:nx-intx+1,1:order_y+adjy],Internal,Left,Δx,Δy,
        cx[intx:nx-intx+1,1:order_y+adjy],
        cy[intx:nx-intx+1,1:order_y+adjy],
        nx-2(order_x+retx),order_y+rety,order_x=order_x,order_y=order_y)
    # y right (right)
    uₓₓ[order_x+retx+1:nx-order_y-retx,ny-(order_y+rety-1):ny] = 
        Stencil2D!(uₓₓ[order_x+retx+1:nx-order_y-retx,ny-(order_y+rety-1):ny],
            u[intx:nx-intx+1,ny-order_y-adjy:ny],Internal,Right,Δx,Δy,
            cx[intx:nx-intx+1,ny-order_y-adjy:ny],
            cy[intx:nx-intx+1,ny-order_y-adjy:ny],
            nx-2(order_x+retx),order_y+rety,order_x=order_x,order_y=order_y)

    # Corners - do entire corner as block
    # x left y left
    uₓₓ[1:order_x+retx,1:order_y+rety]          = 
        Stencil2D!(uₓₓ[1:order_x+retx,1:order_y+rety],u[1:order_x+adjx,1:order_y+adjy],Left,Left,Δx,Δy,
            cx[1:order_x+adjx,1:order_y+adjy],
            cy[1:order_x+adjx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,order_x=order_x,order_y=order_y)
    # x left y right
    uₓₓ[1:order_x+retx,ny-order_y-rety+1:ny]    =
        Stencil2D!(uₓₓ[1:order_x+retx,1:order_y+rety],u[1:order_x+adjx,ny-order_y-adjy:ny],Left,Right,Δx,Δy,
            cx[1:order_x+adjx,ny-order_y-adjy:ny],
            cy[1:order_x+adjx,ny-order_y-adjy:ny],
            order_x+adjx,order_y+adjy,order_x=order_x,order_y=order_y)
    # x right y left
    uₓₓ[nx-order_x-retx+1:nx,1:order_y+rety]      = 
        Stencil2D!(uₓₓ[1:order_x+retx,1:order_y+rety],u[nx-order_x-adjx:nx,1:order_y+adjy],Right,Left,Δx,Δy,
            cx[nx-order_x-adjx:nx,1:order_y+adjy],
            cy[nx-order_x-adjx:nx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,order_x=order_x,order_y=order_y)
    # x right y right
    uₓₓ[nx-order_x-retx+1:nx,ny-order_y-rety+1:ny]  = 
        Stencil2D!(uₓₓ[1:order_x+retx,1:order_y+rety],u[nx-order_x-adjx:nx,ny-order_y-adjy:ny],Right,Right,Δx,Δy,
            cx[nx-order_x-adjx:nx,ny-order_y-adjy:ny],
            cy[nx-order_x-adjx:nx,ny-order_y-adjy:ny],
            order_x+adjx,order_y+adjy,order_x=order_x,order_y=order_y)
    # return uₓₓ
end



"""
    Stencil2D!
For internal nodes use
    1. `Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;order_x=2,order_y=order_x)`
For boundary nodes use either
    2. `Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy;order_x=2,order_y=order_x)`
    3. `Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy;order_x=2,order_y=order_x)`
For corners use
    4. `Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy;order_x=2,order_y=order_x)`

Use 1. for internal nodes, 2. and 3. for non-corner boundaries, 4. for corners.

Computes the 2D finite difference stencil at a given node
"""
function Stencil2D end

function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},::NodeType{:Internal},Δx,Δy,cx,cy;
    order_x=2,order_y=order_x)

    SecondDerivativeInd!(uₓₓ,u,cx,Δx,Internal)
    SecondDerivativeIndAdd!(uₓₓ,u,cy,Δy,Internal)
    uₓₓ
end

### X boundaries
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x)
    # side x internal y
    order_x == 2 ? adjx = 1 : adjx = order_x + Int64(order_x/2)
    halfy = Int64(order_y/2) # half way

    # println(xnode)
    for j = 1:ny
        # println("before",uₓₓ[1:nx,j])
        # uₓₓ[1:nx,j] = SecondDerivative(u[:,j+halfy],cx[:,j+halfy],Δx,xnode,order=order_x)
        SecondDerivative!(uₓₓ[1:nx,j],u[:,j+halfy],cx[:,j+halfy],Δx,xnode,order=order_x)
        # println("after",uₓₓ[1:nx,j])
        for i = 1:adjx
            uₓₓ[i,j] += SecondDerivative(u[i,j:j+order_y],cy[i,j:j+order_y],Δy,Internal,order=order_y)
        end
    end
    return uₓₓ
end
### Y boundaries
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x)
    # internal x side y
    halfx = Int64(order_x/2) # half way
    order_y == 2 ? adjy = 1 : adjy = order_y + Int64(order_y/2)

    for i = 1:nx
        # uₓₓ[i,1:ny] = SecondDerivative(u[i+halfx,:],cy[i+halfx,:],Δy,ynode,order=order_y)
        SecondDerivative!(uₓₓ[i,:],u[i+halfx,:],cy[i+halfx,:],Δy,ynode,order=order_y)
        for j = 1:adjy
            uₓₓ[i,j] += SecondDerivative(u[i:i+order_y,j],cx[i:i+order_y,j],Δx,Internal,order=order_x)
        end
    end
    return uₓₓ
end
### Corners
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x)
    # corner x-y
    order_x == 2 ? adjx = -1 : adjx = Int64(order_x/2)
    order_y == 2 ? adjy = -1 : adjy = Int64(order_y/2)

    for j = 1:order_y+adjy # x direction - top/bottom
        # uₓₓ[:,j] = SecondDerivative(u[:,j],cx[:,j],Δx,xnode,order=order_x)
        SecondDerivative!(uₓₓ[:,j],u[:,j],cx[:,j],Δx,xnode,order=order_x)
    end
    for i = 1:order_x+adjx # y direction - left/right
        # uₓₓ[i,:] .+= SecondDerivative(u[i,:],cy[i,:],Δy,ynode,order=order_y)
        SecondDerivative!(uₓₓ[i,:],u[i,:],cy[i,:],Δy,ynode,order=order_y)
    end
    return uₓₓ
end

