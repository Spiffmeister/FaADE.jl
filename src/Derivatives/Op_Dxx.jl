
"""
    Dₓₓ

    Dₓₓ(u::AbstractVector,c::Vector,n::Int64,Δx::Float64;order::Int64=2)
    Dₓₓ(u::Matrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix;dim::Int64=1,order::Int64=2)
    Dₓₓ(u::AbstractMatrix,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,cx::AbstractMatrix,cy::AbstractMatrix;
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
        Diff(uₓₓ,u,cx,cy) = Dₓₓ!(uₓₓ,u,cx,cy, nx,ny,Δx,Δy,order)
        return Diff
    end
end


"""
    Dₓₓ!

    Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},n::Int64,Δx::Float64;order::Int64=2)
    Dₓₓ!(uₓₓ::AbstractMatrix{Float64},u::AbstractMatrix{Float64},nx::Int64,ny::Int64,Δ::Float64,c::AbstractMatrix{Float64};dim::Int64=1,order::Int64=2)
    Dₓₓ!(uₓₓ::AbstractVector{Float64},u::AbstractVector{Float64},c::AbstractVector{Float64},m::Int64,Δ::Float64;order::Int64=2)

Mutable function for 1D and 2D second derviative SBP operator
"""
function Dₓₓ! end
### 1D second derviative mutable function
@inline function Dₓₓ!(uₓₓ::AbstractVector,u::AbstractVector,c::AbstractVector,n::Int64,Δx::Float64,order::Int64)
    #TODO HIGHER ORDER METHODS
    # half = halforder(order)
    # adj = BoundaryNodeOutput(order)
    adj = BoundaryNodeInput(order)
    ret = BoundaryNodeOutput(order)
    
    SecondDerivativeInternal1D!(uₓₓ,u,c,Δx,n,order)

    # println("Internal   ",uₓₓ)
    SecondDerivative!(uₓₓ,u[1:adj],c[1:adj],Δx,Left,order=order)
    # println("left   ",uₓₓ)
    SecondDerivative!(uₓₓ,u[n-adj+1:n],c[n-adj+1:n],Δx,Right,order=order)
    # println("right    ",uₓₓ)
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

    rx = BoundaryNodeOutput(order_x)
    ry = BoundaryNodeOutput(order_y)

    ix = BoundaryNodeInput(order_x)
    iy = BoundaryNodeInput(order_y)


    SecondDerivativeInternal2D!(uₓₓ,u,cx,cy,Δx,Δy,nx,ny,order_x)

    ### Boundary nodes - avoiding corners
    # left and right x boundaries, for a matrix zeros(nx,ny) this is the 'top' and 'bottom' boundaries
    # x left (top)
    uₓₓ[1:rx, ry+1:ny-ry] = 
        Stencil2D!(uₓₓ[1:rx, ry+1:ny-ry],
            u[1:ix, inty:ny-inty+1],
            Left,Internal,Δx,Δy,
            cx[1:order_x+adjx, inty:ny-inty+1],
            cy[1:order_x+adjx, inty:ny-inty+1],
            rx,ny-2ry,
            order_x=order_x,order_y=order_y,
            Xrng=1:rx)
    # x right (bottom)
    # println("u:",(nx-ix+1:nx,inty:ny-inty+1),"  uxx:",(nx-rx+1:nx,ry+1:ny-ry))
    uₓₓ[nx-rx+1:nx, ry+1:ny-ry] =
        Stencil2D!(uₓₓ[nx-rx+1:nx, ry+1:ny-ry],
            u[nx-ix+1:nx,  inty:ny-inty+1],
            Right,Internal,Δx,Δy,
            cx[nx-ix+1:nx, inty:ny-inty+1],
            cy[nx-ix+1:nx, inty:ny-inty+1],
            rx,ny-2ry,
            order_x=order_x,order_y=order_y,
            Yrng=(1:ny-2ry))

    # left and right y boundaries, for a matrix zeros(nx,ny) this is the 'left' and 'right' boundaries
    # y left (left)
    uₓₓ[rx+1:nx-rx, 1:ry] = 
        Stencil2D!(uₓₓ[rx+1:nx-rx, 1:ry],
            u[intx:nx-intx+1, 1:iy],
            Internal,Left,Δx,Δy,
            cx[intx:nx-intx+1,1:iy],
            cy[intx:nx-intx+1,1:iy],
            nx-2rx,ry,
            order_x=order_x,order_y=order_y)
    # y right (right)
    uₓₓ[rx+1:nx-rx, ny-ry+1:ny] = 
        Stencil2D!(uₓₓ[rx+1:nx-rx, ny-ry+1:ny],
            u[intx:nx-intx+1, ny-iy+1:ny],
            Internal,Right,Δx,Δy,
            cx[intx:nx-intx+1,ny-iy+1:ny],
            cy[intx:nx-intx+1,ny-iy+1:ny],
            nx-2rx,ry,
            order_x=order_x,order_y=order_y)

    # Corners - do entire corner as block
    # x left y left
    uₓₓ[1:rx, 1:order_y+rety]          = 
        Stencil2D!(uₓₓ[1:rx, 1:order_y+rety],
            u[1:order_x+adjx, 1:order_y+adjy],Left,Left,Δx,Δy,
            cx[1:order_x+adjx,1:order_y+adjy],
            cy[1:order_x+adjx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x left y right
    uₓₓ[1:rx, ny-ry+1:ny]    =
        Stencil2D!(uₓₓ[1:rx, ny-ry+1:ny],
            u[1:order_x+adjx, ny-iy+1:ny],
            Left,Right,Δx,Δy,
            cx[1:order_x+adjx,ny-iy+1:ny],
            cy[1:order_x+adjx,ny-iy+1:ny],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x right y left
    uₓₓ[nx-ry+1:nx, 1:ry]      = 
        Stencil2D!(uₓₓ[nx-ry+1:nx, 1:ry],
            u[nx-ix+1:nx, 1:order_y+adjy],
            Right,Left,Δx,Δy,
            cx[nx-ix+1:nx,1:order_y+adjy],
            cy[nx-ix+1:nx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x right y right
    uₓₓ[nx-ry+1:nx, ny-ry+1:ny]  = 
        Stencil2D!(uₓₓ[nx-ry+1:nx, ny-ry+1:ny],
        u[nx-ix+1:nx, ny-iy+1:ny],
            Right,Right,Δx,Δy,
            cx[nx-ix+1:nx,ny-iy+1:ny],
            cy[nx-ix+1:nx,ny-iy+1:ny],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # println()
    # println()
    # return uₓₓ

    uₓₓ
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
function Stencil2D! end
### X boundaries
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=1:BoundaryNodeOutput(order_x), 
            Yrng=(BoundaryNodeInput(order_y),ny-BoundaryNodeInput(order_y)+1))
    # side x internal y

    halfy = halforder(order_y) # half way
    (xnode == Right) & (order_x > 2) ? offset = halforder(order_x) : offset = 0 #offset for right hand side
    Inx = BoundaryNodeInput(order_x) # for order > 2, u has more rows

    for j = 1:ny #boundaries nxodes[1]:nxodes[2]
        SecondDerivative!(uₓₓ[Xrng,j],u[1:Inx,j+halfy],cx[1:Inx,j+halfy],Δx,xnode,order=order_x)
        for i = Xrng #cross terms
            """ j:j+order_y in u because u[,1] is uₓₓ[,-2] """
            uₓₓ[i,j] += SecondDerivative(u[i+offset,j:j+order_y],cy[i+offset,j:j+order_y],Δy,Internal,order=order_y)
        end
    end
    return uₓₓ
end
### Y boundaries
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=BoundaryNodeInput(order_x):(nx-BoundaryNodeInput(order_x)+1),
            Yrng::UnitRange=1:BoundaryNodeOutput(order_y))
        
    # internal x side y
    halfx = halforder(order_x) # half way
    (ynode == Right) & (order_y > 2) ? offset = halforder(order_y) : offset = 0
    Iny = BoundaryNodeInput(order_y)

    for i = 1:nx
        SecondDerivative!(uₓₓ[i,Yrng],u[i+halfx,1:Iny],cy[i+halfx,1:Iny],Δy,ynode,order=order_y)
        for j = Yrng
            uₓₓ[i,j] += SecondDerivative(u[i:i+order_y,j+offset],cx[i:i+order_y,j+offset],Δx,Internal,order=order_x)
        end
    end
    return uₓₓ
end
### Corners
@views function Stencil2D!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=1:BoundaryNodeOutput(order_x),
            Yrng::UnitRange=1:BoundaryNodeOutput(order_y))
    # corner x-y

    (xnode == Right)  & (order_x > 2) ? oy = halforder(order_x) : oy = 0 #x node offset for left boundary
    (ynode == Right)  & (order_y > 2) ? ox = halforder(order_y) : ox = 0 #y node offset for left boundary
    for j = Yrng # x direction - top/bottom
        SecondDerivative!(uₓₓ[Xrng,j],u[:,j+ox],cx[:,j+ox],Δx,xnode,order=order_x)
    end
    for i = Xrng # y direction - left/right
        uₓₓ[i,Yrng] += SecondDerivative(u[i+oy,:],cy[i+oy,:],Δy,ynode,order=order_y)
    end
    return uₓₓ
end

