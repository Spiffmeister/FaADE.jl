#======================================#
#====== Second Derivative Methods =====#
#======================================#
# Author: Dean Muir, Kenneth Duru



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
function D₂!(uₓₓ::AbstractMatrix{T},u::AbstractMatrix{T},c::AbstractMatrix{T},n::Integer,Δ::T,order::DerivativeOrder,α::T,dim::Integer) where T
    loopdir = SelectLoopDirection(dim)    
    for (cache,U,C) in zip(loopdir(uₓₓ),loopdir(u),loopdir(c))
        D₂!(cache,U,C,n,Δ,order,α)
    end
    uₓₓ
end
### 2D second derviative function
function D₂!(uₓₓ::AT,u::AT,cx::AT,cy::AT,
        nx::Integer,ny::Integer,Δx::T,Δy::T,
        order_x::Integer=2,order_y::Integer=order_x) where {T,AT}

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


    # SecondDerivativeInternal!(uₓₓ,u,cx,cy,Δx,Δy,nx,ny,order_x,0.0)
    SecondDerivativeInternal!(uₓₓ,u,cx,cy,Δx,Δy,nx,ny,DerivativeOrder{2}(),0.0)

    # for (A,B,C) in zip(eachcol(uₓₓ),eachcol(u),eachcol(cx))
    #     D₂!(A,B,C,nx,Δx,DerivativeOrder{2}(),0.0)
    # end
    # for (A,B,C) in zip(eachrow(uₓₓ),eachrow(u),eachrow(cy))
    #     D₂!(A,B,C,ny,Δy,DerivativeOrder{2}(),1.0)
    # end

    
    ### Boundary nodes - avoiding corners
    # left and right x boundaries, for a matrix zeros(nx,ny) this is the 'top' and 'bottom' boundaries
    # x left (top)
    # uₓₓ[1:rx, ry+1:ny-ry] = 
    #     SecondDerivativeBoundaryStencil!(uₓₓ[1:rx, ry+1:ny-ry],
    #         u[1:ix, inty:ny-inty+1],
    #         Left,Internal,Δx,Δy,
    #         cx[1:order_x+adjx, inty:ny-inty+1],
    #         cy[1:order_x+adjx, inty:ny-inty+1],
    #         rx,ny-2ry,
    #         order_x=order_x,order_y=order_y,
    #         Xrng=1:rx)

    SecondDerivativeBoundaryStencila!(uₓₓ,
        u,cx,cy,
        Δx,Δy,
        rx,ny,
        Left,Internal,
        DerivativeOrder{2}(),DerivativeOrder{2}())
    # x right (bottom)
    # println("u:",(nx-ix+1:nx,inty:ny-inty+1),"  uxx:",(nx-rx+1:nx,ry+1:ny-ry))
    # uₓₓ[nx-rx+1:nx, ry+1:ny-ry] =
    #     SecondDerivativeBoundaryStencil!(uₓₓ[nx-rx+1:nx, ry+1:ny-ry],
    #         u[nx-ix+1:nx,  inty:ny-inty+1],
    #         Right,Internal,Δx,Δy,
    #         cx[nx-ix+1:nx, inty:ny-inty+1],
    #         cy[nx-ix+1:nx, inty:ny-inty+1],
    #         rx,ny-2ry,
    #         order_x=order_x,order_y=order_y,
    #         Yrng=(1:ny-2ry))

    SecondDerivativeBoundaryStencila!(uₓₓ,
        u,cx,cy,
        Δx,Δy,
        rx,ny,
        Right,Internal,
        DerivativeOrder{2}(),DerivativeOrder{2}())

    # left and right y boundaries, for a matrix zeros(nx,ny) this is the 'left' and 'right' boundaries
    # y left (left)
    uₓₓ[rx+1:nx-rx, 1:ry] = 
        SecondDerivativeBoundaryStencil!(uₓₓ[rx+1:nx-rx, 1:ry],
            u[intx:nx-intx+1, 1:iy],
            Internal,Left,Δx,Δy,
            cx[intx:nx-intx+1,1:iy],
            cy[intx:nx-intx+1,1:iy],
            nx-2rx,ry,
            order_x=order_x,order_y=order_y)
    # y right (right)
    uₓₓ[rx+1:nx-rx, ny-ry+1:ny] = 
        SecondDerivativeBoundaryStencil!(uₓₓ[rx+1:nx-rx, ny-ry+1:ny],
            u[intx:nx-intx+1, ny-iy+1:ny],
            Internal,Right,Δx,Δy,
            cx[intx:nx-intx+1,ny-iy+1:ny],
            cy[intx:nx-intx+1,ny-iy+1:ny],
            nx-2rx,ry,
            order_x=order_x,order_y=order_y)

    # Corners - do entire corner as block
    # x left y left
    uₓₓ[1:rx, 1:order_y+rety]          = 
        SecondDerivativeCornerStencil!(uₓₓ[1:rx, 1:order_y+rety],
            u[1:order_x+adjx, 1:order_y+adjy],Left,Left,Δx,Δy,
            cx[1:order_x+adjx,1:order_y+adjy],
            cy[1:order_x+adjx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x left y right
    uₓₓ[1:rx, ny-ry+1:ny]    =
        SecondDerivativeCornerStencil!(uₓₓ[1:rx, ny-ry+1:ny],
            u[1:order_x+adjx, ny-iy+1:ny],
            Left,Right,Δx,Δy,
            cx[1:order_x+adjx,ny-iy+1:ny],
            cy[1:order_x+adjx,ny-iy+1:ny],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x right y left
    uₓₓ[nx-ry+1:nx, 1:ry]      = 
        SecondDerivativeCornerStencil!(uₓₓ[nx-ry+1:nx, 1:ry],
            u[nx-ix+1:nx, 1:order_y+adjy],
            Right,Left,Δx,Δy,
            cx[nx-ix+1:nx,1:order_y+adjy],
            cy[nx-ix+1:nx,1:order_y+adjy],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    # x right y right
    uₓₓ[nx-ry+1:nx, ny-ry+1:ny]  = 
        SecondDerivativeCornerStencil!(uₓₓ[nx-ry+1:nx, ny-ry+1:ny],
        u[nx-ix+1:nx, ny-iy+1:ny],
            Right,Right,Δx,Δy,
            cx[nx-ix+1:nx,ny-iy+1:ny],
            cy[nx-ix+1:nx,ny-iy+1:ny],
            order_x+adjx,order_y+adjy,
            order_x=order_x,order_y=order_y)
    
    uₓₓ
end


"""
    SecondDerivativeBoundaryStencil!

2D second derivative along a boundary.
"""
function SecondDerivativeBoundaryStencil! end
### X boundaries
@views function SecondDerivativeBoundaryStencil!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=1:BoundaryNodeOutput(order_x), 
            Yrng=(BoundaryNodeInput(order_y),ny-BoundaryNodeInput(order_y)+1))
    # side x internal y

    halfy = halforder(order_y) # half way
    (xnode == Right) & (order_x > 2) ? offset = halforder(order_x) : offset = 0 #offset for right hand side
    Inx = BoundaryNodeInput(order_x) # for order > 2, u has more rows

    for j = 1:ny #boundaries nxodes[1]:nxodes[2]
        SecondDerivativeBoundary!(uₓₓ[Xrng,j],u[1:Inx,j+halfy],cx[1:Inx,j+halfy],Δx,xnode,order_x,0.0)
        for i = Xrng #cross terms
            """ j:j+order_y in u because u[,1] is uₓₓ[,-2] """
            uₓₓ[i,j] += SecondDerivativeInternal(u[i+offset,j:j+order_y],cy[i+offset,j:j+order_y],Δy,Internal,order=order_y)
        end
    end
    return uₓₓ
end
function SecondDerivativeBoundaryStencila!(uₓₓ::AbstractArray{T},
        u::AbstractArray{T},cx::AbstractArray{T},cy::AbstractArray{T},
        Δx::T,Δy::T,
        nx::Integer,ny::Integer,
        xnode::NodeType,::NodeType{:Internal},
        order_x=2,order_y=order_x) where T
    # side x internal y

    for (A,B,C) in zip(eachrow(uₓₓ),eachrow(u),eachrow(cx))
        SecondDerivativeBoundary!(A,B,C,Δx,xnode,DerivativeOrder{2}(),0.0)
    end
    for j = 2:ny-1
        # SecondDerivativeBoundary!(view(uₓₓ,:,j),u[1:Inx,j+halfy],cx[1:Inx,j+halfy],Δx,xnode,DerivativeOrder{2}(),0.0)
        for i = 1:nx
            SecondDerivativeInternal!(uₓₓ,u,cy,Δy,Internal,DerivativeOrder{2}(),i,j,1.0)
        end
    end


    uₓₓ
end
function SecondDerivativeBoundaryStencilb!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,::NodeType{:Internal},Δx,Δy,cx,cy,nx,ny;
    order_x=2,order_y=order_x,
        Xrng::UnitRange=1:BoundaryNodeOutput(order_x), 
        Yrng=(BoundaryNodeInput(order_y),ny-BoundaryNodeInput(order_y)+1))
    # side x internal y

    halfy = halforder(order_y) # half way
    (xnode == Right) & (order_x > 2) ? offset = halforder(order_x) : offset = 0 #offset for right hand side
    Inx = BoundaryNodeInput(order_x) # for order > 2, u has more rows
    
    for j = 7:15
        SecondDerivativeBoundary!(view(uₓₓ,:,j),view(u,:,j),view(cx,:,j),Δx,xnode,order_x,0.0)
    end
    for i = 16:21
        SecondDerivativeInternal!(view(uₓₓ,i,:),view(u,i,:),view(cy,i,:),Δy,21,order_y,0.0)
    end
    uₓₓ
end
### Y boundaries
@views function SecondDerivativeBoundaryStencil!(uₓₓ::AbstractMatrix,u::AbstractMatrix,::NodeType{:Internal},ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=BoundaryNodeInput(order_x):(nx-BoundaryNodeInput(order_x)+1),
            Yrng::UnitRange=1:BoundaryNodeOutput(order_y))
        
    # internal x side y
    halfx = halforder(order_x) # half way
    (ynode == Right) & (order_y > 2) ? offset = halforder(order_y) : offset = 0
    Iny = BoundaryNodeInput(order_y)

    for i = 1:nx
        SecondDerivativeBoundary!(uₓₓ[i,Yrng],u[i+halfx,1:Iny],cy[i+halfx,1:Iny],Δy,ynode,order_y,0.0)
        for j = Yrng
            uₓₓ[i,j] += SecondDerivativeInternal(u[i:i+order_y,j+offset],cx[i:i+order_y,j+offset],Δx,Internal,order=order_x)
        end
    end
    return uₓₓ
end
### Corners
"""
    SecondDerivativeCornerStencil!
2D second derivative on a corner.
"""
@views function SecondDerivativeCornerStencil!(uₓₓ::AbstractMatrix,u::AbstractMatrix,xnode::NodeType,ynode::NodeType,Δx,Δy,cx,cy,nx,ny;
        order_x=2,order_y=order_x,
            Xrng::UnitRange=1:BoundaryNodeOutput(order_x),
            Yrng::UnitRange=1:BoundaryNodeOutput(order_y))
    # corner x-y

    (xnode == Right)  & (order_x > 2) ? oy = halforder(order_x) : oy = 0 #x node offset for left boundary
    (ynode == Right)  & (order_y > 2) ? ox = halforder(order_y) : ox = 0 #y node offset for left boundary
    for j = Yrng # x direction - top/bottom
        SecondDerivativeBoundary!(uₓₓ[Xrng,j],u[:,j+ox],cx[:,j+ox],Δx,xnode,order_x,0.0)
    end
    for i = Xrng # y direction - left/right
        uₓₓ[i,Yrng] += SecondDerivative(u[i+oy,:],cy[i+oy,:],Δy,ynode,order=order_y)
    end
    return uₓₓ
end

