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
    D₁!(uₓ,u,length(u),Δx,order)
    return uₓ
end
"""
    D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;order::Integer=2)
2D implementation of ``D_x`` operator
"""
function D₁(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;
        order::Integer=2) where T
    uₓ = zeros(T,size(u))
    D₁!(uₓ,u,nx,ny,Δx,Δy,order)
    return uₓ
end


"""
    D₁!
1D and 2D in place first derivative operator.

See also [`FirstDerivativeBoundary!`](@ref) and [`FirstDerivativeInternal!`](@ref).
"""
function D₁!() end
"""
    D₁!(uₓ::AbstractVector{T},u::AbstractVector{T},n::Integer,Δx::T,order::Integer)
1D [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractVector{T},u::AbstractVector{T},n::Integer,Δx::T,
        order::Integer) where T
    FirstDerivativeBoundary!(uₓ,u,Δx,n,Left,order)
    FirstDerivativeInternal!(uₓ,u,Δx,n,order)
    FirstDerivativeBoundary!(uₓ,u,Δx,n,Right,order)
end
"""
    function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,order::Integer,dim::Integer)
1D implementation for 2D problems for [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,
        order::Integer,dim::Integer) where T
    loopdir = SelectLoopDirection(dim)
    for (cache,U) in zip(loopdir(uₓ),loopdir(u))
        D₁!(cache,U,n,Δ,order)
    end
    uₓ
end
"""
    function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,order::Integer)
2D [`D₁!`](@ref).
"""
function D₁!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,
        order::Integer) where T
    
    order == 2 ? ret = 1 : ret = order
    order == 2 ? nin = 2 : nin = order + halforder(order)
    order == 2 ? ein = 1 : ein = halforder(order)

    

    FirstDerivativeInternal!(uₓ,u,Δx,Δy,nx,ny,order)

    #Boundary nodes
    FirstDerivativeBoundaryStencil!(uₓ[1:ret,order+1:ny-order+1],
        u[1:nin,ein:ny-ein+1],
        Left,Internal,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ[n-ret+1:n,order+1:ny-order+1],
        u[nx-nin+1:nx,ein:ny-ein+1],
        Right,Internal,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ[order+1:nx-order+1,1:ret],
        u[ein:nx-ein+1,1:order],
        Internal,Left,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ[order+1:nx-order+1,n-ret+1:n],
        u[ein:nx-ein+1,ny-nin+1:nx],
        Internal,Right,Δx,Δy,nx,ny,order)

    #Corner nodes
    FirstDerivativeCornerStencil!(uₓ[1:ret,1:ret],
        u[1:nin,1:nin],
        Left,Left,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ[1:ret,ny-ret+1:ny],
        u[1:nin,ny-nin+1:ny],
        Left,Right,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ[nx-ret+1:nx,ny-ret+1:ny],
        u[nx-nin+1:nx,ny-nin+1:ny],
        Right,Right,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ[nx-ret+1:nx,1:ret],
        u[nx-nin+1:nx,1:nin],
        Right,Left,Δx,Δy,nx,ny,order)
    
    uₓ
end


"""
    FirstDerivativeBoundaryStencil!
2D first derivative on a boundary.
"""
function FirstDerivativeBoundaryStencil!() end
"""
    FirstDerivativeBoundaryStencil!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        xnode::NodeType,::NodeType{:Internal},
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer,
        Xrng::UnitRange=1:order,Yrng::UnitRange=1:order)
"""
function FirstDerivativeBoundaryStencil!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        xnode::NodeType,::NodeType{:Internal},
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer,
        Xrng::UnitRange=1:order,Yrng::UnitRange=1:order) where T
    
    halfy = halforder(order_y) # half way
    (xnode == Right) & (order_x > 2) ? offset = halforder(order_x) : offset = 0 #offset for right hand side
    Inx = BoundaryNodeInput(order_x) # for order > 2, u has more rows

    for j = 1:ny #boundaries nxodes[1]:nxodes[2]
        FirstDerivativeBoundary!(uₓ[Xrng,j],u[1:Inx,j+halfy],Δx,nx,xnode,order)
        for i = Xrng #cross terms
            """ j:j+order_y in u because u[,1] is uₓₓ[,-2] """
            uₓₓ[i,j] += FirstDerivativeInternal(u[i+offset,j:j+order_y],Δy,ny,order)
        end
    end
    return uₓₓ
end
"""
    FirstDerivativeBoundaryStencil!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        ::NodeType{:Internal},ynode::NodeType,
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer,
        Xrng::UnitRange=1:order,Yrng::UnitRange=1:order)
"""
function FirstDerivativeBoundaryStencil!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        ::NodeType{:Internal},ynode::NodeType,
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer,
        Xrng::UnitRange=1:order,Yrng::UnitRange=1:order) where T

    halfx = halforder(order_x) # half way
    (ynode == Right) & (order_y > 2) ? offset = halforder(order_y) : offset = 0
    Iny = BoundaryNodeInput(order_y)

    for i = 1:nx
        FirstDerivativeBoundary!(uₓ[i,Yrng],u[i+halfx,1:Iny],Δy,ny,ynode,order)
        for j = Yrng
            uₓₓ[i,j] += FirstDerivativeInternal(u[i:i+order_y,j+offset],Δx,nx,order)
        end
    end
    return uₓₓ
end

"""
    FirstDerivativeCornerStencil!
2D first derivative on a corner.
"""
function FirstDerivativeCornerStencil!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        xnode::NodeType,ynode::NodeType,
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer,
        Xrng::UnitRange=1:order,
        Yrng::UnitRange=1:order) where T

    (xnode == Right)  & (order_x > 2) ? oy = halforder(order_x) : oy = 0 #x node offset for left boundary
    (ynode == Right)  & (order_y > 2) ? ox = halforder(order_y) : ox = 0 #y node offset for left boundary
    for j = Yrng # x direction - top/bottom
        FirstDerivativeBoundary!(uₓ[Xrng,j],u[:,j+ox],Δx,nx,xnode,order)
    end
    for i = Xrng # y direction - left/right
        uₓₓ[i,Yrng] += FirstDerivativeBoundary(u[i+oy,:],Δy,ynode,order)
    end
    return uₓₓ
end


