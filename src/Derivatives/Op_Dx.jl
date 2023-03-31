
"""
    Dₓ

Internally uses [`Dₓ!`](@ref)
"""
function Dₓ end
# 1D First derivative
function Dₓ(u::AbstractVector{T},Δx::T;order::Integer=2) where T
    uₓ = zeros(T,length(u))
    Dₓ!(uₓ,u,length(u),Δx,order)
    return uₓ
end
# 2D First derivative
function Dₓ(u::AbstractMatrix{T},nx::Integer,ny::Integer,Δx::T,Δy::T;order::Integer=2) where T
    
    uₓ = zeros(T,size(u))
    Dₓ!(uₓ,u,nx,ny,Δx,Δy,order)
    
    return uₓ
end


"""
    Dₓ!
In place first derivative operator uses [`FirstDerivativeBoundary!`](@ref) and [`FirstDerivativeInternal!`](@ref)
"""
function Dₓ! end
# 1D
function Dₓ!(uₓ::AbstractVector{T},u::AbstractVector{T},n::Integer,Δx::T,
        order::Integer) where T
    FirstDerivativeBoundary!(uₓ,u,Δx,n,Left,order)
    FirstDerivativeInternal!(uₓ,u,Δx,n,order)
    FirstDerivativeBoundary!(uₓ,u,Δx,n,Right,order)
end
# 1D for 2D
function Dₓ!(uₓ::AbstractArray{T},u::AbstractArray{T},n::Integer,Δ::T,
        dim::Integer,order::Integer) where T
    loopdir = SelectLoopDirection(dim)
    for (cache,U) in zip(loopdir(uₓ),loopdir(u))
        Dₓ!(cache,U,n,Δ,order)
    end
    uₓ
end
# 2D
function Dₓ!(uₓ::AbstractArray{T},u::AbstractArray{T},nx::Integer,ny::Integer,Δx::T,Δy::T,
        order::Integer) where T
    
    rx = 

    FirstDerivativeInternal!(uₓ,u,Δx,Δy,nx,ny,order)

    #Boundary nodes
    FirstDerivativeBoundaryStencil!(uₓ,
        u[1:order],Left,Internal,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ,
        u[1:order],Right,Internal,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ,
        u[1:order],Internal,Left,Δx,Δy,nx,ny,order)
    FirstDerivativeBoundaryStencil!(uₓ,
        u[1:order],Internal,Right,Δx,Δy,nx,ny,order)

    #Corner nodes
    FirstDerivativeCornerStencil!(uₓ,
        u[1:order],Left,Left,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ,
        u[1:order],Left,Right,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ,
        u[1:order],Right,Right,Δx,Δy,nx,ny,order)
    FirstDerivativeCornerStencil!(uₓ,
        u[1:order],Right,Left,Δx,Δy,nx,ny,order)
    
    uₓ
end


"""
    FirstDerivativeBoundaryStencil!
"""
function FirstDerivativeBoundaryStencil! end
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
            uₓₓ[i,j] += FirstDerivativeInternal(u[i:i+order_y,j+offset],Δx,order)
        end
    end
    return uₓₓ
end

"""
    FirstDerivativeCornerStencil!
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


