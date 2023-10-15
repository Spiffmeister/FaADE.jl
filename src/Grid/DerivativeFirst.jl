#=====================================#
#====== FIRST DERIVATIVE METHODS =====#
#=====================================#
# Author: Dean Muir, Kenneth Duru

"""
    D₁!
1D and 2D in place first derivative operator.

See also [`FirstDerivativeBoundary!`](@ref) and [`FirstDerivativeInternal!`](@ref).
"""
function D₁! end
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
    FirstDerivative
Allocating functions for first derivative, useful when need to add value to a matrix.

See also [`FirstDerivativeInternal!`](@ref) and [`FirstDerivativeBoundary`](@ref)
"""
function FirstDerivative end
### Internal nodes
"""
    FirstDerivative(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer)
Internal node function.
"""
function FirstDerivative(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer) where T
    uₓ = zeros(T,length(u))
    FirstDerivativeBoundary!(uₓ,u,Δx,length(u),node,order)
    return uₓ
end
### Boundary nodes
"""
    FirstDerivative(u::AbstractVector{T},Δx::T,node::NodeType,order::Integer)
Boundary node function.
"""
function FirstDerivative(u::AbstractVector{T},Δx::T,node::NodeType,order::Integer) where T
    order == 2 ? nnodes = 1 : nnodes = order
    uₓ = zeros(T,nnodes)
    FirstDerivativeBoundary!(uₓ,u,Δx,length(u),node,order)
    return uₓ
end

"""
    FirstDerivativeInternal(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer)
Single node 1D first derivative function.

See also [`FirstDerivativeInternal!`](@ref) for in place multi-node functions
"""
@inline function FirstDerivativeInternal(u::AbstractVector{T},Δx::T,::NodeType{:Internal},order::Integer) where T
    if order == 2
        j = 2
        return (u[j+1] - u[j-1])/(T(2)*Δx)
    elseif order == 4
        j = 3
        return (T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
    elseif order == 6
        j = 4
        return (-T(1/60)*u[i-3] + T(3/20)*u[j-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
    end
end

"""
    FirstDerivativeInternal!
1D and 2D in place functions for first derivative on internal nodes of a domain
"""
function FirstDerivativeInternal!() end
######### 1D FUNCTION
@views @inline function FirstDerivativeInternal!(uₓ::AbstractVector{T},u::AbstractVector{T},Δx::T,n::Integer,order::Integer) where T
    if order == 2
        for i = 2:n-1
            uₓ[i] = (u[j+1] - u[j-1])/(T(2)*Δx)
        end
    elseif order == 4
        for i = 3:n-2
            uₓ[i] = (T(1/12)*u[i-2] - T(2/3)*u[i-1] + T(2/3)*u[i+1] - T(1/12)*u[i+2])/Δx
        end
    elseif order == 6
        for i = 4:n-3
            uₓ[i] = (-T(1/60)*u[i-3] + T(3/20)*u[j-2] - T(3/4)*u[i-1] + T(3/4)*u[i+1] - T(3/20)*u[i+2] + T(1/60)*u[i+3])/Δx
        end
    end
end
######### 2D FUNCTION
@views @inline function FirstDerivativeInternal!(uₓ::AbstractArray{T},
        u::AbstractArray{T},
        Δx::T,Δy::T,nx::Integer,ny::Integer,order::Integer) where T

    if order == 2
        # SECOND ORDER METHOD
        for j=2:ny-1
            for i = 2:nx-1
                uₓ[i,j] = (u[i+1,j] - u[i-1,j])/(T(2)*Δx) + (u[i,j+1] - u[i,j-1])/(T(2)*Δy)
            end
        end
    elseif order == 4
        # FOURTH ORDER METHOD
        for j = 3:ny-2
            for i = 3:nx-2
                uₓ[i,j] = (T(1/12)*u[i-2,j] - T(2/3)*u[i-1,j] + T(2/3)*u[i+1,j] - T(1/12)*u[i+2,j])/Δx + 
                (T(1/12)*u[i,j-2] - T(2/3)*u[i,j-1] + T(2/3)*u[i,j+1] - T(1/12)*u[i,j+2])/Δy
            end
        end
    elseif order == 6
        # SIXTH ORDER METHOD
        for j = 4:ny-3
            for i = 4:nx-3
                uₓ[i,j] = (-T(1/60)*u[i-3,j] + T(3/20)*u[i-2,j] - T(3/4)*u[i-1,j] + T(3/4)*u[i+1,j] - T(3/20)*u[i+2,j] + T(1/60)*u[i+3,j])/Δx +
                            (-T(1/60)*u[i,j-3] + T(3/20)*u[i,j-2] - T(3/4)*u[i,j-1] + T(3/4)*u[i,j+1] - T(3/20)*u[i,j+2] + T(1/60)*u[i,j+3])/Δy
            end
        end
    end
end


"""
    FirstDerivativeBoundary!
1D in place function for first derivative on boundary nodes
"""
@views function FirstDerivativeBoundary!(uₓ::AbstractVector{T},
        u::AbstractVector{T},Δx::T,n::Integer,NT::NodeType,order::Integer) where T
    NT <: Left ? i = 1 : i = -1
    NT <: Left ? j = 1 : j = n
    if order == 2
        uₓ[j] = (u[j+i] - u[j])/Δx
        uₓ[j] = uₓ[j]/Δx
    elseif order == 4
        uₓ[j]       = T(i)*(T(24/17)*u[j] + T(59/34)*u[j+i]   - T(4/17)*u[j+2i]   - T(3/34)*u[j+3i])/Δx
        uₓ[j+i]     = T(i)*(T(1/2)*u[j]   + T(1/2)*u[j+2i])/Δx
        uₓ[j+2i]    = T(i)*(T(4/43)*u[j]  - T(59/86)*u[j+i]   + T(59/86)*u[j+3i]  - T(4/43)*u[j+4i])/Δx
        uₓ[j+3i]    = T(i)*(T(3/98)*u[j]  - T(59/98)*u[j+2i]  + T(32/49)*u[j+4i]  - T(4/49)*u[j+5i])/Δx
        uₓ[j:i:j+3i] = uₓ[j:i:j+5i]/Δx
    elseif order == 6
        uₓ[j]       = T(i)*( T(-1.582533518939116)*u[j] + T(2.033378678700676)*u[j+i] - T(0.141512858744873)*u[j+2i] + T(-0.450398306578272)*u[j+3i] + T(0.104488069284042)*u[j+4i] + T(0.036577936277544)*u[j+5i] )
        uₓ[j+i]     = T(i)*( T(-0.462059195631158)*u[j] + T(0.287258622978251)*u[j+2i] + T(0.258816087376832)*u[j+3i] + T(-0.069112065532624)*u[j+4i] - T(0.014903449191300)*u[j+5i] )
        uₓ[j+2i]    = T(i)*( T(0.071247104721830)* u[j] - T(0.636451095137907)*u[j+i] + T(0.606235523609147)*u[j+3i] + T(-0.022902190275815)*u[j+4i] - T(0.018129342917256)*u[j+5i] )
        uₓ[j+3i]    = T(i)*( T(0.114713313798970)* u[j] - T(0.290087484386815)*u[j+i] - T(0.306681191361148)*u[j+2i] + T(0.520262285050482)*u[j+4i]  - T(0.051642265516119)*u[j+5i] + T(0.013435342414630)*u[j+6i] )
        uₓ[j+4i]    = T(i)*( T(-0.036210680656541)*u[j] + T(0.105400944933782)*u[j+i] + T(0.015764336127392)*u[j+2i] + T(-0.707905442575989)*u[j+3i] + T(0.769199413962647)*u[j+5i] - T(0.164529643265203)*u[j+6i] + T(0.018281071473911)*u[j+7i] )
        uₓ[j+5i]    = T(i)*( T(-0.011398193015050)*u[j] + T(0.020437334208704)*u[j+i] + T(0.011220896474665)*u[j+2i] + T( 0.063183694641876)*u[j+3i] - T(0.691649024426814)*u[j+4i] + T(0.739709139060752)*u[j+6i] + T(-0.147941827812150)*u[j+7i] + T(0.016437980868017)*u[j+8i] )
        uₓ[j:i:j+5i] = uₓ[j:i:j+5i]/Δx
    end
end
