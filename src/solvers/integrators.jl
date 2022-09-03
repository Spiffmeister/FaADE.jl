#= EXPLICIT METHODS =#

"""
    forward_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,g)
or
    function forward_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,x,y,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::Matrix,ky::Matrix,gx,gy)

Inbuilt forward euler method,
``u^{n+1} = u^n + \\Delta t F(u^{n+1},u^n,x,\\Delta x, t,\\Delta t,k,g)``
- `n` (`nx`, `ny`): is the number of grid points
- `Δx` (`Δy`): is the grid spacing
- `k` (`kx`, `gy`): is a matrix or function containing the diffusion coefficients at the grid points
- `g` (`gx`, `gy`): is a vector or function containing the boundary conditions

"""
function forward_euler end
function forward_euler(uₙ::AbstractVector,uₒ::AbstractVector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::AbstractVector,t::Float64,x::AbstractVector,g)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,n,x,Δx,t,Δt,k,g)
    return uₙ
end
# Matrix version
function forward_euler(uₙ::AbstractMatrix,uₒ::AbstractMatrix,RHS::Function,nx::Int,ny::Int,x,y,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::AbstractMatrix,ky::AbstractMatrix,gx,gy)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
    return uₙ
end


"""
    RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
"""
function RK4 end
function RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
    k1 = RHS(uₙ,uₒ        ,n,x,Δx,t,Δt,k,       boundary)
    k2 = RHS(uₙ,uₒ+Δt/2*k1,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k3 = RHS(uₙ,uₒ+Δt/2*k2,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k4 = RHS(uₙ,uₒ+Δt*k3  ,n,x,Δx,t+Δt,Δt,k,    boundary)
    uₙ = uₒ + Δt/6 * (k1 + 2k2 + 2k3 + k4)
    return uₙ
end



#= IMPLICIT METHODS =#
"""
    implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
or
    function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
"""
function implicit_euler end
# Vector form
function implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary)
    return uⱼ
end
# Matrix form
function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y)
    return uⱼ
end
    


"""
    conj_grad(b::Vector,uⱼ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,H::Array,boundary;tol::Float64=1e-5,maxIT::Int=10,warnings=false)
or
    conj_grad(b::AbstractMatrix,uⱼ::AbstractMatrix,RHS::Function,nx::Int,ny::Int,x::Vector,y::Vector,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::Matrix,ky::Matrix,gx,gy,Hx::Vector{Float64},Hy::Vector{Float64}
    ;tol=1e-5,rtol=1e-10,maxIT=10,warnings=true)

See also [`build_H`](@ref), [`A`](@ref)
"""
function conj_grad end
# VECTOR FORM
function conj_grad(b::Vector,uⱼ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,H::Array,boundary;tol::Float64=1e-5,maxIT::Int=10,warnings=false)
    xₖ = uⱼ #Initial guess
    rₖ = A(uⱼ,RHS,n,Δx,x,Δt,t,k,boundary) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,H,rₖ))
    while (rnorm > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,n,Δx,x,Δt,t,k,boundary)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,x,Δt,t,k,boundary) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,x,Δt,t,k,boundary))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,H,rₖ))
        # rnorm = norm(rₖ)
        i += 1
    end
    if (norm(rₖ)>tol) & warnings
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end
function conj_grad(b::AbstractMatrix,uⱼ::AbstractMatrix,RHS::Function,nx::Int,ny::Int,x::Vector,y::Vector,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::Matrix,ky::Matrix,Hx::Vector{Float64},Hy::Vector{Float64}
    ;tol=1e-5,rtol=1e-10,maxIT=10,warnings=true)
    # MATRIX FORM
    converged = true
    xₖ = uⱼ #Initial guess
    rₖ = A(uⱼ,RHS,x,y,Δx,Δy,Δt,kx,ky) - b
    dₖ = -rₖ
    i = 0

    Adₖ = zeros(Float64,(nx,ny))

    rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ,nx,ny))
    unorm = sqrt(innerH(uⱼ,Hx,Hy,uⱼ,nx,ny))
    while (rnorm > rtol*unorm) & (i < maxIT)
        Adₖ = A!(Adₖ,dₖ,RHS,x,y,Δx,Δy,Δt,kx,ky)
        dₖAdₖ = innerH(dₖ,Hx,Hy,Adₖ,nx,ny)
        αₖ = -innerH(rₖ,Hx,Hy,dₖ,nx,ny)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A!(rₖ,xₖ,RHS,x,y,Δx,Δy,Δt,kx,ky) .- b
        βₖ = innerH(rₖ,Hx,Hy,A(rₖ,RHS,x,y,Δx,Δy,Δt,kx,ky),nx,ny)/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ,nx,ny))
        i += 1
    end
    # println(i," ",rnorm/unorm,"   ",rtol,"  ",t,"   ",Δt)
    if (rnorm>rtol*unorm) & warnings
        converged = false
        warnstr = string("CG did not converge at t=",t,"    Δt=",Δt,"   # iters=",i)
        @warn warnstr
        # xₖ = uⱼ #return to previous time
        # t = t-Δt
    end
    # if adaptive
    #     Δt = adapt_time(converged,Δt,5.0*min(Δx^2,Δy^2))
    # end
    # println("Δt=",Δt)
    return xₖ, converged
end





#= SUPPORT =#
"""
    build_H(n::Int64,order::Int64)
"""
function build_H(n::Int64,order::Int64)
    H = ones(n)
    if order == 2
        H[1] = H[end] = 0.5
    elseif order == 4
        H[1] = H[end] = 17.0/48.0
        H[2] = H[end-1] = 59.0/48.0
        H[3] = H[end-2] = 43.0/48.0
        H[4] = H[end-3] = 49.0/48.0
    elseif order == 6
        H[1] = H[end] = 13649.0/43200.0
        H[2] = H[end-1] = 12013.0/8640.0
        H[3] = H[end-2] = 2711.0/4320.0
        H[4] = H[end-3] = 5359.0/4320.0
        H[5] = H[end-4] = 7877.0/8640.0
        H[6] = H[end-5] = 43801.0/43200.0
    else
        error("Order must be 2,4 or 6")
    end
    return H
end

"""
    A

    1. A(uⱼ::AbstractVector,PDE::Function,n::Int64,Δx::Float64,x::Vector,Δt::Float64,t::Float64,k::Vector{Float64},g)
    2. A(uⱼ::AbstractMatrix,PDE::Function,x,y,Δx,Δy,t,Δt,kx,ky)
"""
function A end
function A(uⱼ::AbstractVector,PDE::Function,n::Int64,Δx::Float64,x::Vector,Δt::Float64,t::Float64,k::Vector,g)
    # tmp can be any vector of length(uⱼ)
    tmp = zeros(Float64,length(uⱼ))
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,n,x,Δx,t,Δt,k,g)
    return tmp
end
function A(uⱼ::AbstractMatrix,PDE::Function,x,y,Δx,Δy,Δt,kx,ky)
    # A for 2D arrays
    tmp = zeros(Float64,size(uⱼ))
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,kx,ky)
    return tmp
end
function A!(tmp,uⱼ::AbstractMatrix,PDE::Function,x,y,Δx,Δy,Δt,kx,ky)
    # A for 2D arrays
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,kx,ky)
    return tmp
end


"""
    innerH

    1. innerH(u::AbstractVector,H::AbstractArray,v::AbstractVector)
    2. innerH(u::AbstractMatrix,Hx::AbstractVector,Hy::AbstractVector,v::AbstractMatrix)
"""
function innerH end
# H inner product for 1D problems
function innerH(u::AbstractVector,H::AbstractArray,v::AbstractVector)
    return dot(u,H*v)
end
# H inner product for 2D problems
function innerH(u::AbstractMatrix,Hx::AbstractVector,Hy::AbstractVector,v::AbstractMatrix)
    nx,ny = size(u)
    tmp = 0.0
    for j = 1:ny
        for i = 1:nx
            tmp += u[i,j]*Hx[i]*Hy[j]*v[i,j]
        end
    end
    return tmp
end
function innerH(u::AbstractMatrix,Hx::AbstractVector,Hy::AbstractVector,v::AbstractMatrix,nx,ny)
    tmp = 0.0
    for j = 1:ny
        for i = 1:nx
            tmp += u[i,j]*Hx[i]*Hy[j]*v[i,j]
        end
    end
    return tmp
end