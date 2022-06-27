#= EXPLICIT METHODS =#
function forward_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,g)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,n,x,Δx,t,Δt,k,g)
    return uₙ
end
function forward_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,x,y,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::Matrix,ky::Matrix,gx,gy)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
    return uₙ
end




function RK4(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary)
    k1 = RHS(uₙ,uₒ        ,n,x,Δx,t,Δt,k,       boundary)
    k2 = RHS(uₙ,uₒ+Δt/2*k1,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k3 = RHS(uₙ,uₒ+Δt/2*k2,n,x,Δx,t+0.5Δt,Δt,k, boundary)
    k4 = RHS(uₙ,uₒ+Δt*k3  ,n,x,Δx,t+Δt,Δt,k,    boundary)
    uₙ = uₒ + Δt/6 * (k1 + 2k2 + 2k3 + k4)
    return uₙ
end



#= IMPLICIT METHODS =#
function implicit_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,boundary;maxIT::Int=100,α::Float64=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,n,x,Δx,t,Δt,k,boundary)
    return uⱼ
end
function implicit_euler(uₙ::Matrix,uₒ::Matrix,RHS::Function,nx::Int,ny::Int,Δx::Float64,Δy::Float64,Δt::Float64,kx::Matrix,ky::Matrix,t::Float64,x::Vector,y::Vector,boundary_x,boundary_y;maxIT::Int=100,α::Float64=1.5)
    uⱼ = uₒ
    for j = 1:maxIT
        uⱼ = uⱼ - α * Δt * (uⱼ - uₒ - Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y))
    end
    uⱼ = uₒ + Δt*RHS(uₙ,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,boundary_x,boundary_y)
    return uⱼ
end
    





function conj_grad(b::Vector,uⱼ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,H::Array,boundary;tol::Float64=1e-5,maxIT::Int=10,warnings=false)
    # VECTOR FORM
    xₖ = zeros(length(b)) #Initial guess
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
function conj_grad(b::Matrix,uⱼ::Matrix,RHS::Function,nx::Int,ny::Int,x::Vector,y::Vector,Δx::Float64,Δy::Float64,t::Float64,Δt::Float64,kx::Matrix,ky::Matrix,gx,gy,Hx::Vector{Float64},Hy::Vector{Float64};tol=1e-5,maxIT=10,warnings=false)
    # MATRIX FORM
    xₖ = uⱼ #Initial guess
    rₖ = A(uⱼ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
    while (rnorm > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
        dₖAdₖ = innerH(dₖ,Hx,Hy,Adₖ)
        αₖ = -innerH(rₖ,Hx,Hy,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy) - b
        βₖ = innerH(rₖ,Hx,Hy,A(rₖ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
        i += 1
    end
    if (norm(rₖ)>tol) & warnings
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end





#= SUPPORT =#
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


function A(uⱼ,PDE::Function,n,Δx,x,Δt,t,k,g)
    # tmp can be any vector of length(uⱼ)
    tmp = zeros(length(uⱼ))
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,n,x,Δx,t,Δt,k,g)
    return tmp
end
function A(uⱼ,PDE::Function,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
    # A for 2D arrays
    tmp = zeros(size(uⱼ))
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
    return tmp
end


function innerH(u::Vector,H::Array,v::Vector)
    # H inner product for 1D problems
    return dot(u,H*v)
end
function innerH(u::Matrix,Hx::Vector,Hy::Vector,v::Matrix)
    # H inner product for 2D problems
    nx,ny = size(u)
    tmp = 0.0
    for i = 1:nx
        for j = 1:ny
            tmp += u[i,j]*Hx[i]*Hy[j]*v[i,j]
        end
    end
    return tmp
end
