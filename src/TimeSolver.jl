

struct prob
    RHS             :: Function
    boundary_left   :: Symbol
    boundary_right  :: Symbol
    boundary        :: Union{Vector{Float64},Function}
    order           :: Int64
    k               :: Union{Vector{Float64},Float64}
end

mutable struct solution
    u       :: Array{Float64}
    x       :: Array{Float64}
    Δx      :: Float64
    n       :: Int64
    t_f     :: Float64
    Δt      :: Float64
    solver  :: Symbol

    function solution(u₀,x,Δx,t_f,Δt,solver)
        N = ceil(Int64,t_f/Δt)
        n = length(x)

        u = zeros(Float64,n,N)
        u[:,1] = u₀(x)

        new(u,x,Δx,n,t_f,Δt,solver)
    end
end



"""
    time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2)
"""
function time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;
    boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2,α::Float64=1.5,tol::Float64=1e-5,maxIT::Int64=-1,warnings::Bool=false)

    # Initialise solution
    soln = solution(u₀,x,Δx,t_f,Δt,method)

    # Get the length of the time array
    N = ceil(Int64,t_f/Δt)

    if method != :cgie # Not using conjugate gradient
        function RHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            # Combines the PDE and SATs (forcing term included)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)
            uₓₓ[1:order] .+= SAT_left(boundary_left,u,Δx,g(t),c=k,order=order)
            uₓₓ[end-order+1:end] .+= SAT_right(boundary_right,u,Δx,g(t),c=k,order=order)
            return uₓₓ
        end
    end

    if method == :euler
        # Eulers method
        for i = 1:N-1
            t = i*Δt
            soln.u[:,i+1] = forward_euler(soln.u[:,i+1],soln.u[:,i],RHS,n,Δx,Δt,k,t,x,boundary)
        end
    elseif method == :rk4
        # Runge-Kutta 4
        for i = 1:N-1
            t = i*Δt
            soln.u[:,i+1] = RK4(soln.u[:,i+1],soln.u[:,i],RHS,n,Δx,Δt,k,t,x,boundary)
        end
    elseif method == :impliciteuler
        # Implicit euler
        if maxIT == -1
            maxIT = 100
        end
        for i = 1:N-1
            t = i*Δt
            soln.u[:,i+1] = implicit_euler(soln.u[:,i+1],soln.u[:,i],RHS,n,Δx,Δt,k,t,x,boundary,α=α,maxIT=maxIT)
        end
    elseif method == :cgie
        if maxIT == -1
            maxIT = 15
        end
        function cgRHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)
            SATₗ, = SAT_left(boundary_left,u,Δx,g(t),order=order,seperate_forcing=true)
            SATᵣ, = SAT_right(boundary_left,u,Δx,g(t),order=order,seperate_forcing=true)
            uₓₓ[1:order] += SATₗ
            uₓₓ[end-order+1:end] += SATᵣ
            return uₓₓ
        end

        if order == 2
            H = diagm(ones(length(x)))
            H[1,1] = H[end,end] = 0.5
        elseif order == 4
            H = diagm(ones(length(x)))
            H[1,1] = H[end,end] = 0.5
            H[2,2] = H[end-1,end-1] = 0.5
            H[3,3] = H[end-2,end-2] = 0.5
        elseif order == 6
            H = diagm(ones(length(x)))
        end
        H = Δx*H

        for i = 1:N-1
            t = i*Δt
            uⱼ = soln.u[:,i]
            SATₗ,Fₗ = SAT_left(boundary_left,uⱼ,Δx,boundary(t),order=order,seperate_forcing=true)
            SATᵣ,Fᵣ = SAT_right(boundary_right,uⱼ,Δx,boundary(t),order=order,seperate_forcing=true)
            uⱼ[1:order] += Δt*Fₗ
            uⱼ[end-order+1:end] += Δt*Fᵣ
            soln.u[:,i+1] = conj_grad(uⱼ,uⱼ,cgRHS,n,Δx,Δt,k,t,x,H,boundary;tol=1e-5,maxIT=20)
        end
    else
        error("Method must be :euler, :rk4, :impliciteuler or :cgie")
    end

    return soln
end


# function time_solver(u,PDE,SAT,nx,ny,Δx,Δy,x,y,Δt,k,boundary_x,boundary_y;method=:euler,order=2,maxIT=5)
# end



#= EXPLICIT METHODS =#

function forward_euler(uₙ::Vector,uₒ::Vector,RHS::Function,n::Int,Δx::Float64,Δt::Float64,k::Vector,t::Float64,x::Vector,g)
    # Simple forward euler method
    uₙ = uₒ + Δt*RHS(uₙ,uₒ,n,x,Δx,t,Δt,k,g)
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




#= SUPPORT =#
function A(uⱼ,PDE::Function,n,Δx,x,Δt,t,k,g)
    # tmp can be any vector of length(uⱼ)
    tmp = zeros(length(uⱼ))
    tmp = uⱼ - Δt*PDE(tmp,uⱼ,n,x,Δx,t,Δt,k,g)
    return tmp
end

function A(uⱼ,PDE::Function,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y)
    # A for 2D arrays
    tmp = uⱼ
    tmp -= Δt*PDE(tmp,uⱼ,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y)
    return tmp
end


function innerH(u::Vector,H::Array,v::Vector)
    # H inner product for 1D problems
    return dot(u,H*v)
end

function innerH(u::Matrix,Hx::Vector,Hy::Vector,v::Matrix)
    # H inner product for 2D problems
    nx,ny = size(u)
    tmp = 0
    for i = 1:nx
        for j = 1:ny
            tmp += u[i,j]*v[i,j]*Hx[i]*Hy[j]
        end
    end
    return tmp
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

function conj_grad(b,uⱼ,RHS,nx,ny,Δx,x,Δy,y,Δt,t,k,Hx,Hy;order=2,tol=1e-5,maxIT=10)
    xₖ = zeros(length(b)) #Initial guess
    rₖ = A(uⱼ,RHS,n,Δx,Δt,k,t,x,order=order) - b
    dₖ = -rₖ
    i = 0
    rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
    while (rnorm > tol) & (i < maxIT)
        Adₖ = A(dₖ,RHS,nx,ny,Δx,x,Δy,y,Δt,t,k,boundary_x,boundary_y,order=order)
        dₖAdₖ = innerH(dₖ,H,Adₖ)
        αₖ = -innerH(rₖ,H,dₖ)/dₖAdₖ
        xₖ = xₖ + αₖ*dₖ
        rₖ = A(xₖ,RHS,n,Δx,Δt,k,t,x,order=order) - b
        βₖ = innerH(rₖ,H,A(rₖ,RHS,n,Δx,Δt,k,t,x,order=order))/dₖAdₖ
        dₖ = - rₖ + βₖ*dₖ
        rnorm = sqrt(innerH(rₖ,Hx,Hy,rₖ))
        i += 1
    end
    if (norm(rₖ)>tol)
        warnstr = string("CG did not converge at t=",t)
        @warn warnstr
    end
    return xₖ
end


