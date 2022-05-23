

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

struct grid
    x   :: Array{Float64}
    Δx  :: Float64
    n   :: Int64

    function grid(𝒟,n)
        Δx = (𝒟[2] - 𝒟[1])/(n-1)
        x = collect(range(𝒟[1],𝒟[2],step=Δx))
        new(x,Δx,n)
    end
end

"""
    time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2)
"""
function time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;
    boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2,α::Float64=1.5,tol::Float64=1e-5,maxIT::Int64=-1,warnings::Bool=false)
    #===== 1D TIME SOLVER =====#

    # Initialise solution
    soln = solution(u₀,x,Δx,t_f,Δt,method)
    # Get the length of the time array
    N = ceil(Int64,t_f/Δt)

    soln.u[1,1] = boundary(0)[1]
    soln.u[end,1] = boundary(0)[end]

    if method != :cgie # Not using conjugate gradient
        function RHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            # Combines the PDE and SATs (forcing term included)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)

            if boundary_left != :Periodic
                uₓₓ[1:order] .+= SAT_left(boundary_left,u,Δx,g(t),c=k,order=order)
                uₓₓ[end-order+1:end] .+= SAT_right(boundary_right,u,Δx,g(t),c=k,order=order)
            else
                SATₗ,SATᵣ = SAT_Periodic(u,Δx,k,order=order)
                uₓₓ[1:order] += SATₗ
                uₓₓ[end-order+1:end] += SATᵣ
            end
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
            if boundary_left != :Periodic
                SATₗ, = SAT_left(boundary_left,u,Δx,g(t),order=order,separate_forcing=true)
                SATᵣ, = SAT_right(boundary_left,u,Δx,g(t),order=order,separate_forcing=true)
                uₓₓ[1:order] += SATₗ
                uₓₓ[end-order+1:end] += SATᵣ
            else
                SATₗ,SATᵣ, = SAT_Periodic(u,Δx,k,order=order,separate_forcing=true)
                uₓₓ[1:order] += SATₗ
                uₓₓ[end-order+1:end] += SATᵣ
            end
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
            if boundary_left != :Periodic
                SATₗ,Fₗ = SAT_left(boundary_left,uⱼ,Δx,boundary(t),order=order,separate_forcing=true)
                SATᵣ,Fᵣ = SAT_right(boundary_right,uⱼ,Δx,boundary(t),order=order,separate_forcing=true)
            else
                SATₗ,SATᵣ,Fₗ,Fᵣ = SAT_Periodic(uⱼ,Δx,k,order=order,separate_forcing=true)
            end
            uⱼ[1:order] += Δt*Fₗ
            uⱼ[end-order+1:end] += Δt*Fᵣ
            soln.u[:,i+1] = conj_grad(uⱼ,uⱼ,cgRHS,n,Δx,Δt,k,t,x,H,boundary;tol=1e-5,maxIT=20)
        end
    else
        error("Method must be :euler, :rk4, :impliciteuler or :cgie")
    end

    return soln
end

function time_solver(PDE::Function,u₀::Function,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,x::Vector{Float64},y::Vector{Float64},t_f::Float64,Δt::Float64,kx::Matrix{Float64},ky::Matrix{Float64},gx,gy,boundary_x::Symbol,boundary_y::Symbol;
    method=:euler,order_x=2,order_y=order_x,α::Float64=1.5,maxIT::Int64=-1,warnings::Bool=false)
    #===== 2D TIME SOLVER =====#

    # Preallocate and set initial
    N = ceil(Int64,t_f/Δt)
    u = zeros(Float64,nx,ny,N)
    for i = 1:nx
        for j = 1:ny
            u[i,j,1] = u₀(x[i],y[j])
        end
    end

    if method != :cgie #Not using CG
        function RHS(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            if boundary_x != :Periodic
                for i = 1:ny #x boundaries
                    SATₗ = SAT_left(boundary_x,u[:,i],Δx,gx(t),c=kx[i,:],order=order_x)
                    SATᵣ = SAT_right(boundary_x,u[:,i],Δx,gx(t),c=kx[i,:],order=order_x)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δx,kx[:,i],order=order_x)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx #y boundaries
                    SATₗ = SAT_left(boundary_y,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                    SATᵣ = SAT_right(boundary_y,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δx,ky[i,:],order=order_y)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_x+1:end] += SATᵣ
                end
            end
            return uₓₓ
        end
    end
    
    if method == :euler
        for i = 1:N-1
            t = i*Δt
            u[:,:,i+1] = forward_euler(u[:,:,i+1],u[:,:,i],RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
        end
    elseif method == :impliciteuler
        for i = 1:N-1
            t = i*Δt
            u[:,:,i+1] = implicit_euler(u[:,:,i+1],u[:,:,i],RHS,nx,ny,Δx,Δy,Δt,kx,ky,t,x,y,gx,gy)
        end
    elseif method == :cgie
        maxIT == -1 ? maxIT = 15 : nothing

        function cgRHS(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            if boundary_x != :Periodic
                for i = 1:ny
                    SATₗ, = SAT_left(boundary_x,u[:,i],Δx,gx(t),c=kx[:,i],order=order_x,separate_forcing=true)
                    SATᵣ, = SAT_right(boundary_x,u[:,i],Δx,gx(t),c=kx[:,i],order=order_x,separate_forcing=true)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δx,kx[:,i],order=order_x,separate_forcing=true)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx
                    SATₗ, = SAT_left(boundary_y,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y,separate_forcing=true)
                    SATᵣ, = SAT_right(boundary_y,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y,separate_forcing=true)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δy,ky[i,:],order=order_y,separate_forcing=true)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            end
            return uₓₓ
        end

        Hx = build_H(nx,order_x)*Δx
        Hy = build_H(ny,order_y)*Δy

        for i = 1:N-1
            t = i*Δt
            uⱼ = u[:,:,i]
            if boundary_x != :Periodic
                for i = 1:ny
                    SATx,Fₗ = SAT_left(boundary_x,uⱼ[:,i],Δx,gx(t),c=kx[:,i],order=order_x,separate_forcing=true)
                    SATx,Fᵣ = SAT_right(boundary_x,uⱼ[:,i],Δx,gx(t),c=kx[:,i],order=order_x,separate_forcing=true)
                    uⱼ[1:order_x,i] += Δt*Fₗ
                    uⱼ[end-order_x+1:end,i] += Δt*Fᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx
                    SATy,Fₗ = SAT_left(boundary_y,uⱼ[i,:],Δy,gy(t),c=ky[i,:],order=order_y,separate_forcing=true)
                    SATy,Fᵣ = SAT_right(boundary_y,uⱼ[i,:],Δy,gy(t),c=ky[i,:],order=order_y,separate_forcing=true)
                    uⱼ[i,1:order_y] += Δt*Fₗ
                    uⱼ[i,end-order_y+1:end] += Δt*Fᵣ
                end
            end
            # println("t",i)
            u[:,:,i+1] = conj_grad(uⱼ,uⱼ,cgRHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy,Hx,Hy,tol=1e-5,maxIT=20)
        end
    end
    return u
end



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
    xₖ = zeros(size(b)) #Initial guess
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


