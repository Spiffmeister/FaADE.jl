

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
    method=:euler,order_x=2,order_y=order_x,α::Float64=1.5,maxIT::Int64=-1,warnings::Bool=false,samplefactor::Int64=1)
    #===== 2D TIME SOLVER =====#

    # Preallocate and set initial
    N = ceil(Int64,t_f/Δt)
    soln = zeros(Float64,ny,nx,ceil(Int64,N/samplefactor))
    uₙ = zeros(Float64,ny,nx)
    uₒ = zeros(Float64,ny,nx)
    for i = 1:nx
        for j = 1:ny
            uₒ[j,i] = u₀(x[i],y[j])
        end
    end

    function storage!(soln,tmp,sample,k)
        if mod(sample,samplefactor) == 0
            soln[:,:,k] = tmp
            k += 1
        end
        return soln, k
    end

    k = 1

    if method != :cgie #Not using CG
        function RHS(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            if boundary_x != :Periodic
                for i = 1:ny #x boundaries
                    SATₗ = SAT_left(boundary_x,u[i,:],Δx,gx(t),c=kx[i,:],order=order_x)
                    SATᵣ = SAT_right(boundary_x,u[i,:],Δx,gx(t),c=kx[i,:],order=order_x)
                    uₓₓ[i,1:order_x] += SATₗ
                    uₓₓ[i,end-order_x+1:end] += SATᵣ
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δx,kx[i,:],order=order_x)
                    uₓₓ[i,1:order_x] += SATₗ
                    uₓₓ[i,end-order_x+1:end] += SATᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx #y boundaries
                    SATₗ = SAT_left(boundary_y,u[:,i],Δy,gy(t),c=ky[:,i],order=order_y)
                    SATᵣ = SAT_right(boundary_y,u[:,i],Δy,gy(t),c=ky[:,i],order=order_y)
                    uₓₓ[1:order_y,i] += SATₗ
                    uₓₓ[end-order_y+1:end,i] += SATᵣ
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δx,ky[:,i],order=order_y)
                    uₓₓ[1:order_y,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            end
            return uₓₓ
        end
    end

    if method == :euler
        for i = 1:N-1
            t = i*Δt
            uₙ = forward_euler(uₙ,uₒ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            soln,k = storage!(soln,uₙ,i,k)
            uₒ = uₙ
        end
    elseif method == :impliciteuler
        for i = 1:N-1
            t = i*Δt
            u[:,:,i+1] = implicit_euler(u[:,:,i+1],u[:,:,i],RHS,nx,ny,Δx,Δy,Δt,kx,ky,t,x,y,gx,gy)
            soln,k=storage!(soln,uₙ,i,k)
            uₒ = uₙ
        end
    elseif method == :cgie
        maxIT == -1 ? maxIT = 15 : nothing

        function cgRHS(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            if boundary_x != :Periodic
                for i = 1:ny
                    SATₗ, = SAT_left(boundary_x,u[i,:],Δx,gx(t),c=kx[i,:],order=order_x,separate_forcing=true)
                    SATᵣ, = SAT_right(boundary_x,u[i,:],Δx,gx(t),c=kx[i,:],order=order_x,separate_forcing=true)
                    uₓₓ[i,1:order_x] += SATₗ
                    uₓₓ[i,end-order_x+1:end] += SATᵣ
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δx,kx[i,:],order=order_x,separate_forcing=true)
                    uₓₓ[i,1:order_x] += SATₗ
                    uₓₓ[i,end-order_x+1:end] += SATᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx
                    SATₗ, = SAT_left(boundary_y,u[:,i],Δy,gy(t),c=ky[:,i],order=order_y,separate_forcing=true)
                    SATᵣ, = SAT_right(boundary_y,u[:,i],Δy,gy(t),c=ky[:,i],order=order_y,separate_forcing=true)
                    uₓₓ[1:order_y,i] += SATₗ
                    uₓₓ[end-order_y+1:end,i] += SATᵣ
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δy,ky[:,i],order=order_y,separate_forcing=true)
                    uₓₓ[1:order_y,i] += SATₗ
                    uₓₓ[end-order_y+1:end,i] += SATᵣ
                end
            end
            return uₓₓ
        end

        Hx = build_H(nx,order_x)*Δx
        Hy = build_H(ny,order_y)*Δy

        for i = 1:N-1
            t = i*Δt
            uⱼ = uₒ
            if boundary_x != :Periodic
                for i = 1:ny
                    SATx,Fₗ = SAT_left(boundary_x,uⱼ[i,:],Δx,gx(t),c=kx[i,:],order=order_x,separate_forcing=true)
                    SATx,Fᵣ = SAT_right(boundary_x,uⱼ[i,:],Δx,gx(t),c=kx[i,:],order=order_x,separate_forcing=true)
                    uⱼ[i,1:order_x] += Δt*Fₗ
                    uⱼ[i,end-order_x+1:end] += Δt*Fᵣ
                end
            end
            if boundary_y != :Periodic
                for i = 1:nx
                    SATy,Fₗ = SAT_left(boundary_y,uⱼ[:,i],Δy,gy(t),c=ky[:,i],order=order_y,separate_forcing=true)
                    SATy,Fᵣ = SAT_right(boundary_y,uⱼ[:,i],Δy,gy(t),c=ky[:,i],order=order_y,separate_forcing=true)
                    uⱼ[1:order_y,i] += Δt*Fₗ
                    uⱼ[end-order_y+1:end,i] += Δt*Fᵣ
                end
            end
            # println("t",i)
            uₙ = conj_grad(uⱼ,uⱼ,cgRHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy,Hx,Hy,tol=1e-5,maxIT=20)
            soln,k = storage!(soln,uₙ,i,k)
            uₒ = uₙ
        end
    end
    return soln
end
