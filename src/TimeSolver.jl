





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



mutable struct solution
    u       :: Vector{Vector{Float64}}
    x       :: Vector{Float64}
    Δt      :: Union{Float64,Vector{Float64}}
    t       :: Vector{Float64}

    function solution(u₀,x,Δx,t,Δt;preallocate=false)
        if preallocate
            N = ceil(Int64,t/Δt)
            n = length(x)
            u = [zeros(Float64,n) for _ in 1:N]

            u[1] = u₀

            new(u,x,Δt,collect(range(0.0,t,length=N)))
        else
            u = u₀
            new([u],x,[Δt],[t])
        end

    end
end

mutable struct solution_2d
    u   :: Vector{Matrix{Float64}}
    x   :: Vector{Float64}
    y   :: Vector{Float64}
    Δt  :: Vector{Float64}
    t   :: Vector{Float64}
    function solution_2d(u₀,x,y,t,Δt)
        new([u₀],x,y,[Δt],[t])
    end
end


"""
    Struct for storing checkpoints for 2D simulations
"""
mutable struct checkpoint_2d
    # Solution info
    soln        :: solution_2d
    Δx          :: Float64
    Δy          :: Float64
    t_f         :: Float64
    # Diffusion coefficient matricies
    kx          :: Matrix{Float64}
    ky          :: Matrix{Float64}
    # Boundary functions
    gx          :: Function
    gy          :: Function
    # Parallel penalty function if provided
    parallel    :: Bool
    penalty_fn  :: Union{Function,Nothing}
    # Simulation parameters
    order_x     :: Int64
    order_y     :: Int64
    method      :: Symbol
    maxIT       :: Int64
    samplefactor:: Float64
    tol         :: Float64
    rtol        :: Float64
    adaptive    :: Bool
end




"""
    time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::Symbol;
        boundary_right::Symbol=boundary_left,method::Symbol=:euler,order::Int64=2)
or
    function time_solver(PDE::Function,u₀::Function,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,x::Vector{Float64},y::Vector{Float64},t_f::Float64,Δt::Float64,kx::Matrix{Float64},ky::Matrix{Float64},gx,gy,boundary_x::Symbol,boundary_y::Symbol;
        method=:euler,order_x=2,order_y=order_x,maxIT::Int64=15,warnings::Bool=false,samplefactor::Int64=1,tol=1e-5,adaptive=true,penalty_fn=nothing)

Inbuilt function for integrating in time.

See [`forward_euler`](@ref), [`RK4`](@ref), [`implicit_euler`](@ref), [`conj_grad`](@ref)
"""
function time_solver end
#===== 1D TIME SOLVER =====#
function time_solver(PDE::Function,u₀::Function,n::Int64,x::Vector{Float64},Δx::Float64,t_f::Float64,Δt::Float64,k::Vector{Float64},boundary::Function,boundary_left::BoundaryCondition;
    boundary_right::BoundaryCondition=boundary_left,method::Symbol=:euler,order::Int64=2,α::Float64=1.5,tol::Float64=1e-5,maxIT::Int64=-1,warnings::Bool=false,samplefactor=1.0)

    uₙ = zeros(Float64,n)

    uₒ = u₀.(x)

    # Initialise solution
    soln = solution(uₒ,x,Δx,0.0,Δt)
    # Get the length of the time array
    N = ceil(Int64,t_f/Δt)
    t = Δt
    Δt₀ = Δt

    if method != :cgie # Not using conjugate gradient
        soln = solution(uₒ,x,Δx,t_f,Δt,preallocate=true)

        function RHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            # Combines the PDE and SATs (forcing term included)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)

            if boundary_left != Periodic
                uₓₓ[1:order]    += SAT(boundary_left,Left,u,Δx,g(t),c=k,order=order)
                uₓₓ[n-order+1:n]+= SAT(boundary_right,Right,u,Δx,g(t),c=k,order=order)
            else
                SATₗ,SATᵣ = SAT_Periodic(u,Δx,k,order=order)
                uₓₓ[1:order]        += SATₗ
                uₓₓ[end-order+1:end]+= SATᵣ
            end
            return uₓₓ
        end
    else
        soln = solution(uₒ,x,Δx,0.0,Δt)
    end

    if method == :euler
        # Eulers method
        for i = 1:N-1
            t = i*Δt
            soln.u[i+1] = forward_euler(soln.u[i+1],soln.u[i],RHS,n,Δx,Δt,k,t,x,boundary)
        end
    elseif method == :rk4
        # Runge-Kutta 4
        for i = 1:N-1
            t = i*Δt
            soln.u[i+1] = RK4(soln.u[i+1],soln.u[i],RHS,n,Δx,Δt,k,t,x,boundary)
        end
    elseif method == :impliciteuler
        # Implicit euler
        if maxIT == -1
            maxIT = 100
        end
        for i = 1:N-1
            t = i*Δt
            soln.u[i+1] = implicit_euler(soln.u[i+1],soln.u[i],RHS,n,Δx,Δt,k,t,x,boundary,α=α,maxIT=maxIT)
        end
    elseif method == :cgie
        if maxIT == -1
            maxIT = 15
        end
        function cgRHS(uₓₓ,u,n,x,Δx,t,Δt,k,g)
            uₓₓ = PDE(uₓₓ,u,n,x,Δx,t,Δt,k,order=order)
            if boundary_left != :Periodic
                # uₓₓ[1:order]        += SAT(boundary_left,Left,u,Δx,order=order,forcing=false)
                uₓₓ[1:order]        .= SAT_Dirichlet!(uₓₓ[1:order],Left,u,Δx,order=order,forcing=false)
                # uₓₓ[end-order+1:end]+= SAT(boundary_right,Right,u,Δx,order=order,forcing=false)
                println(uₓₓ[end-order+1:end])
                uₓₓ[end-order+1:end].= SAT_Dirichlet!(uₓₓ[end-order+1:end],Right,u,Δx,order=order,forcing=false)
                println(uₓₓ[end-order+1:end])
                println()
            else
                SATₗ,SATᵣ, = SAT_Periodic(u,Δx,k,order=order)
                uₓₓ[1:order]        += SATₗ
                uₓₓ[end-order+1:end]+= SATᵣ
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

        while t ≤ t_f
            uⱼ = copy(uₒ)
            if boundary_left != :Periodic
                # Fₗ = SAT(boundary_left,Left,boundary(t),Δx,order=order,forcing=true)

                Fᵣ = SAT(boundary_right,Right,boundary(t),Δx,order=order,forcing=true)
                println("yes",Fᵣ)

                uⱼ[1:order] .= Δt*SAT_Dirichlet!(uⱼ[1:order],Left,boundary(t),Δx,order=order,forcing=true)
                println("yes",uⱼ[n-order+1:n])
                # uⱼ[n-order+1:n] .= Δt*SAT_Dirichlet!(uⱼ[n-order+1:n],Right,boundary(t),Δx,order=order,forcing=true)
                # println("yes",uⱼ[n-order+1:n])
                # uⱼ[1:order] += Δt*Fₗ
                uⱼ[end-order+1:end] += Δt*Fᵣ
                println("yes",uⱼ[end-order+1:end])
                println()

            end
            uₙ = conj_grad(uⱼ,uⱼ,cgRHS,n,Δx,Δt,k,t,x,H,boundary;tol=1e-5,maxIT=20)

            # if (samplefactor - t) ≤ 0.0
            #     samplefactor += samplefactor
            # end
            t += Δt
            append!(soln.t,t)
            append!(soln.Δt,Δt)
            push!(soln.u,copy(uₙ))
            uₒ = uₙ
        end
    else
        error("Method must be :euler, :rk4, :impliciteuler or :cgie")
    end

    return soln
end
#===== 2D TIME SOLVER =====#
function time_solver(PDE::Function,u₀::Function,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,x::AbstractVector{Float64},y::AbstractVector{Float64},t_f::Float64,Δt::Float64,kx::AbstractMatrix{Float64},ky::AbstractMatrix{Float64},gx,gy,boundary_x::BoundaryCondition,boundary_y::BoundaryCondition;
    method=:euler,order_x=2,order_y=order_x,maxIT::Int64=15,warnings::Bool=false,samplefactor::Float64=0.0,tol=1e-5,rtol=1e-14,adaptive=true,penalty_fn=nothing,checkpointing=false)

    # Preallocate and set initial
    N = ceil(Int64,t_f/Δt)
    # soln = zeros(Float64,nx,ny,ceil(Int64,N))#/samplefactor))
    uₙ = zeros(Float64,nx,ny)
    uₒ = zeros(Float64,nx,ny)
    if nprocs() > 1
        uₙ = SharedArray(uₙ)
        uₒ = SharedArray(uₒ)
    end

    for i = 1:nx
        for j = 1:ny
            uₒ[i,j] = u₀(x[i],y[j])
        end
    end
    soln = solution_2d(copy(uₒ),x,y,0.0,Δt)

    umw = [0.0]


    parallel_penalty = false
    if typeof(penalty_fn) <: Function
        parallel_penalty = true
    end


    function storage(soln,tmp,samplefactor,t)
        if (samplefactor - t) ≤ 0.0
            push!(soln,tmp)
            samplefactor += samplefactor
        end
        return soln, samplefactor
    end

    t = Δt
    Δt₀ = Δt
    samplefactor_base = samplefactor


    if method != :cgie #Not using CG
        function RHS(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            if boundary_x != Periodic
                for i = 1:ny #x boundaries
                    uₓₓ[1:order_x,i]            += SAT(boundary_x,Left,u[:,i],Δx,gx(t),c=kx[:,i],order=order_x)
                    uₓₓ[end-order_x+1:end,i]    += SAT(boundary_x,Right,u[:,i],Δx,gx(t),c=kx[:,i],order=order_x)
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δx,kx[:,i],order=order_x)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            end
            if boundary_y != Periodic
                for i = 1:nx #y boundaries
                    uₓₓ[i,1:order_y]            += SAT(boundary_y,Left,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                    uₓₓ[i,end-order_y+1:end]    += SAT(boundary_y,Right,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δy,ky[i,:],order=order_y)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            end
            return uₓₓ
        end
    end
    if method == :euler
        for i = 1:N-1
            t = i*Δt
            uₙ = forward_euler(uₙ,uₒ,RHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy)
            # soln.u, samplefactor = storage(soln.u,uₙ,samplefactor,t)
            push!(soln.u,uₙ)
            uₒ = uₙ
        end
    elseif method == :impliciteuler
        for i = 1:N-1
            t = i*Δt
            u[:,:,i+1] = implicit_euler(u[:,:,i+1],u[:,:,i],RHS,nx,ny,Δx,Δy,Δt,kx,ky,t,x,y,gx,gy)
            soln = storage!(soln,uₙ,i)
            uₒ = uₙ
        end
    elseif method == :cgie
        maxIT < 1 ? maxIT = 10 : nothing

        function cgRHS(uₓₓ,u)
            uₓₓ = PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            ### SATs
            if boundary_x != Periodic
                for i = 1:ny
                    # uₓₓ[1:order_x,i]        += SAT(boundary_x,Left,u[:,i],Δx,c=kx[:,i],order=order_x,forcing=false)
                    uₓₓ[1:order_x,i]        .= SAT_Dirichlet!(uₓₓ[1:order_x,i],Left,u[1:order_x,i],Δx,c=kx[1:order_x,i],order=order_x,forcing=false)
                    # uₓₓ[end-order_x+1:end,i]+= SAT(boundary_x,Right,u[:,i],Δx,c=kx[:,i],order=order_x,forcing=false)
                    uₓₓ[nx-order_x+1:nx,i]  .= SAT_Dirichlet!(uₓₓ[nx-order_x+1:nx,i],Right,u[:,i],Δx,c=kx[nx-order_x+1:nx,i],order=order_x,forcing=false)
                end
            else
                for i = 1:ny
                    SATₗ,SATᵣ = SAT_Periodic(u[:,i],Δx,kx[:,i],order=order_x)
                    uₓₓ[1:order_x,i] += SATₗ
                    uₓₓ[end-order_x+1:end,i] += SATᵣ
                end
            end
            if boundary_y != Periodic
                for i = 1:nx
                    uₓₓ[i,1:order_y]        += SAT(boundary_y,Left,gy(t),Δy,c=ky[i,:],order=order_y,forcing=true)
                    uₓₓ[i,end-order_y+1:end]+= SAT(boundary_y,Right,gy(t),Δy,c=ky[i,:],order=order_y,forcing=true)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            else
                for i = 1:nx
                    SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δy,ky[i,:],order=order_y)
                    uₓₓ[i,1:order_y] += SATₗ
                    uₓₓ[i,end-order_y+1:end] += SATᵣ
                end
            end
            return uₓₓ
        end

        Hx = build_H(nx,order_x)*Δx
        Hy = build_H(ny,order_y)*Δy
        i = 1
        converged = true
        chkpnt = t_f/4.0

        uⱼ = zeros(nx,ny)
        if nprocs() > 1
            uⱼ = SharedArray(uⱼ)
        end

        while t ≤ t_f
            # uⱼ = copy(uₒ)
            uⱼ .= uₒ
            ### FORCING TERMS
            if boundary_x != Periodic
                for i = 1:ny
                    # Fₗ = SAT(boundary_x,Left,gx(t),Δx,c=kx[:,i],order=order_x,forcing=true)
                    uⱼ[1:order_x,i] .= Δt*SAT_Dirichlet!(uⱼ[1:order_x,i],Left,gx(t),Δx,c=kx[1:order_x,i],order=order_x,forcing=true)
                    # Fᵣ = SAT(boundary_x,Right,gx(t),Δx,c=kx[:,i],order=order_x,forcing=true)
                    uⱼ[nx-order_x+1:nx,i] .= Δt*SAT_Dirichlet!(uⱼ[nx-order_x+1:nx,i],Right,gx(t),Δx,c=kx[nx-order_x+1:nx,i],order=order_x,forcing=true)
                    # uⱼ[1:order_x,i] += Δt*Fₗ
                    # uⱼ[nx-order_x+1:nx,i] += Δt*Fᵣ
                end
            end
            if boundary_y != Periodic
                for i = 1:nx
                    Fₗ = SAT(boundary_y,Left,gy(t),Δy,c=ky[i,:],order=order_y,forcing=true)
                    Fᵣ = SAT(boundary_y,Right,gy(t),Δy,c=ky[i,:],order=order_y,forcing=true)
                    uⱼ[i,1:order_y] += Δt*Fₗ
                    uⱼ[i,end-order_y+1:end] += Δt*Fᵣ
                end
            end
            # Advance the solution in time
            uₙ, converged = conj_grad(uⱼ,uⱼ,cgRHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,gx,gy,Hx,Hy,tol=tol,rtol=rtol,maxIT=maxIT)
            # If a global penalty function has been applied, update the solution
            if parallel_penalty
                uₙ,tmp = penalty_fn(uₙ,uₒ,Δt)
            end

            if converged
                # If CG converged store and update the solution, increase the time step if adaptive time stepping on
                # soln.u = cat(soln.u,uₙ,dims=3)
                if (samplefactor - t) ≤ 0.0
                    push!(soln.u,uₙ)
                    samplefactor += samplefactor_base
                    append!(soln.t,t)
                    append!(soln.Δt,Δt)
                end
                uₒ = uₙ
                i += 1
                t += Δt
                # append!(umw,tmp)
                if adaptive #if adaptive time stepping is turned on
                    if Δt < 300*Δt₀
                        Δt *= 1.05
                    end
                end
            else
                # If CG fails, restart the step. If the step has been tried before abort the time solve, returning the current solution.
                if Δt ≤ Δt₀
                    warnstr = string("Could not progress with minimum time step Δt₀=",Δt₀,", exiting.")
                    @warn(warnstr)
                    break
                end
                # Δt = Δt₀
                Δt = Δt/2.0
            end

            if checkpointing
                if chkpnt < t
                    push!(soln.u,uₙ)
                    append!(soln.t,t)
                    append!(soln.Δt,Δt)

                    chkpnt += t_f/4.0
                    save_object("chkpnt.jld2",soln)
                    println("checkpointing at t=",t)
                end
            end

        end
    end

    # if checkpoint
    #     checkpt = checkpoint_2d(soln,Δx,Δy,t_f,kx,ky,gx,gy,parallel_penalty,penalty_fn,order_x,order_y,method,maxIT,samplefactor,tol,rtol,adaptive)
    #     return checkpt
    # end

    if length(soln.u) == 1
        push!(soln.u,uₙ)
        append!(soln.t,t)
        append!(soln.Δt,Δt)
    end

    return soln, umw
end
