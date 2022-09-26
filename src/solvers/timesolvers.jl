







function solve end

"""
    solve

    1. Construction
        a. DataBlock
        b. SATs
        c. Derivative function
    2. Solving
"""
function solve(Prob::VariableCoefficientPDE1D,grid::GridType,Δt,t_f,solver;adaptive=false)

    DBlock = DataBlock{Float64}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.K)
    CGBlock = ConjGradBlock{Float64}(grid.n)
    
    soln = solution{Float64}(grid,0.0,Δt,Prob)

    DBlock.u .= soln.u[1]
    
    
    # Build functions
    _,SAT_Left    = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    _,SAT_Right   = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    Diff = generate_Derivative(grid.n,grid.Δx,Prob.order)
    
    # Replace cache with the derivative of v + SATs
    function CGRHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray)
        Diff(cache,u,k)
        SAT_Left(cache,u,k,SolutionMode)
        SAT_Right(cache,u,k,SolutionMode)
    end

    t = Δt
    DBlock.uₙ₊₁ .= DBlock.u

    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    # copyUtoSAT!(DBlock.boundary.SAT_Left,DBlock.u,Left,Prob.order)
    # copyUtoSAT!(DBlock.boundary.SAT_Right,DBlock.u,Right,Prob.order)
    
    while t < t_f
        SAT_Left(DBlock.boundary.SAT_Left,Δt*Prob.BoundaryConditions.Left.RHS(t),DBlock.K,DataMode)
        SAT_Right(DBlock.boundary.SAT_Right,Δt*Prob.BoundaryConditions.Right.RHS(t),DBlock.K,DataMode)
        
        copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary,Prob.order)

        # copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary.SAT_Left,Left,Prob.order)
        # copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary.SAT_Right,Right,Prob.order)

        conj_grad!(DBlock,CGBlock,CGRHS!,Δt,Prob.order)

        if CGBlock.converged #If CG converges
            DBlock.u .= DBlock.uₙ₊₁
            copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
            # copyUtoSAT!(DBlock.boundary.SAT_Left,DBlock.u,Left,Prob.order)
            # copyUtoSAT!(DBlock.boundary.SAT_Right,DBlock.u,Right,Prob.order)
            if adaptive
                Δt *= 1.05
            end
            t += Δt
        else #If CG fails, reset and retry step
            DBlock.uₙ₊₁ .= DBlock.u
            if DBlock.Δt < soln.Δt[1]/10.0
                error("CG could not converge, aborting at t=",t," with Δt=",DBlock.Δt)
            end
        end
    end

    push!(soln.u,DBlock.u)

    return soln

end
# 2D Mode
function solve(Prob::VariableCoefficientPDE2D)

    DBlock = DataBlock{Float64}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.K)
    CGBlock = ConjGradBlock{Float64}(grid.nx,grid.ny)

    soln = solution{Float64}(grid,0.0,Δt,Prob)

    DBlock.u .= soln.u[1]

    _,SAT_Left  = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    _,SAT_Right = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    _,SAT_Up    = SAT(Prob.BoundaryConditions.Up,grid,Prob.order,solver)
    _,SAT_Down  = SAT(Prob.BoundaryConditions.Down,grid,Prob.order,solver)

    Diff = generate_Derivative(grid.nx,grid.ny,grid.Δx,grid.Δy,Prob.order)

    function CGRHS!(cache::AbstractArray,u::AbstractArray,kx::AbstractArray,ky::AbstractArray)
        Diff(cache,u,kx,ky)
        SAT_Left(cache,u,K[1],K[2],SolutionMode)
        SAT_Right(cache,u,K[1],K[2],SolutionMode)
        SAT_Up(cache,u,K[1],K[2],SolutionMode)
        SAT_Down(cache,u,K[1],K[2],SolutionMode)
    end

    t = Δt
    DBlock.uₙ₊₁ .= DBlock.u
    
    copyUtoSAT!(DBlock,Prob.order)

    while t < t_f

        SAT_Left(DBlock.boundary.SAT_Left, Δt*Prob.BoundaryConditions.Left.RHS(t), DBlock.K,DataMode)
        SAT_Right(DBlock.boundary.SAT_Right, Δt*Prob.BoundaryConditions.Right.RHS(t), DBlock.K,DataMode)
        SAT_Up(DBlock.boundary.SAT_Up, Δt*Prob.BoundaryConditions.Up.RHS(t), DBlock.K,DataMode)
        SAT_Down(DBlock.boundary.SAT_Down, Δt*Prob.BoundaryConditions.Down.RHS(t), DBlock.K,DataMode)

        copySATtoU!(DBlock,Prob.order)

        conj_grad!()

        if CGBlock
            DBlock.u .= DBlock.uₙ₊₁
            copyUtoSAT!(DBlock,Prob.order)
            if adaptive
                Δt *= 1.05
            end
            t += Δt
        else
            DBlock.uₙ₊₁ .= DBlock.u
            if DBlock.Δt < soln.Δt[1]/10.0
                error("CG could not converge, aborting at t=",t," with Δt=",DBlock.Δt)
            end
        end

    end
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
function time_solver(u₀::Function,t_f::Float64,Δt::Float64,grid::Grid1D,BoundaryTerms;
        method::Symbol=:euler,α::Float64=1.5,tol::Float64=1e-5,maxIT::Int64=-1,warnings::Bool=false,samplefactor=1.0)

    uₙ = zeros(Float64,n)

    uₒ = u₀.(grid.domain)

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
            uₓₓ = Dₓₓ(uₓₓ,u,n,x,Δx,t,Δt,k,order=grid.order)

            # for SAT in BoundaryTerms.SATs
            #     if type != Periodic
            #         SAT(uₓₓ)
            #     end
            # end

            if boundary_left != Periodic
                uₓₓ[1:order]    += SAT(boundary_left,Left,u,Δx,g(t),c=k,order=order)
                uₓₓ[n-order+1:n]+= SAT(boundary_right,Right,u,Δx,g(t),c=k,order=order)
            else
                SATₗ,SATᵣ = SAT_Periodic(u,Δx,k,order=order)
                uₓₓ[1:order]        += SATₗ
                uₓₓ[end-order+1:end]+= SATᵣ
                # uₓₓ = SAT_Periodic!(uₓₓ,u,Δx,k,order=order)
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
            if boundary_left != Periodic
                SAT!(uₓₓ,boundary_left,Left,u[1:order],Δx,order=order,forcing=false)
                SAT!(uₓₓ,boundary_right,Right,u[n-order+1:n],Δx,order=order,forcing=false)
            else
                uₓₓ = SAT_Periodic!(uₓₓ,u,Δx,c=k,order=order)
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
            if boundary_left != Periodic
                SAT!(uⱼ,boundary_left,Left,Δt*boundary(t),Δx,order=order,c=k[1:order],forcing=true)
                SAT!(uⱼ,boundary_right,Right,Δt*boundary(t),Δx,order=order,c=k[n-order+1:n],forcing=true)
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
function time_solver(PDE::Function,u₀::Function,nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,x::AbstractVector{Float64},y::AbstractVector{Float64},t_f::Float64,Δt::Float64,kx::AbstractMatrix{Float64},ky::AbstractMatrix{Float64},BoundaryTerms,gy,boundary_y::BoundaryConditionType;
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
                    uₓₓ[:,i] = SAT_Periodic!(uₓₓ[:,i],u[:,i],Δx,kx[:,i],order=order_x)
                end
            end
            if boundary_y != Periodic
                for i = 1:nx #y boundaries
                    uₓₓ[i,1:order_y]            += SAT(boundary_y,Left,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                    uₓₓ[i,end-order_y+1:end]    += SAT(boundary_y,Right,u[i,:],Δy,gy(t),c=ky[i,:],order=order_y)
                end
            else
                for i = 1:nx
                    uₓₓ[i,:] = SAT_Periodic!(uₓₓ[i,:],u[i,:],Δy,c=ky[i,:],order=order_y)
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


        SATFns = construct_SAT(BoundaryTerms)

        #=
        DPenalties = BoundaryTerms.SATs[1].penalties
        BOp = BoundaryTerms.SATs[1].BDₓᵀ
        gx = BoundaryTerms.SATs[1].RHS

        FDL(A,B,C) = SAT_Dirichlet_implicit!(A,Left,B,C,Δx,DPenalties.α,DPenalties.τ,BOp,order_x,eachcol)
        FFDL(A,B,C) = SAT_Dirichlet_implicit_data!(A,Left,B,C,Δx,DPenalties.α,DPenalties.τ,BOp,order_x,eachcol)

        FDR(A,B,C) = SAT_Dirichlet_implicit!(A,Right,B,C,Δx,DPenalties.α,DPenalties.τ,BOp,order_x,eachcol)
        FFDR(A,B,C) = SAT_Dirichlet_implicit_data!(A,Right,B,C,Δx,DPenalties.α,DPenalties.τ,BOp,order_x,eachcol)
        =#

        function cgRHS(uₓₓ,u,kx,ky)
            PDE(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,order_x=order_x,order_y=order_y)
            ### SATs
            for ApproxTerm in SATFns
                ApproxTerm(uₓₓ,u,kx,SolutionMode)
            end

            
            # if boundary_x != Periodic
                # FDL(uₓₓ,u,kx)
                # FDR(uₓₓ,u,kx)
            # else
            if boundary_y != Periodic
                # for i = 1:nx
                    # uₓₓ[i,1:order_y]        += SAT(boundary_y,Left,u[i,1:order_y],Δy,c=ky[i,1:order_y],order=order_y)
                    # uₓₓ[i,ny-order_y+1:ny]  += SAT(boundary_y,Right,u[i,ny-order_y+1:ny],Δy,c=ky[i,ny-order_y+1:ny],order=order_y)
                    # SAT!(@views(uₓₓ[i,1:order_y]),boundary_y,Left,u[i,1:order_y],Δy,c=ky[i,1:order_y],order=order_y)
                    # SAT!(@views(uₓₓ[i,ny-order_y+1:ny]),boundary_y,Right,u[i,ny-order_y+1:ny],Δy,c=ky[i,ny-order_y+1:ny],order=order_y)
                # end
                SATAdd!(uₓₓ,boundary_y,Left,u,Δx,c=ky,order=order_y)
                SATAdd!(uₓₓ,boundary_y,Right,u,Δx,c=ky,order=order_y)
            else
                for i = 1:nx
                    # SATₗ,SATᵣ = SAT_Periodic(u[i,:],Δy,c=ky[i,:],order=order_y)
                    # uₓₓ[i,1:order_y] +=SATₗ
                    # uₓₓ[i,ny-order_y+1:ny] += SATᵣ
                    SAT_Periodic!(@views(uₓₓ[i,:]),u[i,:],Δy,c=ky[i,:],order=order_y)
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
            uⱼ .= uₒ
            ### FORCING TERMS
            for ApproxTerm in SATFns
                ApproxTerm(uⱼ,Δt*BoundaryTerms.SATs[1].RHS(t),kx,DataMode)
            end
            if boundary_y != Periodic
                for i = 1:nx
                    SAT!(@views(uⱼ[1:order_y,i]),boundary_y,Left,Δt*gy(t),Δy,c=ky[i,1:order_y],order=order_y,forcing=true)
                    SAT!(@view(uⱼ[ny-order_y+1:ny,i]),boundary_y,Right,Δt*gy(t),Δy,c=ky[ny-order_y+1:ny,i],order=order_y,forcing=true)
                end
            end
            # Advance the solution in time
            uₙ, converged = conj_grad(uⱼ,uⱼ,cgRHS,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky,Hx,Hy,tol=tol,rtol=rtol,maxIT=maxIT)
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


    if length(soln.u) == 1
        push!(soln.u,uₙ)
        append!(soln.t,t)
        append!(soln.Δt,Δt)
    end

    return soln, umw
end
