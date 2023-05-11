








"""
    solve

1. `solve(Prob::VariableCoefficientPDE1D{T},grid::GridType{T,1},Δt::T,t_f::T,solver::Symbol;...)`
2. `solve(Prob::VariableCoefficientPDE2D{T},grid::GridType{T,2},Δt::T,t_f::T,solver::Symbol;...)`

Solve the 1D or 2D PDE on the given grid required inputs:
- Prob: [`VariableCoefficientPDE1D`](@ref) or [`VariableCoefficientPDE2D`](@ref)
- grid: [`Grid1D`](@ref) or [`Grid2D`](@ref)
- `Δt` - the time step
- `t_f` - the final time
- `solver::Symbol` - must be `:cgie` currently.

Optional inputs:
- `adaptive::Bool=false` `true` or `false` - Adaptive time stepping
- `source::Union{Nothing,Function}=nothing` - Include a source term
- `penalty_func::Union{Nothing,Function}=nothing` - Include a 3D penalty function - 2D PROBLEMS ONLY


"""
function solve end
#= 1D SOLVER =#
function solve(Prob::VariableCoefficientPDE1D,grid::GridType{T,1},Δt::T,t_f::T,solver::Symbol;adaptive::Bool=false,source::Union{Nothing,Function}=nothing,penalty_func::Union{Nothing,Function}=nothing,Pgrid::Union{Nothing,ParallelGrid}=nothing,interpfn::Union{Nothing,Function}=nothing,sample_rate::Float64=0.0) where T

    target_state = 0.0
    if t_f == Inf
        target_state = 1e-5
        println("Going for steady state at rel-error Δu=",target_state)
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
    end

    DBlock = DataBlock{T}(Prob,grid,Δt)
    CGBlock = ConjGradBlock{T}(grid,Prob.order)
    soln = solution{T}(grid,0.0,Δt,Prob)
    BoundaryConditions = Prob.BoundaryConditions

    if typeof(Pgrid) <: ParallelGrid
        penalty_func = generate_parallel_penalty(Pgrid,grid,Prob.order)
    end
    typeof(penalty_func) <: Nothing ? penalty_function_enabled = false : penalty_function_enabled = true


    DBlock.u .= soln.u[1]
    
    
    # Build functions
    if Prob.BoundaryConditions.Left.type == Periodic
        _,SAT_LR    = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    elseif Prob.BoundaryConditions.Left.type != Periodic
        _,SAT_Left  = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
        _,SAT_Right = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    end
    Diff = generate_SecondDerivative(grid.n,grid.Δx,Prob.order)
    
    if solver == :cgie
        # Replace cache with the derivative of v + SATs
        function CGRHS!(cache::AbstractArray{T},u::AbstractArray{T},k::AbstractArray{T})
            Diff(cache,u,k)
            if Prob.BoundaryConditions[1].type != Periodic
                SAT_Left(cache,u,k,SolutionMode)
                SAT_Right(cache,u,k,SolutionMode)
            else
                SAT_LR(cache,u,k)
            end
        end
    elseif solver != :cgie
        integrate = ExplicitBlock{T}(grid,Δt,:RK4)
        function ExplicitRHS!(cache::AbstractArray{T},u::AbstractArray{T},k::AbstractArray{T},t::T) where T
            Diff(cache,u,k)
            if Prob.BoundaryConditions[1].type != Periodic
                SAT_Left(cache,u,k,t)
                SAT_Right(cache,u,k,t)
            # else
                # SAT_LR(cache,u,k)
            end
        end
    end
    
    
    t = Δt
    Δt₀ = Δt
    DBlock.uₙ₊₁ .= DBlock.u
    CGBlock.b .= DBlock.u
    
    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    
    
    while t < t_f

        
        
        if solver == :cgie
            if Prob.BoundaryConditions[1].type != Periodic
                setBoundary!(BoundaryConditions.Left.RHS,DBlock.boundary.RHS_Left,t,Δt)
                setBoundary!(BoundaryConditions.Right.RHS,DBlock.boundary.RHS_Right,t,Δt)
    
                SAT_Left(CGBlock.b, DBlock.boundary.RHS_Left, DBlock.K, DataMode)
                SAT_Right(CGBlock.b, DBlock.boundary.RHS_Right, DBlock.K, DataMode)
            end
            if typeof(source) <: Function
                addSource!(source,CGBlock.b,grid,t,Δt)
            end

            conj_grad!(CGRHS!,DBlock,CGBlock,Δt)

            if CGBlock.converged | !adaptive #If CG converges

                if penalty_function_enabled
                    penalty_func(DBlock.uₙ₊₁,DBlock.u,Δt)
                end

                # USED FOR DETERMINING EQUILIBRIUM
                DBlock.Δu = norm(DBlock.u .- DBlock.uₙ₊₁)/norm(DBlock.u)
                if (DBlock.Δu ≤ target_state) & (t_f == Inf)
                    t_f = t
                end

                DBlock.u .= DBlock.uₙ₊₁
                CGBlock.b .= DBlock.uₙ₊₁
                copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
                if adaptive & (Δt<300Δt₀)
                    Δt *= 1.05
                end
                t += Δt

                

            else #If CG fails, reset and retry step
                DBlock.uₙ₊₁ .= DBlock.u
                Δt /= 2.0
                CGBlock.converged = true
                if Δt < Δt₀/10.0
                    error("CG could not converge, aborting at t=",t," with Δt=",Δt)
                end
            end
        elseif solver == :RK4
            integrate(ExplicitRHS!,DBlock,t)
            t += Δt
        end
        # if sample < t
        #     sample += sample
        #     UpdateSolution!(soln,DBlock.u,t,Δt)
        # end 
    end

    push!(soln.u,DBlock.u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)

    return soln

end
#= 2D SOLVER =#
function solve(Prob::VariableCoefficientPDE2D,grid::GridType{T,2},Δt::T,t_f::T,solver::Symbol;adaptive::Bool=false,penalty_func::Union{Nothing,Function}=nothing,Pgrid::Union{Nothing,ParallelGrid}=nothing,source::Union{Nothing,Function}=nothing) where T

    target_state = 0.0
    if t_f == Inf
        target_state = 1e-5
        println("Going for steady state at rel-error Δu=",target_state)
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
    end

    DBlock = DataBlock{T}(Prob,grid,Δt,Prob.Kx,Prob.Ky)
    CGBlock = ConjGradBlock{T}(grid,Prob.order)
    soln = solution{T}(grid,0.0,Δt,Prob)

    if typeof(Pgrid) <: ParallelGrid
        penalty_func = generate_parallel_penalty(Pgrid,grid,Prob.order)
    end
    typeof(penalty_func) <: Nothing ? penalty_function_enabled = false : penalty_function_enabled = true

    DBlock.u .= soln.u[1]

    if Prob.BoundaryConditions.Left.type != Periodic
        _,SAT_Left  = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
        _,SAT_Right = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    else
        _,SAT_LR    = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    end
    if Prob.BoundaryConditions.Up.type != Periodic
        _,SAT_Up    = SAT(Prob.BoundaryConditions.Up,grid,Prob.order,solver)
        _,SAT_Down  = SAT(Prob.BoundaryConditions.Down,grid,Prob.order,solver)
    else
        _,SAT_UD    = SAT(Prob.BoundaryConditions.Up,grid,Prob.order,solver)
    end

    Diff = generate_SecondDerivative(grid.nx,grid.ny,grid.Δx,grid.Δy,Prob.order)

    function CGRHS!(cache::AbstractArray,u::AbstractArray,K::AbstractArray)
        Diff(cache,u,K[1],K[2])
        if Prob.BoundaryConditions.Left.type != Periodic #Left/Right boundaries
            SAT_Left(cache,u,K[1],SolutionMode)
            SAT_Right(cache,u,K[1],SolutionMode)
        else
            SAT_LR(cache,u,K[1]) #Periodic SAT
        end
        if Prob.BoundaryConditions.Up.type != Periodic #Up/Down boundaries
            SAT_Up(cache,u,K[2],SolutionMode)
            SAT_Down(cache,u,K[2],SolutionMode)
        else
            SAT_UD(cache,u,K[2]) #Periodic SAT
        end
    end

    t = Δt
    Δt₀ = Δt
    DBlock.uₙ₊₁ .= DBlock.u
    CGBlock.b .= DBlock.u

    # tmpu = zeros(T,(grid.nx,grid.ny))

    # kfin = t_f/Δt
    
    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    while t ≤ t_f
    # for k = 1:nf
        # t = k*Δt
        # t += Δt

        if Prob.BoundaryConditions.Left.type != Periodic #Left/Right boundaries
            setBoundary!(Prob.BoundaryConditions.Left.RHS,DBlock.boundary.RHS_Left,grid.gridy,grid.ny,t,Δt)
            setBoundary!(Prob.BoundaryConditions.Right.RHS,DBlock.boundary.RHS_Right,grid.gridy,grid.ny,t,Δt)

            SAT_Left(CGBlock.b, DBlock.boundary.RHS_Left, DBlock.K[1],DataMode)
            SAT_Right(CGBlock.b, DBlock.boundary.RHS_Right, DBlock.K[1],DataMode)
        end
        if Prob.BoundaryConditions.Up.type != Periodic #Up/Down boundaries
            setBoundary!(Prob.BoundaryConditions.Up.RHS,DBlock.boundary.RHS_Up,grid.gridx,grid.nx,t,Δt)
            setBoundary!(Prob.BoundaryConditions.Down.RHS,DBlock.boundary.RHS_Down,grid.gridx,grid.nx,t,Δt)

            SAT_Up(CGBlock.b, DBlock.boundary.RHS_Up, DBlock.K[2],DataMode)
            SAT_Down(CGBlock.b, DBlock.boundary.RHS_Down, DBlock.K[2],DataMode)
        end
        if typeof(source) <: Function
            addSource!(source,CGBlock.b,grid,t,Δt)
        end
        conj_grad!(CGRHS!,DBlock,CGBlock,Δt,warnings=false)
        
        if CGBlock.converged | !adaptive
            # If CG converges OR adaptive time stepping is off
            if penalty_function_enabled # Add parallel penalty
                # println(DBlock.uₙ₊₁[1,10])
                # tmpu .= DBlock.uₙ₊₁
                penalty_func(DBlock.uₙ₊₁,DBlock.u,Δt)
                # println(DBlock.uₙ₊₁[1,10])
            end
            # USED FOR DETERMINING EQUILIBRIUM
            DBlock.Δu = norm(DBlock.u .- DBlock.uₙ₊₁)/norm(DBlock.u)
            if (DBlock.Δu ≤ target_state) & (t_f == Inf)
                t_f = t
            end
            DBlock.u .= DBlock.uₙ₊₁
            CGBlock.b .= DBlock.uₙ₊₁
            # copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
            if adaptive & (Δt<300Δt₀)
                Δt *= 1.05
            end
            t += Δt
        else
            # If adaptive time stepping is turned on and CG fails
            DBlock.uₙ₊₁ .= DBlock.u
            CGBlock.b .= DBlock.u
            Δt = Δt/2.0
            CGBlock.converged = true
            if (Δt < Δt₀/10.0) #If Δt ~ 𝒪(Δt₀/10)
                error("CG could not converge, aborting at t=",t," with Δt=",DBlock.Δt)
            end
        end
    end

    push!(soln.u,DBlock.u)
    push!(soln.Δt,Δt)
    push!(soln.t,t-Δt)
    soln.Δu = DBlock.Δu
    return soln
end


