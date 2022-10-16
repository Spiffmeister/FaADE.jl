








"""
    solve

Solve the 1D or 2D PDE on the given grid
"""
function solve end
#= 1D SOLVER =#
function solve(Prob::VariableCoefficientPDE1D{T},grid::GridType{T,1},Δt::T,t_f::T,solver::Symbol;adaptive::Bool=false,source::Union{Nothing,Function}=nothing) where T

    DBlock = DataBlock{T}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.K)
    CGBlock = ConjGradBlock{T}(grid)
    soln = solution{T}(grid,0.0,Δt,Prob)

    DBlock.u .= soln.u[1]
    
    
    # Build functions
    if Prob.BoundaryConditions.Left.type == Periodic
        _,SAT_LR    = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    elseif Prob.BoundaryConditions.Left.type != Periodic
        _,SAT_Left  = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
        _,SAT_Right = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    end
    Diff = generate_Derivative(grid.n,grid.Δx,Prob.order)
    
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

    t = Δt
    Δt₀ = Δt
    DBlock.uₙ₊₁ .= DBlock.u

    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    
    while t < t_f

        if Prob.BoundaryConditions[1].type != Periodic
            DBlock.boundary.RHS_Left = Prob.BoundaryConditions.Left.RHS(t)
            DBlock.boundary.RHS_Right = Prob.BoundaryConditions.Right.RHS(t)

            SAT_Left(DBlock.boundary.SAT_Left, Δt*DBlock.boundary.RHS_Left,DBlock.K,DataMode)
            SAT_Right(DBlock.boundary.SAT_Right, Δt*DBlock.boundary.RHS_Right,DBlock.K,DataMode)
        end
        
        copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary,Prob.order)

        conj_grad!(DBlock,CGBlock,CGRHS!,Δt,Prob.order)

        if CGBlock.converged #If CG converges
            DBlock.u .= DBlock.uₙ₊₁
            copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
            if adaptive & (Δt<300Δt₀)
                Δt *= 1.05
            end
            t += Δt
        else #If CG fails, reset and retry step
            DBlock.uₙ₊₁ .= DBlock.u
            if Δt < Δt₀/10.0
                error("CG could not converge, aborting at t=",t," with Δt=",Δt)
            end
        end
    end

    push!(soln.u,DBlock.u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)

    return soln

end
#= 2D SOLVER =#
function solve(Prob::VariableCoefficientPDE2D{T},grid::GridType{T,2},Δt::T,t_f::T,solver::Symbol;adaptive::Bool=false,penalty_func::Union{Nothing,Function}=nothing,source::Union{Nothing,Function}=nothing) where T

    DBlock = DataBlock{T}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.Kx,Prob.Ky)
    CGBlock = ConjGradBlock{T}(grid)
    soln = solution{T}(grid,0.0,Δt,Prob)

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

    Diff = generate_Derivative(grid.nx,grid.ny,grid.Δx,grid.Δy,Prob.order)

    function CGRHS!(cache::AbstractArray,u::AbstractArray,K::AbstractArray)
        Diff(cache,u,K[1],K[2])
        if Prob.BoundaryConditions.Left.type != Periodic
            SAT_Left(cache,u,K[1],SolutionMode)
            SAT_Right(cache,u,K[1],SolutionMode)
        else
            SAT_LR(cache,u,K[1])
        end
        if Prob.BoundaryConditions.Up.type != Periodic
            SAT_Up(cache,u,K[2],SolutionMode)
            SAT_Down(cache,u,K[2],SolutionMode)
        else
            SAT_UD(cache,u,K[2])
        end
    end

    t = Δt
    Δt₀ = Δt
    DBlock.uₙ₊₁ .= DBlock.u
    
    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)

    while t < t_f

        if Prob.BoundaryConditions.Left.type != Periodic
            DBlock.boundary.RHS_Left .= Δt*Prob.BoundaryConditions.Left.RHS(t)
            DBlock.boundary.RHS_Right .= Δt*Prob.BoundaryConditions.Right.RHS(t)

            SAT_Left(DBlock.uₙ₊₁, DBlock.boundary.RHS_Left, DBlock.K[1],DataMode)
            SAT_Right(DBlock.uₙ₊₁, DBlock.boundary.RHS_Right, DBlock.K[1],DataMode)
        end
        if Prob.BoundaryConditions.Up.type != Periodic
            DBlock.boundary.RHS_Up .= Δt*Prob.BoundaryConditions.Up.RHS(t)
            DBlock.boundary.RHS_Down .= Δt*Prob.BoundaryConditions.Down.RHS(t)

            SAT_Up(DBlock.uₙ₊₁, DBlock.boundary.RHS_Up, DBlock.K[2],DataMode)
            SAT_Down(DBlock.uₙ₊₁, DBlock.boundary.RHS_Down, DBlock.K[2],DataMode)
        end

        # copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary,Prob.order)

        conj_grad!(DBlock,CGBlock,CGRHS!,Δt,Prob.order)

        if penalty_function_enabled
            penalty_func(DBlock.uₙ₊₁,DBlock.u,Δt)
        end

        if CGBlock.converged
            DBlock.u .= DBlock.uₙ₊₁
            # copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
            if adaptive & (Δt<300Δt₀)
                Δt *= 1.05
            end
            t += Δt
        else
            DBlock.uₙ₊₁ .= DBlock.u
            Δt = Δt/2.0
            CGBlock.converged = true
            if (Δt < Δt₀/10.0) | !adaptive
                error("CG could not converge, aborting at t=",t," with Δt=",DBlock.Δt)
            end
        end

    end

    push!(soln.u,DBlock.u)
    push!(soln.Δt,Δt)
    push!(soln.t,t)
    return soln
end


