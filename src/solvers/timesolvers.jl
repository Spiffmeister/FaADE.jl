








"""
solve

1. Construction
a. DataBlock
b. SATs
c. Derivative function
2. Solving
"""
function solve end
#= 1D SOLVER =#
function solve(Prob::VariableCoefficientPDE1D,grid::GridType,Δt,t_f,solver;adaptive=false)

    DBlock = DataBlock{Float64}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.K)
    CGBlock = ConjGradBlock{Float64}(grid.n)
    
    soln = solution{Float64}(grid,0.0,Δt,Prob)

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
    function CGRHS!(cache::AbstractArray,u::AbstractArray,k::AbstractArray)
        Diff(cache,u,k)
        if Prob.BoundaryConditions[1].type != Periodic
            SAT_Left(cache,u,k,SolutionMode)
            SAT_Right(cache,u,k,SolutionMode)
        else
            SAT_LR(cache,u,k)
        end
    end

    t = Δt
    DBlock.uₙ₊₁ .= DBlock.u

    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    
    while t < t_f

        if Prob.BoundaryConditions[1].type != Periodic
            SAT_Left(DBlock.boundary.SAT_Left,Δt*Prob.BoundaryConditions.Left.RHS(t),DBlock.K,DataMode)
            SAT_Right(DBlock.boundary.SAT_Right,Δt*Prob.BoundaryConditions.Right.RHS(t),DBlock.K,DataMode)
        end
        
        copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary,Prob.order)

        conj_grad!(DBlock,CGBlock,CGRHS!,Δt,Prob.order)

        if CGBlock.converged #If CG converges
            DBlock.u .= DBlock.uₙ₊₁
            copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
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
#= 2D SOLVER =#
function solve(Prob::VariableCoefficientPDE2D,grid::GridType,Δt,t_f,solver;adaptive=false,penalty_fn!::Union{Nothing,Function}=nothing)

    DBlock = DataBlock{Float64}(Prob.BoundaryConditions,grid,Δt,Prob.order,Prob.Kx,Prob.Ky)
    CGBlock = ConjGradBlock{Float64}(grid.nx,grid.ny)

    soln = solution{Float64}(grid,0.0,Δt,Prob)

    typeof(penalty_fn!) <: Nothing ? penalty_function_enabled = false : penalty_function_enabled = true

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
    DBlock.uₙ₊₁ .= DBlock.u
    
    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)

    while t < t_f

        if Prob.BoundaryConditions.Left.type != Periodic
            SAT_Left(DBlock.uₙ₊₁, Δt*Prob.BoundaryConditions.Left.RHS(t), DBlock.K[1],DataMode)
            SAT_Right(DBlock.uₙ₊₁, Δt*Prob.BoundaryConditions.Right.RHS(t), DBlock.K[1],DataMode)
        end
        if Prob.BoundaryConditions.Up.type != Periodic
            SAT_Up(DBlock.uₙ₊₁, Δt*Prob.BoundaryConditions.Up.RHS(t), DBlock.K[2],DataMode)
            SAT_Down(DBlock.uₙ₊₁, Δt*Prob.BoundaryConditions.Down.RHS(t), DBlock.K[2],DataMode)
        end

        # copySATtoU!(DBlock.uₙ₊₁,DBlock.boundary,Prob.order)

        conj_grad!(DBlock,CGBlock,CGRHS!,Δt,Prob.order)

        if penalty_function_enabled
            penalty_fn!(DBlock.uₙ₊₁,DBlock.u,Δt)
        end

        if CGBlock.converged
            DBlock.u .= DBlock.uₙ₊₁
            # copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
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

    push!(soln.u,DBlock.u)
    return soln
end


