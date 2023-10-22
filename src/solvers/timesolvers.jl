








"""
    solve
Solve the 1D or 2D PDE on the given grid.

Calls:
```julia
solve(Prob::VariableCoefficientPDE1D{T},grid::GridType{T,1},Δt::T,t_f::T,solver::Symbol;...)
solve(Prob::VariableCoefficientPDE2D{T},grid::GridType{T,2},Δt::T,t_f::T,solver::Symbol;...)
```

Inputs:
- Prob: [`VariableCoefficientPDE1D`](@ref) or [`VariableCoefficientPDE2D`](@ref)
- grid: [`Grid1D`](@ref) or [`Grid2D`](@ref)
- `Δt` - the time step
- `t_f` - the final time
- `solver::Symbol` - must be `:cgie` currently.

Optional inputs:
- `adaptive::Bool=false` `true` or `false` - Adaptive time stepping
- `source::Union{Nothing,Function}=nothing` - Include a source term
- `penalty_func::Union{Nothing,Function}=nothing` - penalty function for mapping points along field lines

TODO: Heavy optimisation required
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
    order = DerivativeOrder{Prob.order}()

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
    Diff = generate_SecondDerivative(grid.n,grid.Δx,order)
    
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
    # nt = 1000; for i = 1:nt
        
        
        if solver == :cgie
            if Prob.BoundaryConditions[1].type != Periodic
                setBoundaryConditions!(BoundaryConditions.Left.RHS,DBlock.boundary.RHS_Left,t,Δt)
                setBoundaryConditions!(BoundaryConditions.Right.RHS,DBlock.boundary.RHS_Right,t,Δt)
    
                SAT_Left(CGBlock.b, DBlock.boundary.RHS_Left, DBlock.K, DataMode)
                SAT_Right(CGBlock.b, DBlock.boundary.RHS_Right, DBlock.K, DataMode)
            end
            if typeof(source) <: Function
                addSource!(source,CGBlock.b,grid,t,Δt)
            end

            conj_grad!(CGRHS!,DBlock,CGBlock,Δt)

            if CGBlock.scalar.converged | !adaptive #If CG converges

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
function solve(Prob::VariableCoefficientPDE2D,grid::GridType{T,2},Δt::T,t_f::T,solver::Symbol;
        adaptive::Bool=false,penalty_func::Union{Nothing,Function}=nothing,Pgrid::Union{Nothing,ParallelGrid,ParallelData}=nothing,source::Union{Nothing,Function}=nothing,warnings=false,target=1e-8,nf=1) where T

    local btype1 :: BoundaryConditionType
    local btype2 :: BoundaryConditionType
    local penalty_function_enabled :: Bool

    target_state = 0.0
    if t_f == Inf
        target_state = target
        println("Going for steady state at rel-error Δu=",target_state)
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
        @warn "MAX ITERATIONS NOT SET"
    end

    DBlock = DataBlock{T}(Prob,grid,Δt)
    CGBlock = ConjGradBlock{T}(grid,Prob.order)
    soln = solution{T}(grid,T(0),Δt,Prob)

    if typeof(Pgrid) <: ParallelGrid
        # penalty_func = generate_parallel_penalty(Pgrid,grid,Prob.order)
        println("hi")
    end
    typeof(penalty_func) <: Nothing ? penalty_function_enabled = false : penalty_function_enabled = true

    DBlock.u .= soln.u[1]

    if Prob.BoundaryConditions.Left.type != Periodic
        _,SAT_Left!  = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
        # SAT_Left!(SATBL,grid,Prob.order,solver)
        _,SAT_Right! = SAT(Prob.BoundaryConditions.Right,grid,Prob.order,solver)
    else
        _,SAT_LR!    = SAT(Prob.BoundaryConditions.Left,grid,Prob.order,solver)
    end
    if Prob.BoundaryConditions.Up.type != Periodic
        _,SAT_Up!    = SAT(Prob.BoundaryConditions.Up,grid,Prob.order,solver)
        _,SAT_Down!  = SAT(Prob.BoundaryConditions.Down,grid,Prob.order,solver)
    else
        BSUD,SAT_UD!    = SAT(Prob.BoundaryConditions.Up,grid,Prob.order,solver)
    end

    Diff! = generate_SecondDerivative(grid.nx,grid.ny,grid.Δx,grid.Δy,Prob.order)

    btype1 = Prob.BoundaryConditions.Left.type :: BoundaryConditionType
    btype2 = Prob.BoundaryConditions.Up.type :: BoundaryConditionType


    CGRHS! = let btype1 = Prob.BoundaryConditions.Left.type,
            btype2 = Prob.BoundaryConditions.Up.type,
            mode = SolutionMode
        function CGRHS!(cache::AT,u::AT,K::KT) where {AT,KT}
            Diff!(cache,u,K[1],K[2])
            if (btype1 != Periodic) #Left/Right boundaries
                SAT_Left!(cache,u,K[1],mode)
                SAT_Right!(cache,u,K[1],mode)
            else
                SAT_LR!(cache,u,K[1]) #Periodic SAT
            end
            if btype2 != Periodic #Up/Down boundaries
                SAT_Up!(cache,u,K[2],mode)
                SAT_Down!(cache,u,K[2],mode)
            else
                SAT_UD!(cache,u,K[2]) #Periodic SAT
                # BSUD(cache,u,K[2]) #Periodic SAT
            end
            cache
        end
    end

    t = Δt
    Δt₀ = Δt
    DBlock.uₙ₊₁ .= DBlock.u
    CGBlock.b .= DBlock.u
    if t_f != Inf
        tprint₀ = tprint = t_f/10.0
    else
        tprint₀ = tprint = 50.0
    end

    # tmpu = zeros(T,(grid.nx,grid.ny))

    # kfin = t_f/Δt
    
    copyUtoSAT!(DBlock.boundary,DBlock.u,Prob.order)
    while t ≤ t_f
    # for k = 1:nf
        # t = k*Δt

        if btype1 != Periodic #Left/Right boundaries
            setBoundaryConditions!(Prob.BoundaryConditions.Left.RHS,DBlock.boundary.RHS_Left,grid.gridy,grid.ny,t,Δt)
            setBoundaryConditions!(Prob.BoundaryConditions.Right.RHS,DBlock.boundary.RHS_Right,grid.gridy,grid.ny,t,Δt)

            SAT_Left!(CGBlock.b, DBlock.boundary.RHS_Left, DBlock.K[1],DataMode)
            SAT_Right!(CGBlock.b, DBlock.boundary.RHS_Right, DBlock.K[1],DataMode)
        end
        if btype2 != Periodic #Up/Down boundaries
            setBoundaryConditions!(Prob.BoundaryConditions.Up.RHS,DBlock.boundary.RHS_Up,grid.gridx,grid.nx,t,Δt)
            setBoundaryConditions!(Prob.BoundaryConditions.Down.RHS,DBlock.boundary.RHS_Down,grid.gridx,grid.nx,t,Δt)

            SAT_Up!(CGBlock.b, DBlock.boundary.RHS_Up, DBlock.K[2],DataMode)
            SAT_Down!(CGBlock.b, DBlock.boundary.RHS_Down, DBlock.K[2],DataMode)
        end
        if typeof(source) <: Function
            addSource!(source,CGBlock.b,grid,t,Δt)
        end
        conj_grad!(CGRHS!,DBlock,CGBlock,Δt,warnings=warnings)
        
        if CGBlock.scalar.converged | !adaptive
            # If CG converges OR adaptive time stepping is off
            # println(Pgrid)
            if typeof(Pgrid) <: ParallelData
                # println("hi")
                applyParallelPenalty!(DBlock.uₙ₊₁,DBlock.u,Δt,Pgrid)
            elseif penalty_function_enabled
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
            # if tprint < t
            #     println("t=",t," out of t_f=",t_f,"     Δu=",DBlock.Δu)
            #     tprint += tprint₀
            # end
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


"""
    SolverData{method,adaptive}
"""
struct SolverData{method,TT} 
    adaptive    :: Bool
    target      :: TT
    function SolverData{TT}(;adaptive=false,target=0,method=:cgie) where TT
        new{method,TT}(adaptive,target)
    end
end

function solve(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM},Δt::TT,t_f::TT;
        solver::Symbol=:cgie,adaptive::Bool=false) where {TT,DIM}

        
    if solver == :cgie
        SD = SolverData{TT}(adaptive=adaptive,method=solver)
        DBlock = DataMultiBlock(P,G,Δt,0.0)
        soln = solution(G,0.0,Δt,P)

        if typeof(G) <: LocalGridType
            DBlock[1].u .= soln.u[1]
        elseif typeof(G) <: GridMultiBlock
            for I in eachblock(DBlock)
                DBlock[I].u .= soln.u[1][I]
            end
        end

        # implicitsolve(P,G,Δt,t_f,SD)
        implicitsolve(soln,DBlock,G,Δt,t_f,SD)
    end
end

# function implicitsolve(P::newPDEProblem{TT,DIM},G::GridType,Δt::TT,t_f::TT,solverconfig::SolverData{:cgie}) where {TT,DIM}
function implicitsolve(soln,DBlock,G,Δt::TT,t_f::TT,solverconfig::SolverData) where {TT}

    # target_state = 0.0
    # if t_f == Inf
    #     target_state = 1e-5
    #     println("Going for steady state at rel-error Δu=",target_state)
    #     @warn "MAX ITERATIONS NOT SET"
    #     @warn "MAX ITERATIONS NOT SET"
    #     @warn "MAX ITERATIONS NOT SET"
    # end
    target_state = TT(1)

    penalty_function_enabled = DBlock.parallel
    # if typeof(DBlock.parallel) <: Nothing
    #     penalty_function_enabled = false
    # end

    # DBlock = DataMultiBlock(P,G,Δt,0.0)
    # soln = solution(G,0.0,Δt,P)

    # if typeof(G) <: LocalGridType
    #     DBlock[1].u .= soln.u[1]
    # elseif typeof(G) <: GridMultiBlock
    #     for I in eachblock(DBlock)
    #         DBlock[I].u .= soln.u[1][I]
    #     end
    # end
    
    t = TT(0)
    Δt₀ = Δt

    copyto!(:uₙ₊₁,  :u, DBlock)
    copyto!(:b,     :u, DBlock)

    while t < t_f
        # for I in eachblock(DBlock)
        #     DBlock[I].SC.t = t
        #     DBlock[I].SC.Δt = Δt
        # end
        
        # setBoundaryConditions!(DBlock)
        # setDiffusionCoefficient!(DBlock)
        # applySATs(:b,DBlock,DataMode)
        # addSource!(:b,DBlock)

        if typeof(solverconfig) <: SolverData{:cgie}
            implicit_euler(DBlock,t,Δt)
        elseif typeof(solverconfig) <: SolverData{:cn}
            crank_nicholson(DBlock,t,Δt)
        end
        # conj_grad!(DBlock,Δt)
        # crank_nicholson(DBlock)
        # println(DBlock.SC.converged)
        if DBlock.SC.converged | !solverconfig.adaptive #If CG converges
            if penalty_function_enabled
                penalty_func(DBlock.uₙ₊₁,DBlock.u,Δt)
            end

            # USED FOR DETERMINING EQUILIBRIUM
            # DBlock.SC.Δu = TT(0)
            for i in eachblock(DBlock)
                # DBlock.SC.Δu += norm(DBlock[i].u .- DBlock[i].uₙ₊₁)/norm(DBlock[i].u)
                relerr(DBlock[i])
            end
            collectΔu(DBlock)

            # DBlock.Δu = norm(DBlock.u .- DBlock.uₙ₊₁)/norm(DBlock.u)
            if (DBlock.SC.Δu ≤ target_state) & (t_f == Inf)
                t_f = t
            end
            copyto!(:u,:uₙ₊₁,DBlock)
            copyto!(:b,:uₙ₊₁,DBlock)
            if solverconfig.adaptive & (Δt<300Δt₀)
                Δt *= 1.05
            end
            t += Δt

            

        else #If CG fails, reset and retry step
            # DBlock.uₙ₊₁ .= DBlock.u
            setValue(DBlock,DBlock,:uₙ₊₁,:u)
            DBlock.SC.Δt /= 2.0
            DBlock.SC.converged = true
            if DBlock.SC.Δt < Δt₀/10.0
                error("CG could not converge, aborting at t=",t," with Δt=",Δt)
            end
            # t += Δt
        end
        # if sample < t
        #     sample += sample
        #     UpdateSolution!(soln,DBlock.u,t,Δt)
        # end 
    end

    if typeof(G) <: LocalGridType
        push!(soln.u,DBlock[1].u)
    else
        outu = [DBlock[I].u for I in eachblock(DBlock)]
        push!(soln.u,outu)
    end
    push!(soln.t,t)
    push!(soln.Δt,Δt)

    return soln

end
