







"""
    SolverData{method,adaptive}
"""
struct SolverData{method,TT} 
    adaptive    :: Bool
    target      :: TT
    parallel    :: Bool
    function SolverData{TT}(;adaptive=false,target=0,method=:cgie,parallel=false) where TT
        new{method,TT}(adaptive,target,parallel)
    end
end






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
function solve(P::newPDEProblem{TT,DIM},G::GridType{TT,DIM},Δt::TT,t_f::TT;
        solver::Symbol=:cgie,adaptive::Bool=false,θ=TT(1)) where {TT,DIM}

    P.Parallel === nothing ? parallel = false : parallel = true 
    if solver ∈ [:cgie,:cn,:theta]
        if solver == :cgie
            θ == TT(1)
        elseif solver == :cn
            θ = TT(1/2)
        elseif solver == :theta
            θ = θ
        end
        SD = SolverData{TT}(adaptive=adaptive,method=solver,parallel=parallel)
        DBlock = DataMultiBlock(P,G,Δt,0.0,θ=θ)
        soln = solution(G,0.0,Δt,P)

        if typeof(G) <: LocalGridType
            DBlock[1].u .= soln.u[1]
            # DBlock[1].uₙ₊₁ .= soln.u[1]
        elseif typeof(G) <: GridMultiBlock
            for I in eachblock(DBlock)
                DBlock[I].u .= soln.u[1][I]
                # DBlock[I].uₙ₊₁ .= soln.u[1][I]
            end
        end

        
        if parallel
            for I in eachblock(DBlock)
                compute_parallel_operator(DBlock[I].Parallel.w_b,DBlock[I].u,DBlock[I].Parallel)
            end
        end


        # implicitsolve(P,G,Δt,t_f,SD)
        implicitsolve(soln,DBlock,G,Δt,t_f,SD)
    else
        error("solver not defined, try :cgie or :cn")
    end
end

# function implicitsolve(P::newPDEProblem{TT,DIM},G::GridType,Δt::TT,t_f::TT,solverconfig::SolverData{:cgie}) where {TT,DIM}
function implicitsolve(soln,DBlock,G,Δt::TT,t_f::TT,solverconfig::SolverData) where {TT}

    target_state = TT(1)

    t = TT(0)
    Δt₀ = Δt

# @show size(soln.u[1])

    copyto!(:uₙ₊₁,  :u, DBlock)
    
    # while t < t_f+Δt/2
    nt = round(t_f/Δt)
    # nt = 1;
    for i = 0:nt-1
        t = i*Δt

        theta_method(DBlock,t,Δt)

        if DBlock.SC.converged | !solverconfig.adaptive #If CG converges
            if solverconfig.parallel
                # println("b4",norm(DBlock[1].uₙ₊₁))
                applyParallelPenalty!(DBlock[1].uₙ₊₁,DBlock[1].u,DBlock.SC.Δt,DBlock.SC.θ,DBlock[1].Parallel)
                # println("afta",norm(DBlock[1].uₙ₊₁))
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
            if solverconfig.adaptive & (Δt<300Δt₀)
                Δt *= 1.05
            end
            DBlock.SC.t += DBlock.SC.Δt
            t = DBlock.SC.t

            # push!(soln.τ_hist,DBlock[1].Parallel.τ_i[1])

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
        push!(soln.u,DBlock[1].uₙ₊₁)
    else
        outu = [DBlock[I].u for I in eachblock(DBlock)]
        push!(soln.u,outu)
    end
    push!(soln.t,t)
    push!(soln.Δt,Δt)

    return soln

end
