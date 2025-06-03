







"""
    SolverData{method,TT}
"""
struct SolverData{TT} 
    adaptive    :: Bool
    target      :: TT
    parallel    :: Bool

    @doc """
        SolverData{TT}(;adaptive=false,target=0,method=:cgie,parallel=false) where TT
    """
    function SolverData{TT}(;adaptive=false,target=0,parallel=false) where TT
        new{TT}(adaptive,target,parallel)
    end
end






"""
    solve
Solve the 1D or 2D PDE on the given grid.
"""
function solve end
"""
    solve(P::PDEProblem{TT,DIM},G::GridType{TT,DIM},Δt::TT,t_f::TT;
        solver::Symbol=:cn,adaptive::Bool=false,θ=TT(1),target=TT(0)) where {TT,DIM}

Inputs:
- Prob: [`Problem1D`](@ref) or [`Problem2D`](@ref)
- grid: [`Grid1D`](@ref) or [`Grid2D`](@ref)
- `Δt` - the time step
- `t_f` - the final time
- `solver::Symbol` - must be `:cgie` currently.

Optional inputs:
- `adaptive::Bool=false` `true` or `false` - Adaptive time stepping
"""
function solve(P::PDEProblem{TT,DIM},G::GridType{TT,DIM},Δt::TT,t_f::TT;adaptive::Bool=false,θ=TT(1/2),target=TT(0)) where {TT,DIM}

    P.Parallel === nothing ? parallel = false : parallel = true 
    SD = SolverData{TT}(adaptive=adaptive,parallel=parallel,target=target)
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

    
    # if parallel
    #     for I in eachblock(DBlock)
    #         compute_parallel_operator(DBlock[I].Parallel.w_b,DBlock[I].u,DBlock[I].Parallel)
    #     end
    # end


    # implicitsolve(P,G,Δt,t_f,SD)
    implicitsolve(soln,DBlock,G,Δt,t_f,SD)
end

"""
    restart(S::solution{TT,GT,PT};adaptive=false,Δt=TT(0),t_f=TT(0),K=(TT(0),TT(0)),K_para=TT(0),θ=TT(0.5)) where {TT,GT,PT}


"""
function restart(S::solution{TT,GT,PT};adaptive=false,Δt=TT(0),t_f=TT(0),K=(TT(0),TT(0)),K_para=TT(0),θ=TT(0.5)) where {TT,GT,PT}

    S.problem.Parallel == nothing ? parallel = false : parallel = true

    P = S.problem

    
    if K_para == TT(0)
        K_para = P.Parallel.κ
    end
    nParallel = ParallelData(P.Parallel.PGrid,S.grid,P.order,κ=K_para)
    
    if K[1] == TT(0)
        K = (P.Kx,P.Ky)
    end
    nP = Problem2D(P.order,P.InitialCondition,K[1],K[2],S.grid,P.BoundaryConditions,P.source.source,nParallel)

    if t_f == TT(0)
        t_f = S.t[end]
    end
    if Δt == TT(0)
        Δt = S.Δt[end]
    end
    

    SD = SolverData{TT}(adaptive=adaptive,method=:cn,parallel=parallel,target=TT(0))
    DBlock = DataMultiBlock(nP,S.grid,Δt,t_f,θ=θ)
    soln = solution(S.grid,0.0,Δt,nP)

    DBlock[1].u .= S.u[end]
    soln.u[1] .= S.u[end]


    


    implicitsolve(soln,DBlock,S.grid,Δt,t_f,SD)
end

# function implicitsolve(P::PDEProblem{TT,DIM},G::GridType,Δt::TT,t_f::TT,solverconfig::SolverData{:cgie}) where {TT,DIM}
function implicitsolve(soln,DBlock,G,Δt::TT,t_f::TT,solverconfig::SolverData) where {TT}

    if typeof(G) <: LocalGridType
        # uglobal = [zeros(size(G))] # TESTING
        # τglobal = zeros(length(uglobal))
        # Par = [DBlock[1].Parallel]

    else
        uglobal = [zeros(size(G.Grids[I])) for I in eachgrid(G)]
        # τglobal = zeros(length(uglobal))

        # Par = [DBlock[I].Parallel for I in eachblock(DBlock)]
    end
    # uglobal = [zeros(TT,size(G.Grids[I])) for I in eachgrid(G)]

    target_state = solverconfig.target

    t = TT(0)
    Δt₀ = Δt

    copyto!(:uₙ₊₁,  :u, DBlock)

    # while t < t_f+Δt/2
    nt = round(t_f/Δt)
    # nt = 1;
    for i = 0:nt-1
        t = i*Δt

        theta_method(DBlock,t,Δt)

        if DBlock.SC.converged | !solverconfig.adaptive #If perpendicular push converged
            if solverconfig.parallel
                if typeof(G) <: LocalGridType
                    if !(typeof(DBlock[1].Parallel.Interpolant) <: Nothing)
                        _update_ParallelMap(DBlock)
                    end
                    applyParallelPenalty!(DBlock[1].uₙ₊₁,DBlock.SC.t,DBlock.SC.Δt,DBlock[1].Parallel,DBlock[1].grid)
                    τ = DBlock[1].Parallel.τ
                else
                    # @show DBlock.ParallelData.Interpolant[1].G.Z[50,50]
                    # error()
                    _update_ParallelMap(DBlock)
                    # setglobalu!(uglobal,DBlock) # Circle case
                    # _setglobalu!(uglobal,DBlock)
                    # _updateCHSinterp(DBlock) # When CHS interpolation is used
                    # _setglobalu!(DBlock)
                    computeglobalw!(DBlock.ParallelData,uglobal,t,Δt)

                    τ = maximum(DBlock.ParallelData.τ)
                    for I in eachblock(DBlock)
                        applyParallelPenalty!(DBlock[I].uₙ₊₁,τ,DBlock.SC.Δt,DBlock.ParallelData.PData,DBlock[I].grid,I)
                    end
                end
            end

            # USED FOR DETERMINING EQUILIBRIUM
            # DBlock.SC.Δu = TT(0)
            for i in eachblock(DBlock)
                # DBlock.SC.Δu += norm(DBlock[i].u .- DBlock[i].uₙ₊₁)/norm(DBlock[i].u)
                relerr(DBlock[i])
            end
            collectΔu(DBlock)
            soln.Δu = norm(DBlock[1].u .- DBlock[1].uₙ₊₁)/norm(DBlock[1].u)
            DBlock.SC.Δu = soln.Δu

            # DBlock.Δu = norm(DBlock.u .- DBlock.uₙ₊₁)/norm(DBlock.u)
            if (DBlock.SC.Δu ≤ target_state) & (t_f == Inf)
                t_f = t
            end
            
            copyto!(:u,:uₙ₊₁,DBlock)
            if solverconfig.adaptive & (Δt<300Δt₀)
                DBlock.SC.Δt = min(DBlock.SC.Δt*1.05,300Δt₀)
            end
            DBlock.SC.t += DBlock.SC.Δt
            Δt = DBlock.SC.Δt
            t += Δt

            if !isnothing(DBlock[1].Parallel)
                # push!(soln.τ_hist,DBlock[1].Parallel.τ_i[1])
                push!(soln.τ_hist,τ)
            end
        else #If CG fails, reset and retry step
            # DBlock.uₙ₊₁ .= DBlock.u
            setValue(DBlock,DBlock,:uₙ₊₁,:u)
            DBlock.SC.Δt /= 2.0
            DBlock.SC.converged = true
            if DBlock.SC.Δt < Δt₀/10.0
                error("CG could not converge, aborting at t=",t," with Δt=",Δt)
            end
            t = DBlock.SC.t
            Δt = DBlock.SC.Δt
            t += Δt
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
