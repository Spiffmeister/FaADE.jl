

"""
    solution
Solution data structure, contains the initial condition and the solution at time of
 simulation termination. Also contains the grid structure and the PDE problem that was provided
 originally by the user.

Fields:
 - `u`, `grid`, `Δt`, `t`, `problem`, `Δu`
"""
mutable struct solution{TT,
        AT,
        GT<:GridType,
        PT<:PDEProblem}
    u       :: Vector{AT}
    grid    :: GT
    Δt      :: Union{TT,Vector{TT}}
    t       :: Vector{TT}
    problem :: PT
    Δu      :: TT

    τ_hist  :: Vector{TT}
end
"""
    solution{TT}(grid::LocalGridType,t::TT,Δt::TT,prob::PDEProblem) where TT
"""
function solution(grid::LocalGridType{TT},t::TT,Δt::TT,prob::PDEProblem) where TT
    u = _setInitialCondition(prob.InitialCondition,grid)

    return solution{TT,typeof(u),typeof(grid),typeof(prob)}([u],grid,[Δt],[t],prob,0.0,Vector{TT}())
end
"""
    solution(G::GridMultiBlock{TT,1},t::TT,Δt::TT,prob::PDEProblem) where TT
1 dimensional multiblock problems
"""
function solution(G::GridMultiBlock{TT,1},t::TT,Δt::TT,prob::PDEProblem) where TT
    u = [prob.InitialCondition.(G.Grids[I].grid) for I in eachgrid(G)]
    
    return solution{TT,typeof(u),typeof(G),typeof(prob)}([u],G,[Δt],[t],prob,0.0,Vector{TT}())
end
"""
    solution(G::GridMultiBlock{TT,2},t::TT,Δt::TT,prob::PDEProblem) where TT
2 dimensional multiblock problems
"""
function solution(G::GridMultiBlock{TT,2},t::TT,Δt::TT,prob::PDEProblem) where TT
    u = [zeros(TT,size(G.Grids[I])) for I in eachgrid(G)]

    for I in eachgrid(G)
        LG = G.Grids[I]
        for j = 1:LG.ny
            for i =1:LG.nx
                u[I][i,j] = prob.InitialCondition(LG[i,j]...)
            end
        end
    end

    return solution{TT,typeof(u),typeof(G),typeof(prob)}([u],G,[Δt],[t],prob,0.0,Vector{TT}())
end




function _setInitialCondition end
function _setInitialCondition(IC::Function,G::LocalGridType{TT}) where TT
    u = zeros(TT,size(G))
    for I in eachindex(G)
        u[I] = IC(G[I]...)
    end
    return u
end

function setInitialCondition end
function setInitialCondition(IC,G::GridMultiBlock{TT,1}) where TT
    u = [zeros(TT,length(G.Grids[i])) for i in eachindex(G.Grids)]
    for j in eachindex(G.Grids)
        for i in eachindex(G.Grids[j].grid)
            u[i,j] = IC(G.Grids[j].grid[i])
        end
    end
end


function UpdateSolution!(soln::solution{T,AT},u::AbstractArray{T},t::T,Δt::T) where {T,AT}
    push!(soln.u,u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)
end




