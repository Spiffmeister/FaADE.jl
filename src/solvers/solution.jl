

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
        PT<:newPDEProblem}
    u       :: Vector{AT}
    grid    :: GT
    Δt      :: Union{TT,Vector{TT}}
    t       :: Vector{TT}
    problem :: PT
    Δu      :: TT

    τ_hist  :: Vector{TT}
end
"""
    solution{TT}(grid::GridType,t::TT,Δt::TT,prob::PDEProblem;preallocate=false) where TT
    DEPRECATED
"""
# function solution{TT}(grid::GridType,t::TT,Δt::TT,prob::PDEProblem;preallocate=false) where TT
#     if preallocate
#         N = ceil(Int64,t/Δt)
#         n = length(x)
#         u = [zeros(Float64,n) for _ in 1:N]

#         u[1] = u₀

#         new{TT,typeof(u),typeof(grid),typeof(PT)}(u,grid,Δt,collect(range(0.0,t,length=N)))
#     else #If an adaptive time step is being used, preallocation is impossible

#         if typeof(grid) <: Grid1D
#             u = prob.InitialCondition.(grid.grid)
#         elseif typeof(grid) <: Grid2D
#             u = zeros(TT,size(grid))
#             for I in eachindex(grid)
#                 u[I] = prob.InitialCondition(grid[I]...)
#             end
#         end


#         return solution{TT,typeof(u),typeof(grid),typeof(prob)}([u],grid,[Δt],[t],prob,0.0)
#     end

# end
"""
    solution{TT}(grid::LocalGridType,t::TT,Δt::TT,prob::newPDEProblem) where TT
"""
function solution(grid::LocalGridType{TT},t::TT,Δt::TT,prob::newPDEProblem) where TT
    u = _setInitialCondition(prob.InitialCondition,grid)

    return solution{TT,typeof(u),typeof(grid),typeof(prob)}([u],grid,[Δt],[t],prob,0.0,Vector{TT}())
end
"""
    solution(G::GridMultiBlock{TT,1},t::TT,Δt::TT,prob::newPDEProblem) where TT
1 dimensional multiblock problems
"""
function solution(G::GridMultiBlock{TT,1},t::TT,Δt::TT,prob::newPDEProblem) where TT
    u = [prob.InitialCondition.(G.Grids[I].grid) for I in eachgrid(G)]
    
    return solution{TT,typeof(u),typeof(G),typeof(prob)}([u],G,[Δt],[t],prob,0.0,Vector{TT}())
end
"""
    solution(G::GridMultiBlock{TT,2},t::TT,Δt::TT,prob::newPDEProblem) where TT
2 dimensional multiblock problems
"""
function solution(G::GridMultiBlock{TT,2},t::TT,Δt::TT,prob::newPDEProblem) where TT
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




