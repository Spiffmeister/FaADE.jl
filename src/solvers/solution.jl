

"""
    solution
Solution data structure, contains the initial condition and the solution at time of
 simulation termination. Also contains the grid structure and the PDE problem that was provided
 originally by the user.

Fields:
 - `u`, `grid`, `Δt`, `t`, `problem`, `Δu`
"""
mutable struct solution{TT,
        AT<:AbstractArray{TT},
        GT<:GridType,
        PT<:PDEProblem}
    u       :: Vector{AT}
    grid    :: GT
    Δt      :: Union{TT,Vector{TT}}
    t       :: Vector{TT}
    problem :: PT
    Δu      :: TT
    
    function solution{TT}(grid::GridType,t::TT,Δt::TT,prob::PDEProblem;preallocate=false) where TT
        if preallocate
            N = ceil(Int64,t/Δt)
            n = length(x)
            u = [zeros(Float64,n) for _ in 1:N]

            u[1] = u₀

            new{TT,typeof(u),typeof(grid),typeof(PT)}(u,grid,Δt,collect(range(0.0,t,length=N)))
        else #If an adaptive time step is being used, preallocation is impossible

            if typeof(grid) <: Grid1D
                u = prob.InitialCondition.(grid.grid)
            elseif typeof(grid) <: Grid2D
                u = zeros(TT,(grid.nx,grid.ny))
                for j = 1:grid.ny
                    for i = 1:grid.nx
                        u[i,j] = prob.InitialCondition.(grid.gridx[i],grid.gridy[j])
                    end
                end
            end


            new{TT,typeof(u),typeof(grid),typeof(prob)}([u],grid,[Δt],[t],prob,0.0)
        end

    end
end


setInitialCondition(IC,G::Grid1D)
setInitialCondition(IC,G::Grid1D) = IC.(G.grid)
# setInitialCondition(IC,G::GridMultiBlock{TT,1}) where TT = [IC(G[i]) for i in 1:length(G)]
function setInitialCondition(IC,G::Grid2D{TT}) where TT
    u = zeros(TT,(grid.nx,grid.ny))
    for j = 1:grid.ny
        for i = 1:grid.nx
            u[i,j] = IC(grid.gridx[i],grid.gridy[j])
        end
    end
    return u
end
function setInitialCondition(IC,G::GridMultiBlock{TT,1}) where {TT,1}
    u = zeros(TT,length(G))
    for i in eachindex(G)
        u[i] = IC(G[i])
    end
end




function UpdateSolution!(soln::solution{T,AT},u::AbstractArray{T},t::T,Δt::T) where {T,AT}
    push!(soln.u,u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)
end




