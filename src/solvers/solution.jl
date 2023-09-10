

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
        PT<:Union{PDEProblem,newPDEProblem}}
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
    function solution{TT}(grid::GridType,t::TT,Δt::TT,prob::newPDEProblem;preallocate=false) where TT
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
            elseif typeof(grid) <: GridMultiBlock{TT,1}
                u = [prob.InitialCondition(grid.Grids[I].grid) for I in eachgrid(grid)]
            end


            new{TT,typeof(u),typeof(grid),typeof(prob)}([u],grid,[Δt],[t],prob,0.0)
        end

    end
end


function setInitialCondition end
setInitialCondition(IC,G::Grid1D) = IC.(G.grid)
function setInitialCondition(IC,G::GridMultiBlock{TT,1}) where TT
    u = [zeros(TT,length(G.Grids[i])) for i in eachindex(G.Grids)]
    for j in eachindex(G.Grids)
        for i in eachindex(G.Grids[j].grid)
            u[i,j] = IC(G.Grids[j].grid[i])
        end
    end
end
function setInitialCondition(IC,G::Grid2D{TT}) where TT
    u = [zeros(TT,(G.Grids[i].grid.nx,G.Grids[i].grid.ny)) for i in eachindex(G.Grids)]
    for j = 1:G.inds[2,end]
        for i = 1:G.inds[1,end]
            u[i,j] = IC(G[i,j],G[i,j])
        end
    end
    return u
end



function UpdateSolution!(soln::solution{T,AT},u::AbstractArray{T},t::T,Δt::T) where {T,AT}
    push!(soln.u,u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)
end




