

"""
    solution
Solution data structure, contains the initial condition and the solution at time of
 simulation termination. Also contains the grid structure and the PDE problem that was provided
 originally by the user.

Fields:
 - `u`, `grid`, `Δt`, `t`, `problem`, `Δu`
"""
mutable struct solution{T,AT<:AbstractArray{T}}
    u       :: Vector{AT}
    grid    :: GridType
    Δt      :: Union{T,Vector{T}}
    t       :: Vector{T}
    problem :: PDEProblem
    Δu      :: T
    
    function solution{T}(grid::GridType,t::T,Δt::T,prob::PDEProblem;preallocate=false) where T
        if preallocate
            N = ceil(Int64,t/Δt)
            n = length(x)
            u = [zeros(Float64,n) for _ in 1:N]

            u[1] = u₀

            new{T,typeof(u)}(u,grid,Δt,collect(range(0.0,t,length=N)))
        else #If an adaptive time step is being used, preallocation is impossible

            if typeof(grid) <: Grid1D
                u = prob.InitialCondition.(grid.grid)
            elseif typeof(grid) <: Grid2D
                u = zeros(T,(grid.nx,grid.ny))
                for j = 1:grid.ny
                    for i = 1:grid.nx
                        u[i,j] = prob.InitialCondition.(grid.gridx[i],grid.gridy[j])
                    end
                end
            end


            new{T,typeof(u)}([u],grid,[Δt],[t],prob,0.0)
        end

    end
end






function UpdateSolution!(soln::solution{T,AT},u::AbstractArray{T},t::T,Δt::T) where {T,AT}
    push!(soln.u,u)
    push!(soln.t,t)
    push!(soln.Δt,Δt)
end




