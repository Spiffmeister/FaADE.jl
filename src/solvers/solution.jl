


mutable struct solution{T}
    u       :: Vector{AbstractArray{T}}
    grid    :: GridType
    Δt      :: Union{T,Vector{T}}
    t       :: Vector{T}
    problem :: PDEProblem

    function solution{T}(grid::GridType,t::T,Δt::T,prob::PDEProblem;preallocate=false) where T
        if preallocate
            N = ceil(Int64,t/Δt)
            n = length(x)
            u = [zeros(Float64,n) for _ in 1:N]

            u[1] = u₀

            new(u,grid,Δt,collect(range(0.0,t,length=N)))
        else #If an adaptive time step is being used, preallocation is impossible

            if typeof(grid) <: Grid1D
                u = prob.InitialCondition(grid.grid)
            elseif typeof(grid) <: Grid2D
                u = zeros(T,(grid.nx,grid.ny))
                for j = 1:grid.ny
                    for i = 1:grid.nx
                        u[i,j] = prob.InitialCondition.(grid.gridx[i],grid.gridy[j])
                    end
                end
            end


            new{T}([u],grid,[Δt],[t],prob)
        end

    end
end




