
struct ParGrid
    x   :: AbstractArray
    y   :: AbstractArray
end

"""
    ParallelGrid
Stores the current, forward and backward planes for the parallel tracing.

In 2D arrays are of format ``[(x_1,y_1),(x_1,y_2),...,(x_n,y_n)]``
"""
struct ParallelGrid{D}
    plane   :: AbstractArray
    Bplane  :: AbstractArray
    Fplane  :: AbstractArray
end



"""
    construct_grid(χ::Function,grid::Grid2D,z::Vector)
Constructs the backwards and forward planes for a given plane

Inputs:
- Field line ODE that returns ``[x,y]``
- GridType
- z values of planes to trace to

Outputs:
- ParallelGrid object (see [ParallelGrid](@ref))
"""
function construct_grid(χ::Function,grid::Grid2D,z::Vector;xmode=:stop,ymode=:period)

    modelist = [:stop,:period,:ignore]
    xmode ∈ modelist ? nothing : error("mode unavailable")
    ymode ∈ modelist ? nothing : error("mode unavailable")

    if typeof(grid) <: Grid2D
        xy = [[x,y] for x in grid.gridx for y in grid.gridy]
    end

    BPlane = construct_plane(χ,xy,z[1],(grid.nx,grid.ny))
    FPlane = construct_plane(χ,xy,z[2],(grid.nx,grid.ny))

    postprocess_plane!(BPlane,[grid.gridx[1],grid.gridx[end]],[grid.gridy[1],grid.gridy[end]],xmode,ymode)
    postprocess_plane!(FPlane,[grid.gridx[1],grid.gridx[end]],[grid.gridy[1],grid.gridy[end]],xmode,ymode)
    
    Pgrid = ParallelGrid{2}(reshape(xy,grid.nx,grid.ny),BPlane,FPlane)

    return Pgrid
end


"""
    construct_plane(χ,X,z,n)
Constructs the forward and backward planes for a single solution plane
"""
function construct_plane(χ::Function,X::AbstractArray{Vector{T}},z,n;periods=1) where T

    # plane = fill(zeros(T,2),n)

    function prob_fn(prob,i,repeat)
        remake(prob,u0=X[i])
    end
    P = ODEProblem(χ,X[1],(T(0),T(periods)*z))
    EP = EnsembleProblem(P,prob_func=prob_fn)

    sim = solve(EP,Tsit5(),EnsembleSerial(),trajectories=prod(n),save_on=false,save_end=true)
    
    plane = fill(zeros(T,2),prod(n))
    for i = 1:prod(n)
        plane[i] = sim.u[i].u[2]
    end

#=
    planex = zeros(T,n)
    planey = zeros(T,n)
    for i = 1:prod(n)
        planex[i] = sim.u[i].u[2][1]
        planex[i] = sim.u[i].u[2][2]
    end
    planex = reshape(planex,n[2],n[1])'
    planey = reshape(planey,n[2],n[1])'
=#

    plane = reshape(plane,n[2],n[1])'

    return plane
end


function postprocess_plane!(X,xbound,ybound,xmode,ymode)
    @views out_of_bounds!(X,xbound,xmode,1)
    @views out_of_bounds!(X,ybound,ymode,2)
end
"""
    out_of_bounds!(X,boundx,boundy)
Move out of bounds points to the boundary
"""
@views function out_of_bounds!(X,bound,mode,d)
    
    if mode == :stop
        for i in eachindex(X)
            if X[i][d] ≤ bound[1]
                X[i][d] = bound[1]
            elseif X[i][d] ≥ bound[2]
                X[i][d] = bound[2]
            end
        end
    elseif mode == :period
        for i in eachindex(X)
            X[i][d] = rem2pi(X[i][d],RoundNearest)
        end
    end
end


Base.show(io::IO, PG::ParallelGrid) = print(io, " generated parallel grid.")

