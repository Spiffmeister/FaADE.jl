

"""
    ParallelGrid
Stores the current, forward and backward planes for the parallel tracing.

In 2D arrays are of format ``[(x_1,y_1),(x_1,y_2),...,(x_n,y_n)]``
"""
struct ParallelGrid
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
    
    Pgrid = ParallelGrid(hcat(xy...),BPlane,FPlane)

    return Pgrid
end


"""
    construct_plane(χ,X,z,n)
Constructs the forward and backward planes for a single solution plane
"""
function construct_plane(χ::Function,X::AbstractArray{Vector{T}},z,n;periods=1) where T

    plane = zeros(T,(2,prod(n)))

    function prob_fn(prob,i,repeat)
        remake(prob,u0=X[i])
    end
    P = ODEProblem(χ,X[1],(T(0),T(periods)*z))
    EP = EnsembleProblem(P,prob_func=prob_fn)

    sim = solve(EP,Tsit5(),EnsembleSerial(),trajectories=prod(n),save_on=false,save_end=true)
    
    for i = 1:prod(n)
        plane[:,i] = sim.u[i].u[2]
    end

    return plane
end


function postprocess_plane!(X,xbound,ybound,xmode,ymode)
    @views out_of_bounds!(X[1,:],xbound,xmode)
    @views out_of_bounds!(X[2,:],ybound,ymode)
end
"""
    out_of_bounds!(X,boundx,boundy)
Move out of bounds points to the boundary
"""
@views function out_of_bounds!(X,bound,mode)
    
    if mode == :stop
        for i in eachindex(X)
            if X[i] ≤ bound[1]
                X[i] = bound[1]
            elseif X[i] ≥ bound[2]
                X[i] = bound[2]
            end
        end
    elseif mode == :period
        for i in eachindex(X)
            X[i] = rem2pi(X[i],RoundNearest)
        end
    end
end