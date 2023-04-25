

struct plane
    x   :: AbstractArray
    y   :: AbstractArray
end

struct ParallelGrid
    X   :: AbstractArray
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
    # χₘₙ = 2.1e-3 + 5.0e-3
    # params = (ϵₘₙ=[χₘₙ/2., χₘₙ/3.],m=[2.0, 3.0],n=[1.0, 2.0])
    BPlane = construct_plane(χ,xy,z[1],(grid.nx,grid.ny))
    FPlane = construct_plane(χ,xy,z[2],(grid.nx,grid.ny))

    # Pgrid = ParallelGrid(FPlane,BPlane)

    return BPlane
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

    sim = solve(EP,Tsit5(),EnsembleDistributed(),trajectories=prod(n),save_on=false,save_end=true)
    
    for i = 1:length(sim.u)
        plane[:,i] = sim.u[i].u[:]
    end

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