


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
"""
function construct_grid(χ::Function,grid::Grid2D,z::Vector)

    if typeof(grid) <: Grid2D
        xy = [[x,y] for x in grid.gridx for y in grid.gridy]
    end

    BPlane = construct_plane(χ,xy,z[1],(grid.nx,grid.ny))
    FPlane = construct_plane(χ,xy,z[2],(grid.nx,grid.ny))

    # Pgrid = ParallelGrid(FPlane,BPlane)

    return BPlane
end


"""
    construct_plane(χ,X,z,n)
Constructs a single z plane
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


"""
    out_of_bounds!(X,boundx,boundy)
Move out of bounds points to the boundary
"""
function out_of_bounds!(X,boundx,boundy)
    
    for i in size(X[1])
        if X[i,1] ≤ boundx[1]
            X[i,1] = boundx[1]
        elseif X[1,i] ≥ boundx[2]
            X[i,1] = boundx[2]
        end

        X[:,2]
    end
    X
end