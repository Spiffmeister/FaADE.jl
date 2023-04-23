

struct ParallelGrid
    Plane   :: AbstractArray
    Fplane  :: AbstractArray
    Bplane  :: AbstractArray
end


"""
    construct_grid(χ::Function,grid::Grid2D,z::Vector)
Constructs the backwards and forward planes for a given plane
"""
function construct_grid(χ::Function,grid::Grid2D,z::Vector)

    if typeof(grid) <: Grid2D
        xy = [[x,y] for x in grid.gridx for y in grid.gridy]
    end
    # χₘₙ = 2.1e-3 + 5.0e-3
    # params = (ϵₘₙ=[χₘₙ/2., χₘₙ/3.],m=[2.0, 3.0],n=[1.0, 2.0])
    BPlane = construct_plane(χ,xy,z[1],(grid.nx,grid.ny))
    FPlane = construct_plane(χ,xy,z[2],(grid.nx,grid.ny))

    Pgrid = ParallelGrid(hcat(xy...),FPlane,BPlane)

    return Pgrid
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

    sim = solve(EP,Tsit5(),EnsembleSerial(),trajectories=prod(n),save_on=false,save_end=true)
    for i = 1:length(sim.u)
        plane[:,i] = sim.u[i].u[2]
    end
    return plane
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