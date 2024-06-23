

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
function construct_grid(χ::Function,grid::Grid2D{T},z::Vector{T};xmode=:stop,ymode=:period) where T

    modelist = [:stop,:period,:ignore]
    xmode ∈ modelist ? nothing : error("mode unavailable")
    ymode ∈ modelist ? nothing : error("mode unavailable")

    if typeof(grid) <: Grid2D
        xy = [collect(grid[I]) for I in eachindex(grid)]
    end

    BPlane = construct_plane(χ,xy,z[1],size(grid))
    FPlane = construct_plane(χ,xy,z[2],size(grid))
    # plane = reshape(xy,grid.nx,grid.ny)

    postprocess_plane!(BPlane,[grid.gridx[1],grid.gridx[end]],[grid.gridy[1],grid.gridy[end]],xmode,ymode)
    postprocess_plane!(FPlane,[grid.gridx[1],grid.gridx[end]],[grid.gridy[1],grid.gridy[end]],xmode,ymode)

    Pgrid = ParallelGrid{T,2,typeof(BPlane.x)}(BPlane,FPlane)#,zeros(size(BPlane.x)),zeros(size(BPlane.x)))

    return Pgrid
end


function construct_grid(χ::Function,grid::GridMultiBlock{TT,DIM},z::Vector{TT};xmode=:stop,ymode=:preiod) where {TT,DIM}

    PGridStorage = Dict()

    for I in eachgrid(grid)
        Pgrid = construct_grid(χ,grid.Grids[I],z,xmode=xmode,ymode=ymode)

        sgi = _subgrid_index(grid,Pgrid.Bplane,I)
        Bplane = ParGrid{TT,typeof(Pgrid.Bplane.x)}(Pgrid.Bplane.x,Pgrid.Bplane.y,sgi)
        
        sgi = _subgrid_index(grid,Pgrid.Fplane,I)
        Fplane = ParGrid{TT,typeof(Pgrid.Fplane.x)}(Pgrid.Fplane.x,Pgrid.Fplane.y,sgi)

        PGridStorage[I] = ParallelGrid{TT,DIM,typeof(Bplane.x)}(Bplane,Fplane)
    end

    return PGridStorage
end

"""
    construct_plane(χ,X,z,n)
Constructs the forward and backward planes for a single solution plane
"""
function construct_plane(χ::Function,X::AbstractArray{Vector{T}},z,n;periods=1) where T

    function prob_fn(prob,i,repeat)
        remake(prob,u0=X[i])
    end
    P = ODEProblem(χ,X[1],(T(0),T(periods)*z))
    EP = EnsembleProblem(P,prob_func=prob_fn)

    sim = solve(EP,Tsit5(),EnsembleSerial(),trajectories=prod(n),reltol=1e-6,
        save_on=false,save_end=true)

    planex = zeros(T,n)
    planey = zeros(T,n)
    for i = 1:prod(n)
        planex[i] = sim.u[i].u[2][1]
        planey[i] = sim.u[i].u[2][2]
    end

    return ParGrid{T,typeof(planex)}(planex,planey,zeros(Int,1,1))
end

function _subgrid_index(grid::GridMultiBlock,plane::ParGrid,I::Int)
    containedgrid = zeros(Int,size(plane.x))
    for I in eachgrid(grid) # For each point, check which subgrid it falls within 
        minx = minimum(grid.Grids[I].gridx)
        maxx = maximum(grid.Grids[I].gridx)
        miny = minimum(grid.Grids[I].gridy)
        maxy = maximum(grid.Grids[I].gridy)

        for J in eachindex(plane.x)
            if (minx ≤ plane.x[J] ≤ maxx) && (miny ≤ plane.y[J] ≤ maxy)
                containedgrid[J] = I
            end
        end
    end
    return containedgrid
end

function postprocess_plane!(X,xbound,ybound,xmode,ymode)
    @views out_of_bounds!(X.x,xbound,xmode)
    @views out_of_bounds!(X.y,ybound,ymode)
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
        
        
        if bound[end] ≤ 3π/2
            for i in eachindex(X)
                X[i] = rem2pi(X[i],RoundNearest)
                if X[i] ≥ bound[end]
                    X[i] = bound[end]
                    # println(X[i])
                end
            end
        else
            for i in eachindex(X)
                X[i] = rem2pi(X[i],RoundDown)
            end
        end
    end
end


Base.show(io::IO, PG::ParallelGrid) = print(io, " generated parallel grid.")

