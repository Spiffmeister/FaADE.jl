

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
function construct_grid(χ::Function,grid::Grid2D{T},z::Vector{T};xmode=:stop,ymode=:period,interpmode=:bicubic) where T

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

    if interpmode == :nearest
        ix,iy = _remap_to_nearest_neighbours(grid,BPlane)
        BPlane = ParGrid{Int,typeof(ix)}(ix,iy,ones(Int,size(ix)))
        
        ix,iy = _remap_to_nearest_neighbours(grid,FPlane)
        FPlane = ParGrid{Int,typeof(ix)}(ix,iy,ones(Int,size(ix)))
    end

    Pgrid = ParallelGrid{eltype(BPlane.x),2,typeof(BPlane),typeof(BPlane.x)}(BPlane,FPlane)

    return Pgrid
end
"""
    construct_grid(χ::Function,grid::GridMultiBlock,z::Vector;interpmode=:nearest)
Constructs the backwards and forward planes for a multiblock grid. Returns a dictionary of PGrid objects corresponding to the grids in GridMultiBlock.

By default will return points for nearest neighbour interpolation.
"""
function construct_grid(χ::Function,grid::GridMultiBlock{TT,DIM},z::Vector{TT};xmode=:stop,ymode=:period,interpmode=:nearest) where {TT,DIM}

    PGridStorage = Dict()

    minx = minimum([minimum(grid.Grids[I].gridx) for I in eachgrid(grid)])
    maxx = maximum([maximum(grid.Grids[I].gridx) for I in eachgrid(grid)])
    miny = minimum([minimum(grid.Grids[I].gridy) for I in eachgrid(grid)])
    maxy = maximum([maximum(grid.Grids[I].gridy) for I in eachgrid(grid)])
    
    for I in eachgrid(grid)
        Pgrid = construct_grid(χ,grid.Grids[I],z,xmode=:ignore,ymode=:ignore)
        
        if interpmode == :nearest
            ix,iy,sgi = _remap_to_nearest_neighbours(grid,Pgrid.Bplane)
            Bplane = ParGrid{Int,typeof(ix)}(ix,iy,sgi)
        elseif interpmode == :linear
            postprocess_plane!(Pgrid.Bplane,[minx,maxx],[miny,maxy],xmode,ymode)
            Bplane = _remap_to_linear(grid,Pgrid.Bplane)
        elseif interpmode == :idw
            postprocess_plane!(Pgrid.Bplane,[minx,maxx],[miny,maxy],xmode,ymode)
            Bplane = _remap_to_idw(grid,Pgrid.Bplane)
        else
            postprocess_plane!(Pgrid.Bplane,[minx,maxx],[miny,maxy],xmode,ymode)
            sgi = _subgrid_index(grid,Pgrid.Bplane)
            Bplane = ParGrid{TT,typeof(Pgrid.Bplane.x)}(Pgrid.Bplane.x,Pgrid.Bplane.y,sgi)
        end
        
        if interpmode == :nearest
            ix,iy,sgi = _remap_to_nearest_neighbours(grid,Pgrid.Fplane)
            Fplane = ParGrid{Int,typeof(ix)}(ix,iy,sgi)
        elseif interpmode == :linear
            postprocess_plane!(Pgrid.Fplane,[minx,maxx],[miny,maxy],xmode,ymode)
            Fplane = _remap_to_linear(grid,Pgrid.Fplane)
        elseif interpmode == :idw
            postprocess_plane!(Pgrid.Fplane,[minx,maxx],[miny,maxy],xmode,ymode)
            Fplane = _remap_to_idw(grid,Pgrid.Fplane)
        else
            postprocess_plane!(Pgrid.Fplane,[minx,maxx],[miny,maxy],xmode,ymode)
            sgi = _subgrid_index(grid,Pgrid.Fplane)
            Fplane = ParGrid{TT,typeof(Pgrid.Fplane.x)}(Pgrid.Fplane.x,Pgrid.Fplane.y,sgi)
        end

        # PGridStorage[I] = ParallelGrid{eltype(Bplane.x),DIM,typeof(Bplane),typeof(Bplane.x)}(Bplane,Fplane)

        PGridStorage[I] = ParallelGrid{eltype(Bplane.weight11),DIM,typeof(Bplane),typeof(Bplane.weight11)}(Bplane,Fplane)
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



function _remap_to_linear(grid::GridMultiBlock,plane::ParGrid)
    index = zeros(Int,size(plane.x))
    sgi = zeros(Int,size(plane.x))

    weight11 = zeros(Int,size(plane.x))
    weight12 = zeros(Int,size(plane.x))
    weight21 = zeros(Int,size(plane.x))
    weight22 = zeros(Int,size(plane.x))


    for I in eachindex(sgi)
        @show pt = (plane.x[I], plane.y[I])

        # First we need to find the cell the node belongs to
        @show ix,iy,sgi[I] = findcell(grid,pt)

        i = ix
        j = iy

                
        @show (grid.Grids[sgi[I]].gridx[ix[I]+1] - grid.Grids[sgi[I]].gridx[ix])
        @show (grid.Grids[sgi[I]].gridy[iy[I]+1] - grid.Grids[sgi[I]].gridy[iy])
        
        subgrid = grid.Grids[sgi[I]]

        # Now we need to perform inverse bilinear interpolation to find the weights
        k0 = (subgrid.gridx[ix,iy] - subgrid.gridx[ix+1,iy]) * (pt[2] - subgrid.gridy[ix,iy])

        k1 = subgrid.gridx[ix,iy]*(-pt[2] + subgrid.gridy[ix,iy] - subgrid.gridy[ix+1,iy] + subgrid.gridy[ix,iy+1]) * 
            subgrid.gridy[ix,iy] * subgrid.gridx[ix+1,iy+1] + subgrid.gridx[ix+1,iy] * (pt[2] - subgrid.gridy[ix,iy+1]) +
            pt[2] * (subgrid.gridx[ix+1,iy+1] - subgrid.gridx[ix,iy+1]) + subgrid.gridx[ix,iy+1] * subgrid.gridy[ix+1,iy]
        
        k2 = subgrid.gridx[ix,iy] * (subgrid.gridy[ix+1,iy] - 2*subgrid.gridy[ix,iy+1] + subgrid.gridy[ix+1,iy+1]) + 
            subgrid.gridx[ix,iy] * (-subgrid.gridx[ix+1,iy] + 2*subgrid.gridx[ix,iy+1] - subgrid.gridx[ix+1,iy+1]) +
            subgrid.gridx[ix+1,iy] * subgrid.gridy[ix,iy+1] * subgrid.gridx[ix,iy+1] * (-subgrid.gridy[ix+1,iy] - subgrid.gridy[ix+1,iy+1]) + subgrid.gridx[ix+1,iy+1] * subgrid.gridy[ix,iy+1]


        # distance between the curves p_11 - p_12 and p_21 - p_22
        v = max( (-k1+sqrt(k1^2 - 4*k0*k2))/2k0, (-k1-sqrt(k1^2 - 4*k0*k2))/2k0 )

        u = - v*(subgrid.gridx[i,j] - subgrid.gridx[i,j+1]) / (v*(subgrid.gridx[i,j] - subgrid.gridx[i+1,j] + subgrid.gridx[i,j+1] - subgrid.gridx[i+1,j+1]) - subgrid.gridx[i,j] + subgrid.gridx[i+1,j])

        # TODO: we should check if the point is correct



        ΔxΔy = (subgrid[ix+1,iy] - subgrid[ix,iy])*(subgrid[ix,iy+1] - subgrid[ix,iy])
        
        weight11[I] = (subgrid.gridx[ix+1,iy] - pt[1]) * (subgrid.y[ix,iy+1] - pt[2])/ΔxΔy  # bottom left
        weight12[I] = (subgrid.gridx[ix+1,iy] - pt[1]) * (pt[2] - subgrid.y[ix,iy])/ΔxΔy    # bottom right
        weight21[I] = (pt[1] - subgrid.gridy[ix,iy])   * (subgrid.y[ix,iy] - pt[2])/ΔxΔy    # top left
        weight22[I] = (pt[1] - subgrid.gridy[ix,iy])   * (pt[2] - subgrid.y[ix,iy])/ΔxΔy    # top right

        
        J = ix[I] + grid.Grids[sgi[I]].ny*(iy[I]-1) # linearise index
        index[I] = J
    end

    return ParGridLinear{eltype(weightx),typeof(weightx),:Bilinear}(weight11,weight12,weight21,weight22,index,sgi)
end

"""
    _remap_to_idw(grid::GridMultiBlock,plane::ParGrid)
Remap the interpolation points to their inverse distance weights
"""
function _remap_to_idw(grid::GridMultiBlock{TT},plane::ParGrid) where TT
    index = zeros(Int,size(plane.x))
    sgi = zeros(Int,size(plane.x))

    weight11 = zeros(TT,size(plane.x))
    weight12 = zeros(TT,size(plane.x))
    weight21 = zeros(TT,size(plane.x))
    weight22 = zeros(TT,size(plane.x))

    iindex = zeros(Int,size(plane.x))
    jindex = zeros(Int,size(plane.x))


    for I in eachindex(sgi)
        @show pt = (plane.x[I], plane.y[I])

        # First we need to find the cell the node belongs to
        @show i,j,sgi[I] = findcell(grid,pt)

        iindex[I] = i
        jindex[I] = j

        subgrid = grid.Grids[sgi[I]]

        p1 = subgrid[i,j]
        p2 = subgrid[i+1,j]
        p3 = subgrid[i,j+1]
        p4 = subgrid[i+1,j+1]

        dist = 1/norm(p1 .- pt)^2 + 1/norm(p2 .- pt)^2 + 1/norm(p3 .- pt)^2 + 1/norm(p4 .- pt)^2

        # @show I, pt, dist
        # @show 1/norm(p1 .- pt)/dist

        weight11[I] = TT(1)/norm(p1 .- pt)^2 /dist
        weight12[I] = TT(1)/norm(p2 .- pt)^2 /dist
        weight21[I] = TT(1)/norm(p3 .- pt)^2 /dist
        weight22[I] = TT(1)/norm(p4 .- pt)^2 /dist

    end



    return ParGridLinear{TT,typeof(weight11),:InverseDistanceWeighting}(weight11,weight12,weight21,weight22,iindex,jindex,sgi)
end

"""
    _remap_to_nearest_neighbours(grid::GridMultiBlock,plane::ParGrid)
Remap the points generated by an x,y plane to their nearest neighbour in the full grid
"""
function _remap_to_nearest_neighbours(grid::GridMultiBlock,plane::ParGrid)
    ix = zeros(Int,size(plane.x))
    iy = zeros(Int,size(plane.y))
    sgi = zeros(Int,size(plane.x))

    for I in eachindex(ix)
        _,ind,gridind = nearestpoint(grid,(plane.x[I],plane.y[I]))
        ix[I] = ind
        iy[I] = ind
        sgi[I] = gridind
    end
    return ix,iy,sgi
end
"""
    _remap_to_nearest_neighbours(grid::Grid2D,plane::ParGrid)
Remap the points for a single Grid2D to their nearest neighbours indexes
"""
function _remap_to_nearest_neighbours(grid::Grid2D,plane::ParGrid)
    ix = zeros(Int,size(plane.x))
    iy = zeros(Int,size(plane.y))

    for I in eachindex(ix)
        _, ind = nearestpoint(grid,(plane.x[I],plane.y[I]))
        ix[I] = iy[I] = ind
    end
    return ix,iy
end
"""
    _subgrid_index(grid::GridMultiBlock,plane::ParGrid,I::Int)
Find the subgrid index closest to each point in a ParGrid.
"""
function _subgrid_index(grid::GridMultiBlock,plane::ParGrid)
    containedgrid = zeros(Int,size(plane.x))
    for J in eachindex(plane.x)
        gridind = findgrid(grid,(plane.x[J],plane.y[J]),mode=:nearest)
        containedgrid[J] = gridind
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
                end
            end
        else
            # @show X[1], bound[end], X[1] > bound[end]
            for i in eachindex(X)
                X[i] = rem2pi(X[i],RoundDown)
                # if X[i] > bound[end] # if X[i] ~ 2π it can round to > 2π instead of 0
                #     X[i] = bound[end]
                # end
            end
            # @show X[1] == bound[end]
        end
    elseif mode == :ignore
    end
end


Base.show(io::IO, PG::ParallelGrid) = print(io, " generated parallel grid.")

