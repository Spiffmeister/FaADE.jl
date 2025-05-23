

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
function construct_grid(χ::Function,grid::Grid2D{T},z::Vector{T};xmode=:ignore,ymode=:ignore,interpmode=:bicubic,remap=nothing,show_config=true,gridoptions=Dict()) where T

    modelist = [:stop,:period,:ignore]
    xmode ∈ modelist ? nothing : error("mode unavailable")
    ymode ∈ modelist ? nothing : error("mode unavailable")

    if "xbound" ∈ keys(gridoptions)
        xbounds = gridoptions["xbound"]
    else
        xbounds = [minimum(grid.gridx),maximum(grid.gridx)]
    end
    if "ybound" ∈ keys(gridoptions)
        ybounds = gridoptions["ybound"]
    else
        ybounds = [minimum(grid.gridy),maximum(grid.gridy)]
    end

    inversemap = nothing

    if typeof(grid) <: Grid2D
        if typeof(remap) <: Nothing
            xy = [collect(grid[I]) for I in eachindex(grid)]
        elseif typeof(remap) <: Function
            xy = [remap(grid[I]...) for I in eachindex(grid)]
        elseif typeof(remap) <: Tuple
            xy = [remap[1](grid[I]...) for I in eachindex(grid)]
            inversemap = remap[2]
        end
    end

    if show_config
        println("Constructing grid")
        println("Bounds in x and y are being adjusted in mode: ", xmode, ymode)
        if !isnothing(remap)
            println("Incoming grid points remapped to logical coordinates.")
        end
        if !isnothing(inversemap)
            println("Outgoing forward and backward points will be remapped")
        end
        println("Grid will be remapped for ",interpmode," interpolation mode")
    end

    BPlane = construct_plane(χ,xy,z[1],size(grid))
    FPlane = construct_plane(χ,xy,z[2],size(grid))
    # plane = reshape(xy,grid.nx,grid.ny)


    postprocess_plane!(BPlane,xbounds,ybounds,xmode,ymode,inversemap)
    postprocess_plane!(FPlane,xbounds,ybounds,xmode,ymode,inversemap)

    if interpmode == :nearest
        ix,iy = _remap_to_nearest_neighbours(grid,BPlane)
        BPlane = ParGrid{Int,typeof(ix)}(ix,iy,ones(Int,size(ix)))
        
        ix,iy = _remap_to_nearest_neighbours(grid,FPlane)
        FPlane = ParGrid{Int,typeof(ix)}(ix,iy,ones(Int,size(ix)))
    elseif interpmode == :bilinear
        Bplane = _remap_to_linear(grid,BPlane)
        Fplane = _remap_to_linear(grid,FPlane)

        Pgrid = ParallelGrid{eltype(Bplane.weight11),2,typeof(Bplane),typeof(Bplane.weight11)}(Bplane,Fplane)
    # elseif interpmode == :chs
        # Bplane = _remap_to_chs(grid,BPlane)
    else
        Pgrid = ParallelGrid{eltype(BPlane.x),2,typeof(BPlane),typeof(BPlane.x)}(BPlane,FPlane)
    end


    return Pgrid
end
"""
    construct_grid(χ::Function,grid::GridMultiBlock,z::Vector;interpmode=:nearest)
Constructs the backwards and forward planes for a multiblock grid. Returns a dictionary of PGrid objects corresponding to the grids in GridMultiBlock.

By default will return points for nearest neighbour interpolation.
"""
function construct_grid(χ::Function,grid::GridMultiBlock{TT,DIM,MET},z::Vector{TT};gridoptions=Dict()) where {TT,DIM,MET}

    show_config = true

    options = keys(gridoptions)

    if "coords" ∈ options
        if typeof(gridoptions["coords"][1]) <: Function
            coordinate_map = gridoptions["coords"][1]
            inverse_coordinate_map = gridoptions["coords"][2]
        elseif typeof(gridoptions["coords"][1]) <: Symbol
            if gridoptions["coords"] == (:rθ,:xy)
                coordinate_map = (x,y) -> [sqrt(x^2 + y^2), atan(y,x)]
                inverse_coordinate_map = (r,θ) -> [r*cos(θ), r*sin(θ)]
            end
        end
    else
        coordinate_map = nothing
        inverse_coordinate_map = nothing
    end

    "xbound" ∈ options ? xbounds = gridoptions["xbound"] : xbounds = [-Inf,Inf]
    "ybound" ∈ options ? ybounds = gridoptions["ybound"] : ybounds = [-Inf,Inf]

    remapping = :none
    if "remapping" ∈ options
        remapping = gridoptions["remapping"]
        if remapping ∉ [:nearest,:bilinear,:idw,:none]
            @warn "Options for remapping grid is :nearest, :bilinear, :idw (inverse distance weighting), :none. Since it was not in the list it is being set to :none"
            remapping = :none
        else
            println("Remapping grid to weights for $remapping")
        end
    end

    xmode = :ignore
    ymode = :ignore
    if "xmode" ∈ options
        xmode = gridoptions["xmode"]
        if xmode ∉ [:stop,:period,:ignore]
            error("xmode not available")
        end
    end
    if "ymode" ∈ options
        ymode = gridoptions["ymode"]
        if ymode ∉ [:stop,:period,:ignore]
            error("ymode not available")
        end
    end

    # Print the incoming grid config to make it more obvious when things go wrong...
    if show_config
        println("Constructing grid")
        println("Bounds in x and y are being adjusted in mode: ", xmode,", ", ymode)
        if !(remapping == :none)
            println("Incoming grid points remapped to logical coordinates.")
        end
        # if !isnothing(inversemap)
        #     println("Outgoing forward and backward points will be remapped")
        # end
        # println("Grid will be remapped for ",interpmode," interpolation mode")
    end

    PGridStorage = Dict()

    if remapping == :none
        subgrid_interior = [reduce(vcat, [subgrid[i,j] for i in 1:subgrid.nx, j in 1:subgrid.ny]) for subgrid in grid.Grids]
        subgrid_boundaries = [GetBoundaryCoordinates(subgrid,:Circular) for subgrid in grid.Grids]
        subgrid_boundary_nodes = [_boundary_index(subgrid_interior[I],subgrid_boundaries[I],grid.Grids[I]) for I in eachindex(grid.Grids)]
        grid_triangulation = [ triangulate(subgrid_interior[I]; boundary_nodes=subgrid_boundary_nodes[I] ) for I in eachindex(grid.Grids) ]
    end

    for I in eachgrid(grid)
        # Use the single grid version for the subgrid
        Pgrid = construct_grid(χ,grid.Grids[I],z,xmode=:ignore,ymode=:ignore,remap=coordinate_map,show_config=false)

        # Post processing to ensure grid points are inbounds etc...
        postprocess_plane!(Pgrid.Bplane,xbounds,ybounds,xmode,ymode,inverse_coordinate_map)
        postprocess_plane!(Pgrid.Fplane,xbounds,ybounds,xmode,ymode,inverse_coordinate_map)
        
        # If we are supposed to remap or not
        if remapping == :nearest
            Bplane = _remap_to_nearest_neighbours(grid,Pgrid.Bplane)
            Fplane = _remap_to_nearest_neighbours(grid,Pgrid.Fplane)
            PGridStorage[I] = ParallelGrid{eltype(Bplane.weight11),DIM,typeof(Bplane),typeof(Bplane.weight11)}(Bplane,Fplane)

        elseif remapping == :bilinear
            Bplane = _remap_to_linear(grid,Pgrid.Bplane,(xbounds,ybounds),coordinate_map)
            Fplane = _remap_to_linear(grid,Pgrid.Fplane,(xbounds,ybounds),coordinate_map)
            PGridStorage[I] = ParallelGrid{eltype(Bplane.weight11),DIM,typeof(Bplane),typeof(Bplane.weight11)}(Bplane,Fplane)
        
        elseif remapping == :idw
            Bplane = _remap_to_idw(grid,Pgrid.Bplane)
            Fplane = _remap_to_idw(grid,Pgrid.Fplane)
            PGridStorage[I] = ParallelGrid{eltype(Bplane.weight11),DIM,typeof(Bplane),typeof(Bplane.weight11)}(Bplane,Fplane)
            
        elseif remapping == :none
            sgi = _subgrid_index(grid,Pgrid.Bplane,grid_triangulation)
            Bplane = ParGrid{TT,typeof(Pgrid.Bplane.x)}(Pgrid.Bplane.x,Pgrid.Bplane.y,sgi)
            sgi = _subgrid_index(grid,Pgrid.Fplane,grid_triangulation)
            Fplane = ParGrid{TT,typeof(Pgrid.Fplane.x)}(Pgrid.Fplane.x,Pgrid.Fplane.y,sgi)
            PGridStorage[I] = ParallelGrid{eltype(Bplane.x),DIM,typeof(Bplane),typeof(Bplane.x)}(Bplane,Fplane)

        end

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



# function remap_parallel_grid(grid,parallelgrid,method)
# end
"""
    _remap_to_linear
"""
function _remap_to_linear end
"""
"""
function _remap_to_linear(grid::Grid2D{TT,MET},plane::ParGrid) where {TT,MET}
    i = 0; j = 0;

    iindex = zeros(Int,size(plane.x))
    jindex = zeros(Int,size(plane.x))

    weight11 = zeros(TT,size(plane.x))
    weight12 = zeros(TT,size(plane.x))
    weight21 = zeros(TT,size(plane.x))
    weight22 = zeros(TT,size(plane.x))

    for I in eachindex(iindex)
        pt = (plane.x[I], plane.y[I])

        i,j = findcell(grid,pt)
        if i == j == -1
            firstnear,firstnearind = nearestpoint(grid,pt,:linear)

            cartind = CartesianIndices(grid.gridx)[firstnearind]
            i = cartind[1]
            j = cartind[2]
            pt = firstnear
            i,j = findcell(grid,pt)

            @warn "Point is not in the grid, moving the point and recomputing weights."
        end

        iindex[I] = i; jindex[I] = j

        a = grid[i,j]
        b = grid[i+1,j]
        c = grid[i,j+1]
        d = grid[i+1,j+1]

        u,v = _inverse_bilinear_interpolation(a,b,c,d,pt)

        if isnan(u) || isnan(v)
            @show i,j, pt
        end

        weight11[I] = (1-v)*(1-u)
        weight12[I] = u*(1-v)
        weight21[I] = (1-u)*v
        weight22[I] = u*v

    end
    return ParGridLinear{eltype(weight11),typeof(weight11),:Bilinear}(weight11,weight12,weight21,weight22,iindex,jindex,ones(Int,size(iindex)))
end
"""
Remap grid to bilinear interpolation weights.
"""
function _remap_to_linear(grid::GridMultiBlock{TT,2,CartesianMetric},plane::ParGrid) where {TT}
    iindex = zeros(Int,size(plane.x))
    jindex = zeros(Int,size(plane.x))
    sgi = zeros(Int,size(plane.x))

    weight11 = zeros(TT,size(plane.x))
    weight12 = zeros(TT,size(plane.x))
    weight21 = zeros(TT,size(plane.x))
    weight22 = zeros(TT,size(plane.x))


    grid_triangulation = [ triangulate(hcat(subgrid.gridx[:],subgrid.gridy[:])') for subgrid in grid.Grids ]


    for I in eachindex(sgi)
        pt = (plane.x[I], plane.y[I])

        sgi[I] = _findgrid(grid_triangulation,pt)

        subgrid = grid.Grids[sgi[I]]

        # First we need to find the cell the node belongs to
        i,j = findcell(subgrid,pt)

        iindex[I] = i; jindex[I] = j

        ΔxΔy = (subgrid.gridx[i+1,j] - subgrid.gridx[i,j])*(subgrid.gridy[i,j+1] - subgrid.gridy[i,j])
        
        weight11[I] = (subgrid.gridx[i+1,j] - pt[1]) * (subgrid.gridy[i,j+1] - pt[2])/ΔxΔy  # bottom left
        weight12[I] = (pt[1] - subgrid.gridx[i,j]) * (subgrid.gridy[i,j+1] - pt[2])/ΔxΔy    # bottom right
        weight21[I] = (subgrid.gridx[i+1,j] - pt[1])   * (pt[2] - subgrid.gridy[i,j])/ΔxΔy    # top left
        weight22[I] = (pt[1] - subgrid.gridx[i,j])   * (pt[2] - subgrid.gridy[i,j])/ΔxΔy    # top right

    end

    return ParGridLinear{eltype(weight11),typeof(weight11),:Bilinear}(weight11,weight12,weight21,weight22,iindex,jindex,sgi)
end
"""
    _remap_to_linear(grid::GridMultiBlock,plane::ParGrid)
Remap the grid to its bilinear interpolation weights.

WARNING: This function may move points if they lie outside the domain.
"""
function _remap_to_linear(grid::GridMultiBlock{TT,2,CurvilinearMetric},plane::ParGrid,bounds,mapping) where {TT}
    iindex = zeros(Int,size(plane.x))
    jindex = zeros(Int,size(plane.x))
    subgridindex = zeros(Int,size(plane.x))

    weight11 = zeros(TT,size(plane.x))
    weight12 = zeros(TT,size(plane.x))
    weight21 = zeros(TT,size(plane.x))
    weight22 = zeros(TT,size(plane.x))

    i = 0; j = 0; sgi = 0;

    nearpt = (0,0)

    for I in eachindex(subgridindex)
        pt = (plane.x[I], plane.y[I])

        firstpt, ind, sgi = nearestpoint(grid,pt,:cartesian)

        #TODO : REMOVE ALL THESE TRY CATCH
        # We should check this is the correct grid and index
        if (ind[1] == 1) || (ind[1] == grid.Grids[sgi].nx) || (ind[2] == 1) || (ind[2] == grid.Grids[sgi].ny)
            
            tmpi,tmpj = findcell(grid.Grids[sgi],pt)
            # if -1 then the point is not in the domain and we need to correct it
            if tmpi == tmpj == -1

                joints = grid.Joint[sgi]
                for joint in joints # this will correct most instances
                    i,j = findcell(grid.Grids[joint.index],pt)
                    sgi = joint.index
                    if i != -1
                        break
                    end
                end

                # we need to try again since it may be outside the domain
                # if it is, move it to the nearest point
                i,j = findcell(grid.Grids[sgi],pt)
                if i == j == -1
                    pt = firstpt
                end
                
            end
            
        end

        i,j = findcell(grid.Grids[sgi],pt)
        
        subgrid = grid.Grids[sgi]

        iindex[I] = i; jindex[I] = j; subgridindex[I] = sgi;



        a = subgrid[i,j]
        b = subgrid[i+1,j]
        c = subgrid[i,j+1]
        d = subgrid[i+1,j+1]

        u,v = _inverse_bilinear_interpolation(a,b,c,d,pt)

        weight11[I] = (1-v)*(1-u)
        weight12[I] = u*(1-v)
        weight21[I] = (1-u)*v
        weight22[I] = u*v

    end

    return ParGridLinear{eltype(weight11),typeof(weight11),:Bilinear}(weight11,weight12,weight21,weight22,iindex,jindex,subgridindex)
end



"""
    _inverse_bilinear_interpolation
Given the four points of a bilinear interpolation, find the weights for the point pt

`` S(u,v) = uvB + u(1-v)D + v(1-u)C + (1-u)(1-v)A ``

where `A`,`B`,`C`,`D` are the four points of the quadrilateral;
C -- D
|    |
A -- B

"""
function _inverse_bilinear_interpolation(a::Tuple{TT,TT},b,c,d,q) where TT

    Q0 = - (b[2] - q[2]) * a[1] + (-q[2] + a[2]) * b[1] - q[1] * a[2] + b[2] * q[1]
    Q1 = (-(q[2] - 2 * b[2] + d[2]) * a[1] + (q[2] - 2 * a[2] + c[2]) * b[1] - (-d[1] - q[1]) * a[2] + (-c[1] - q[1]) * b[2] - (c[2] - d[2]) * q[1] - (-c[1] + d[1]) * q[2])
    Q2 = (-(b[2] - d[2]) * a[1] + (a[2] - c[2]) * b[1] + c[1] * b[2] - c[1] * d[2] - d[1] * a[2] + d[1] * c[2])

   

    if abs(Q2) > 1e-12
        vroot = ((-Q1 + sqrt(Q1^2 - 4*Q0*Q2))/(2*Q2), (-Q1 - sqrt(Q1^2 - 4*Q0*Q2))/(2*Q2))
    else
        vroot = 2Q0 ./ ( -Q1 + sqrt(Q1^2 - 4*Q0*Q2)  , -Q1 - sqrt(Q1^2 - 4*Q0*Q2) )
    end

    if 0 ≤ vroot[1] ≤ 1
        v = vroot[1]
    elseif 0 ≤ vroot[2] ≤ 1
        v = vroot[2]
    elseif isapprox(vroot[1],TT(1),atol=1e-10) || isapprox(vroot[2],TT(1),atol=1e-10)
        v = TT(1)
    elseif isapprox(vroot[1],TT(0),atol=1e-10) || isapprox(vroot[2],TT(0),atol=1e-10)
        v = TT(0)
    else
        error("Cannot find node.")
    end


    if !iszero(v)
        u = (q[1] + a[1]*v - c[1]*v - a[1]) / ((a[1] - b[1] - c[1] + d[1])*v - a[1] + b[1])
    else
        q[1] == a[1] ? u = TT(0) : u = TT(1)
    end

    return u,v
end


"""
    _construct_chs
Constructs the `CubicHermiteSpline` object for a single volume, does not perform point location
"""
function _remap_to_chs(grid::Grid2D{TT,MET},plane::ParGrid) where {TT,MET}
    z = zeros(TT,length(plane.x))
    dzdx = zeros(TT,length(plane.x))
    dzdy = zeros(TT,length(plane.x))



    chsinterp = BivariateCHSInterpolation(grid.gridx[:],grid.gridy[:],z,dzdx,dzdy)

    chsret = ParGridCHS{typeof(chsinterp)}(chsinterp,plane.x,plane.y,ones(Int,(1,1)))

    return chsret
end



"""
    _remap_to_idw(grid::GridMultiBlock,plane::ParGrid)
Remap the interpolation points to their inverse distance weights
"""
function _remap_to_idw(grid::GridMultiBlock{TT,2},plane::ParGrid) where {TT}
    sgi = zeros(Int,size(plane.x))

    weight11 = zeros(TT,size(plane.x))
    weight12 = zeros(TT,size(plane.x))
    weight21 = zeros(TT,size(plane.x))
    weight22 = zeros(TT,size(plane.x))

    iindex = zeros(Int,size(plane.x))
    jindex = zeros(Int,size(plane.x))


    for I in eachindex(sgi)
        i = 0; j = 0;

        pt = (plane.x[I], plane.y[I])

        sgi[I] = findgrid(grid,pt)
        
        subgrid = grid.Grids[sgi[I]]

        cartind = CartesianIndices(subgrid.gridx)

        # First we need to find the cell the node belongs to
        i,j = findcell(subgrid,pt)
        if i == j == -1
            firstnear,firstnearind = nearestpoint(subgrid,pt,:linear)

            cartind = CartesianIndices(subgrid.gridx)[firstnearind]
            i = cartind[1]
            j = cartind[2]
            pt = firstnear
            i,j = findcell(subgrid,pt)

            @warn "Point is not in the grid, moving the point and recomputing weights."
        end

        iindex[I] = i
        jindex[I] = j


        p1 = subgrid[i,j]
        p2 = subgrid[i+1,j]
        p3 = subgrid[i,j+1]
        p4 = subgrid[i+1,j+1]

        dist = 1/norm(p1 .- pt)^2 + 1/norm(p2 .- pt)^2 + 1/norm(p3 .- pt)^2 + 1/norm(p4 .- pt)^2

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

    w11 = zeros(Float64,1,1)

    return ParGridLinear{Float64,Matrix{Float64},:NearestNeighbour}(w11,w11,w11,w11,ix,iy,sgi)
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
function _subgrid_index(grid::GridMultiBlock,plane::ParGrid,grid_triangulation)
    containedgrid = zeros(Int,size(plane.x))

    # Form a DelaunayTriangulation of the grid, since this is the same method used in CubicHermiteSpline it should give the correct indices
    # subgrid_interior = [reduce(vcat, [subgrid[i,j] for i in 1:subgrid.nx, j in 1:subgrid.ny]) for subgrid in grid.Grids]
    # subgrid_boundaries = [GetBoundaryCoordinates(subgrid,:Circular) for subgrid in grid.Grids]
    # subgrid_boundary_nodes = [_boundary_index(subgrid_interior[I],subgrid_boundaries[I],grid.Grids[I]) for I in eachindex(grid.Grids)]
    # grid_triangulation = [ triangulate(subgrid_interior[I]; boundary_nodes=subgrid_boundary_nodes[I] ) for I in eachindex(grid.Grids) ]
    
    # Loop over all grid points and use the DelaunayTriangulation to locate the grid
    for J in eachindex(plane.x)
        pt = (plane.x[J],plane.y[J])
        gridind = _findgrid(grid_triangulation,pt)
        if gridind == 0
            nearpt, _, gridind = nearestpoint(grid,pt)
            plane.x[J] = nearpt[1]
            plane.y[J] = nearpt[2]
        end
        containedgrid[J] = gridind

    end
    return containedgrid
end





"""
    postprocess_plane!(X,xbound,ybound,xmode,ymode, remap=(:xy,:xy))
"""
function postprocess_plane!(X,xbound,ybound,xmode,ymode, remap=nothing)
    @views out_of_bounds!(X.x,xbound,xmode)
    @views out_of_bounds!(X.y,ybound,ymode)

    if remap == (:rθ,:xy)
    elseif typeof(remap) <: Function
        for I in eachindex(X.x)
            X.x[I], X.y[I] = remap(X.x[I],X.y[I])
        end
    end

end

"""

"""
function _remap_coordinates_rθ_to_xy(plane)
    for I in eachindex(plane.x)
        x = plane.x[I] * cos(plane.y[I])
        y = plane.x[I] * sin(plane.y[I])
        plane.x[I] = x
        plane.y[I] = y
    end
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


"""
Find the subgrid index using the provided triangulation from `DelaunayTriangulation`.
"""
function _findgrid(triangulated_grids::Vector{TRI},pt::Tuple{TT,TT}) where {TT,TRI <: Triangulation}

    sgi = 0
    for (gridindex,trigrid) in enumerate(triangulated_grids)
        point_tri = find_triangle(trigrid,pt)
        tmpid = all(0 .≤ point_tri)
        if tmpid
            sgi = gridindex
            continue
        end
    end

    return sgi

end

"""
    _boundary_index(nodes,boundary,grid)
Form a vector of the linear indices of the boundary nodes in `nodes` using the points in `boundary`.

This avoids using the methods from `DelaunayTriangulation` which may not be as efficient.
"""
function _boundary_index(nodes,boundary,grid)
    boundary_indices = zeros(Int64,2sum(size(grid))-3)
    for (I,node) in enumerate(boundary)
        for J in eachindex(nodes)
            if node == nodes[J]
                boundary_indices[I] = J
                continue
            end
        end
        boundary_indices[end] = boundary_indices[1]
    end
    return boundary_indices    
end


