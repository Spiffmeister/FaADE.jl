




"""
    nearestpoint
Find the closest point in a grid to a given point
"""
function nearestpoint end
"""
    nearestpoint(grid::Grid2D,pt::Tuple{TT,TT},indextype=:linear)
Takes a 2D grid and finds the nearest point.

Returns the value of the nearest point and the linear index of that point.
"""
function nearestpoint(grid::Grid2D,pt::Tuple{TT,TT},indextype=:linear) where TT

    tmpdist = TT(0)
    dist = TT(1e10)
    linind = 0
    for I in eachindex(grid)
        # Check the distance to each point in the grid and record the smallest distance
        tmpdist = sqrt( (grid.gridx[I] - pt[1])^2 + (grid.gridy[I] - pt[2])^2 )
        if tmpdist < dist
            dist = tmpdist
            linind = I
        end
    end

    if indextype == :linear
        # Return the LinearIndex in the array
        return (grid.gridx[linind],grid.gridy[linind]), linind
    elseif indextype == :cartesian
        # Return the i,j coordinate in the array
        cartind = CartesianIndices(grid.gridx)[linind]
        return (grid.gridx[linind],grid.gridy[linind]), (cartind[1],cartind[2])
    else
        error("Index type is either :linear or :cartesian")
    end
end
"""
    findgrid(grid::GridMultiBlock{TT,DIM,CartesianMetric},pt::Tuple{TT,TT};mode=:inside)
Check which cartesian grid a point is in.
"""
function findgrid(grid::GridMultiBlock{TT},pt::Tuple{TT,TT};mode=:inside) where {TT}
    gridind = 0
    for I in eachgrid(grid)
        minx = grid.Grids[I].gridx[1]
        maxx = grid.Grids[I].gridx[end]
        miny = grid.Grids[I].gridy[1]
        maxy = grid.Grids[I].gridy[end]

        if (minx <= pt[1] <= maxx) & (miny <= pt[2] <= maxy)
            if gridind == 0
                gridind = I
            else
                # its possible a point lands on the boundary of two grids
                # grid.Joint[I]
            end
        end
    end

    if gridind != 0
        # if the point is in a grid
        return gridind
    elseif (gridind == 0) & (mode == :nearest)
        # if the point is not in a grid but we want the nearest grid
        gridind = nearestpoint(grid,pt)[3]
        @warn "Point $(pt) is not in any grid, since nearest grid has been requested, returning $(gridind)."
        return nearestpoint(grid,pt)[3]
    elseif (gridind == 0) & (mode == :inside)
        # if the point is not inside a grid and we wanted it to be then error out
        error("Point $(pt) is not in any grid")
    end
end
"""
    findcell(grid::GridMultiBlock,pt::Tuple{TT,TT})
Find which cell bounds the point
"""
function findcell(grid::GridMultiBlock{TT,DIM,MET},pt::Tuple{TT,TT}) where {TT,DIM,MET}

    ix = 0
    jy = 0
    subgridindex = 0


    # Find the nearest point in the grid
    if MET == CartesianMetric
        ind, subgridindex = findgrid(grid,pt,mode=:nearest)
    else
        _, ind, subgridindex = nearestpoint(grid,pt)
    end

    subgrid = grid.Grids[subgridindex]
    
    # As a reference point we want the indices of the bottom left point of the cell so the point is bounded by
    cartind = CartesianIndices(subgrid.gridx)
    i = cartind[ind][1]
    j = cartind[ind][2]
    @show subgrid[i,j], i, j
    # We need to check if the point is bounded by the grid or not and shift the indices accordingly
    if i+1 > subgrid.nx
        i = i-1
    end
    if j+1 > subgrid.ny
        j = j-1
    end

    #   [i,j] - [i,j+1] - [i+1,j+1] - [i+1,j]
    # point must be close to a1, but unknown which square
    # b0 --b1 --b2
    # |     |   |
    # a0 --a1 --a2
    # |     |   |
    # c0 --c1 --c2 
    # Is the point above or below the line joining the two points a1--a2

    # TODO: this is awful fix this
    inside = _check_inside(subgrid,i,j,pt) #check a1--a2--b2--b1
    if inside
        ix = i
        jy = j
    else
        inside = _check_inside(subgrid,i-1,j,pt) # check a0--a1--b1--b0
        if inside
            ix = i-1
            jy = j
        else
            inside = _check_inside(subgrid,i-1,j-1,pt) # check a0--a1--c1--c0
            if inside
                ix = i-1
                jy = j-1
            else
                inside = _check_inside(subgrid,i,j-1,pt) # check a1--a2--c2--c1
                if inside
                    ix = i
                    jy = j-1
                end
            end
        end
    end

    if ix == iy == 0

    end

    return ix, jy, subgridindex
end

function findcell(grid::Grid2D,pt::Tuple{TT,TT}) where {TT}

    ix = 0
    jy = 0
    subgridindex = 0

    @show _, inds = nearestpoint(grid,pt,:cartesian)
    i = inds[1]
    j = inds[2]

    if i+1 > grid.nx
        i = i-1
    end
    if j+1 > grid.ny
        j = j-1
    end

    #   [i,j] - [i,j+1] - [i+1,j+1] - [i+1,j]
    # point must be close to a1, but unknown which square
    # b0 --b1 --b2
    # |     |   |
    # a0 --a1 --a2
    # |     |   |
    # c0 --c1 --c2 
    # Is the point above or below the line joining the two points a1--a2

    # TODO: this is awful fix this
    inside = _check_inside(grid,i,j,pt) #check a1--a2--b2--b1
    if inside
        ix = i
        jy = j
    end
    if (!inside) && (i != 1)
        inside = _check_inside(grid,i-1,j,pt) # check a0--a1--b1--b0
        if inside
            ix = i-1
            jy = j
        end
    end
    if (!inside) && (j != 1)
        inside = _check_inside(grid,i,j-1,pt) # check a1--a2--c2--c1
        if inside
            ix = i
            jy = j-1
        end
    end
    if (!inside) && (i != 1) && (j != 1)
        inside = _check_inside(grid,i-1,j-1,pt) # check a0--a1--c1--c0
        if inside
            ix = i-1
            jy = j-1
        end
    end

    if ix == jy == 0
        if pt == (grid.gridx[i,j],grid.gridy[i,j])
            ix = i
            jy = j
        else
            error("Point $(pt) is not in any cell.")
        end
    end

    return ix, jy
end




"""
    _check_inside(grid,i,j,pt)
Check if a point is inside a cell with the point `i,j` at the bottom left.
"""
function _check_inside(grid,i,j,pt)
    tf = false

    @show grid[i,j], grid[i+1,j], grid[i+1,j+1], grid[i,j+1]

    ab = _checkside(grid[i,j],grid[i+1,j],pt)
    bd = _checkside(grid[i+1,j],grid[i+1,j+1],pt)
    da = _checkside(grid[i+1,j+1],grid[i,j+1],pt)
    ac = _checkside(grid[i,j+1],grid[i,j],pt)
    if sign(ab) == sign(bd) == sign(da) == sign(ac) == -1
        tf = true
    end
    if (sign(ab) == 0) || (sign(bd)==0) || (sign(da)==0) || (sign(ac)==0) #if the point is on one of the lines this will do
        tf = true
    end
    return tf
end
"""
    _checkside(a,b,c)
Compute the determinant of [a-b | a-c]
"""
function _checkside(a,b,c)
    return (a[1] - b[1]) * (c[2] - a[2]) - (a[2] - b[2]) * (c[1] - a[1])
end



