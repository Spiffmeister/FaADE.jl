




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
function nearestpoint(grid::Grid2D,pt::Tuple{TT,TT},indextype=:linear,skip=(nothing,)) where TT

    tmpdist = TT(0)
    dist = TT(1e10)
    linind = 0
    for I in eachindex(grid)
        # Check the distance to each point in the grid and record the smallest distance
        tmpdist = sqrt( (grid.gridx[I] - pt[1])^2 + (grid.gridy[I] - pt[2])^2 )
        if (tmpdist < dist) && (I ∉ skip)
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
    nearestpoint(grid::GridMultiBlock,pt::Tuple{TT,TT},indextype=:linear)
Takes a 2D grid and finds the nearest point.

Returns the value of the nearest point and the linear index of that point.

`indextype=:linear` returns the linear index of the point in the grid. `indextype=:cartesian` returns the i,j index of the point in the grid.
"""
function nearestpoint(grid::GridMultiBlock,pt::Tuple{TT,TT},indextype=:linear) where TT

    tmpdist = TT(0)
    dist = TT(1e10)
    linind = 0
    sgi = 0
    for I in eachgrid(grid)
        tmppt, tmpind = nearestpoint(grid.Grids[I],pt)
        # Check the distance to each point in the grid and record the smallest distance
        tmpdist = sqrt( (tmppt[1] - pt[1])^2 + (tmppt[2] - pt[2])^2 )
        if (tmpdist < dist) && !(tmpdist ≈ dist) #check approx due to rounding errors
            dist = tmpdist
            linind = tmpind
            sgi = I
        end
    end
    if indextype == :linear
        # Return the LinearIndex in the array
        return (grid.Grids[sgi].gridx[linind],grid.Grids[sgi].gridy[linind]), linind, sgi
    elseif indextype == :cartesian
        # Return the point, CartesianIndex and subgrid index
        cartind = CartesianIndices(grid.Grids[sgi].gridx)[linind]
        return (grid.Grids[sgi].gridx[linind],grid.Grids[sgi].gridy[linind]), (cartind[1],cartind[2]), sgi
    else
        error("Index type is either :linear or :cartesian")
    end
end





"""
    findgrid(grid::GridMultiBlock{TT,2,CartesianMetric},pt::Tuple{TT,TT};mode=:inside)
Check which cartesian grid a point is in.
"""
function findgrid(grid::GridMultiBlock{TT,2,CartesianMetric},pt::Tuple{TT,TT};mode=:inside) where {TT}
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
`findgrid` for multiblock problems with Curvilinear coordinates.
"""
function findgrid(grid::GridMultiBlock{TT,DIM,CurvilinearMetric},pt::Tuple{TT,TT};mode=:inside) where {TT,DIM}
    
    _, (i,j), sgi = nearestpoint(grid,pt,:cartesian) #return CartesianIndex
    

    subgrid = grid.Grids[sgi] #Pick out the Grid

    # If the grid point is near the boundary we have to be careful
    if (i == 1) || (i == subgrid.nx) || (j == 1) || (j == subgrid.ny)
        i, j = findcell(subgrid,pt)
        
        if i == j == -1

            joints = grid.Joint[sgi]

            dist = TT(1e10)
            # this will correct most instances
            for joint in joints # Loop over this blocks joints to find the neighbours
                
                # newpt,newind = nearestpoint(grid.Grids[joint.index],pt)
                
                i,j = findcell(grid.Grids[joint.index],pt)
                
                if !(i == j == -1)
                    sgi = joint.index
                    break
                end

            end
        
        end

        # It possible that the grid point lies outside the domain entirely
        # in this case we should move it to the nearest boundary node
        # subgrid = grid.Grids[sgi]
        # try
        #     findcell(subgrid,pt)
        # catch
        # end

        # if the point is outside of the domain, flag it as out of bounds
        if (i == j == -1)
            sgi = -1
        end
        # xside = 0
        # for ii in 1:subgrid.nx-1
        #     xside += sign(_checkside(subgrid[ii,1],subgrid[ii+1,1],pt))
        #     xside += sign(_checkside(subgrid[ii,subgrid.ny],subgrid[ii+1,subgrid.ny],pt))
        # end
        # yside = 0
        # for jj in 1:subgrid.ny-1
        #     yside += _checkside(subgrid[1,jj],subgrid[1,jj+1],pt)
        #     yside += _checkside(subgrid[subgrid.nx,jj],subgrid[subgrid.nx,jj+1],pt)
        # end
        # if abs(xside == subgrid.nx) && abs(yside == subgrid.ny)
        # else
        #     @show "fuck"
        # end
    end
    # if (pt[1] ≈ 0.24451026416330812) & (pt[2] ≈ 0.255555498015997)
        # @show sgi
    # end

    return sgi
end




"""
    findcell(grid::GridMultiBlock,pt::Tuple{TT,TT})
Find which cell bounds the point by first finding the nearest point.
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
    subgrid[i,j], i, j
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
        # ix = jy = -1
    end

    return ix, jy, subgridindex
end

"""
    findcell(grid::Grid2D,pt::Tuple{TT,TT},inds=(0,0)) where {TT}
Find which cell the point belongs too if the nearest point is provided with the `inds` input
"""
function findcell(grid::Grid2D,pt::Tuple{TT,TT},inds=(0,0)) where {TT}

    ix = 0
    jy = 0

    if inds == (0,0)
        _, inds = nearestpoint(grid,pt,:cartesian)
    end

    # _, inds = nearestpoint(grid,pt,:cartesian)
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
        # Last chance is that the point lies on one of the nodes
        if pt == (grid.gridx[i,j],grid.gridy[i,j])
            ix = i
            jy = j
        elseif pt == (grid.gridx[i+1,j],grid.gridy[i+1,j])
            ix = i
            jy = j
        elseif pt == (grid.gridx[i+1,j+1],grid.gridy[i+1,j+1])
            ix = i
            jy = j
        elseif pt == (grid.gridx[i,j+1],grid.gridy[i,j+1])
            ix = i
            jy = j
        else
            # error("Point $(pt) is not in any cell.")
            ix = jy = -1
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

    grid[i,j], grid[i+1,j], grid[i+1,j+1], grid[i,j+1]

    ab = sign(_checkside(grid[i,j],grid[i+1,j],pt))
    bd = sign(_checkside(grid[i+1,j],grid[i+1,j+1],pt))
    da = sign(_checkside(grid[i+1,j+1],grid[i,j+1],pt))
    ac = sign(_checkside(grid[i,j+1],grid[i,j],pt))
    if sign(ab) == sign(bd) == sign(da) == sign(ac)
        tf = true
    end
    if (ab == 0) || (bd==0) || (da==0) || (ac==0) #if the point is on one of the lines this will do
        if abs(ab+bd+da+ac) == 3
            tf = true
        end
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




function _findsegment(grid::Grid2D,pt::Tuple{TT,TT}) where TT
end

