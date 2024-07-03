


# function Master_grid(Grid::GridMultiBlock;shape=:square,domain=:connected)
# end





# function _form_master_square(Grid::GridMultiBlock{TT},Val{:square},Val{:connected}) where TT

#     # local masterLeft :: Int
#     doubleindex_x = 0
#     doubleindex_y = 0

#     for I in eachgrid(Grid)
#         for Joint in Grid.Joint[I]
#             if Joint.side == Right
#                 doubleindex_x += 1
#             elseif Joint.side == Up
#                 doubleindex_y += 1
#             end
#         end
#     end

#     n = size(Grid)
#     masterx = zeros(TT,n[1]-doubleindex_x,n[2]-doubleindex_y)
#     mastery = zeros(TT,n[1]-doubleindex_x,n[2]-doubleindex_y)

#     for I in eachgrid(Grid)
#         if I == 1
#             masterx[1:Dom.inds[1][1],1:Dom.inds[2][1]] = Grid.Grids[I].gridx
#             mastery[1:Dom.inds[1][1],1:Dom.inds[2][1]] = Grid.Grids[I].gridy
#         else
#             masterx[Dom.inds[1][I]:Dom.inds[1][I+1],Dom.inds[2][I]:Dom.inds[2][I+1]] = Grid.Grids[I].gridx
#             mastery[Dom.inds[1][I]:Dom.inds[1][I+1],Dom.inds[2][I]:Dom.inds[2][I+1]] = Grid.Grids[I].gridy
#         end
#     end

# end

# function _findside(Grid::GridMultiblock)
# end


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
        tmpdist = sqrt( (grid.gridx[I] - pt[1])^2 + (grid.gridy[I] - pt[2])^2 )
        if tmpdist < dist
            dist = tmpdist
            linind = I
        end
    end

    if indextype == :linear
        return (grid.gridx[linind],grid.gridy[linind]), linind
    elseif indextype == :cartesian
        cartind = CartesianIndices(grid.gridx)[linind]
        return (grid.gridx[linind],grid.gridy[linind]), (cartind[1],cartind[2])
    else
        error("Index type is either :linear or :cartesian")
    end
end
"""
    nearestpoint(grid::GridMultiBlock,pt::Tuple{TT,TT})
Takes a multiblock grid and finds the nearest point.

Returns the value of the nearest point, the linear index of that point, and the grid number in GridMultiBlock.
"""
function nearestpoint(grid::GridMultiBlock,pt::Tuple{TT,TT}) where TT
    point = TT(0)
    dist = TT(1e10)
    index = 0
    gridindex = 0
    for I in eachgrid(grid)
        tmppt,tmpindex = nearestpoint(grid.Grids[I],pt)

        tmpdist = sqrt( (tmppt[1] - pt[1])^2 + (tmppt[2] - pt[2])^2 )
        if (tmpdist < dist ) | (point == 0)
            dist = tmpdist
            point = tmppt
            index = tmpindex
            gridindex = I
        end
    end

    return point, index, gridindex
end

function findgrid(grid::GridMultiBlock,pt::Tuple{TT,TT};mode=:inside) where TT
    gridind = 0
    for I in eachgrid(grid)
        minx = grid.Grids[I].gridx[1]
        maxx = grid.Grids[I].gridx[end]
        miny = grid.Grids[I].gridy[1]
        maxy = grid.Grids[I].gridy[end]

        if (minx <= pt[1] <= maxx) & (miny <= pt[2] <= maxy)
            gridind = I
        end
    end
    if gridind != 0
        return gridind
    elseif (gridind == 0) & (mode == :nearest)
        gridind = nearestpoint(grid,pt)[3]
        return nearestpoint(grid,pt)[3]
    elseif (gridind == 0) & (mode == :inside)
        error("Point $(pt) is not in any grid")
    end
end

function findcell(grid::GridMultiBlock,pt::Tuple{TT,TT}) where TT
    gridind = findgrid(grid,pt)
    grid = grid.Grids[gridind]


    for I in eachindex(grid.gridx)
        if (grid.gridx[I] <= pt[1] <= grid.gridx[I+1]) & (grid.gridy[I] <= pt[2] <= grid.gridy[I+1])
            return I
        end
    end



    return pt,index,gridind
end