


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
function nearestpoint(grid::Grid2D,pt::Tuple{TT,TT}) where TT

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

    return linind
end

