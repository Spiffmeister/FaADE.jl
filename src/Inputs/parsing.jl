




function parse_boundaries(B::Dict,G::GridMultiBlock{TT,1}) where TT

    # for J in eachjoint(G)
    for J in 1:length(G.Joint)
        if length(G.Joint[J]) == 1
            haskey(B,J) ? nothing : error("Boundary conditions not found for grid 1.")
        elseif length(G.Joint[J]) == 2
            !haskey(B,J) ? nothing : error("Block $(J) is submerged but boundary conditions have been applied.")
        end
    end

end



function parse_boundaries(B::Dict,G::GridMultiBlock{TT,2}) where TT

    for J in 1:length(G.Joint)
        if haskey(B,J)
            if length(G.Joint[J]) + length(B[J]) - 4 != 0
                error("There are insufficient boundary conditions for grid $(J).")
            end
            for BC in B[J]
                if typeof(BC) <: SAT_Periodic
                    error("Periodic boundary conditions are not allowed for multiblock problems. You should join the boudaries instead.")
                end
            end
        end
    end
    # error("abort")

end





