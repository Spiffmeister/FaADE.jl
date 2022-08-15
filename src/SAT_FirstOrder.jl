#=====================================#
#======== BOUNDARY CONDITIONS ========#
# 
# 
#=====================================#


function SAT_FirstOrder end


function SAT_FirstOrder(type::BoundaryCondition,::NodeType{:Left},u::AbstractVector,Δx::Float64,g;
        order=2,seperate_forcing::Bool=false)
    if type == Dirichlet
    elseif type == Neumann
    elseif type == Robin
    end

    if !seperate_forcing
        SAT += F
        return SAT
    else
        return SAT, F
    end

end

