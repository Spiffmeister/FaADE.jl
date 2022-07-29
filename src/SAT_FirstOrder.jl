#=====================================#
#======== BOUNDARY CONDITIONS ========#
# 
# 
#=====================================#


function SAT_FirstOrder end


function SAT_FirstOrder(type::BoundaryCondition,::NodeType{:Left},u::AbstractVector,Δx::Float64,g;order=2,seperate_focing::Bool=false)
    if type == Dirichlet
    elseif type == Neumann
    elseif type == Robin
    end
end

