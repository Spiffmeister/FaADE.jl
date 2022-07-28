#=====================================#
#======== BOUNDARY CONDITIONS ========#
# 
# 
#=====================================#


function Boundary(condtype::Symbol,order::Int64;cond::Union{Function,Vector{Float64}})

    if condtype ∉ (:Dirichlet,:Periodic,:Neumann,:Robin)
        error("Must be in Dirichlet, Periodic, Neumann or Robin")
    end

end


function Boundary_dirichlet(;order::Int64=2,derivative::Int64=1)
    # Dirichlet boundary operators for 2, 4, 6 and 8th order SPB operators
    if order == 2

        coeff = ([1.0,1.0])
        index = ([1,2])
        
    elseif order == 4

        coeff = ([-24/17, 53/34, -4/17, -3/34],
            [-1/2, 1/2],
            [4/43, -59/86, 59/86, -4/43],
            [3/98, -59/98, 32/49, -4/49])
        index = ([1,2,3,4],
            [1,3],
            [1,2,4,5],
            [1,3,5,6])
        
    elseif order == 6
        
    end

    BC = BoundaryCond(coeff,index)
    return BC
end



function B_perioidic(;order::Int64=2)

    if order == 2
    elseif order == 4
    end

end