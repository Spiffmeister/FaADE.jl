





# struct SAT end
# function ()(uₓₓ,)




struct SAT_DirichletStruct
    BDₓᵀ    :: Vector
    type    :: BoundaryCondition
    side    :: NodeType
    order   :: Int
    Δx      :: Real
    α       :: Real
    τ       :: Real
end





function (::SAT_DirichletStruct)(::NodeType{:Left})
    cache[1:Op.order] .+= Op.α * c[1] * Op.BDₓᵀ * u[1]
    cache[1]        += Op.τ * u[1]
end
function (::SAT_DirichletStruct)(::NodeType{:Right})
    cache[end-Op.order+1:end] .+= Op.α * c[end] * Op.BDₓᵀ * u[end]
    cache[end]        += Op.τ * u[end]
end
function (::SAT_DirichletStruct)()
    cache[1:order] .-= Op.α * c[1] * Op.BDₓᵀ * u[1]
    cache[1]        -= Op.τ * u[1]
end
function (::SAT_DirichletStruct)()
    cache[end-Op.order+1:end] .-= Op.α * c[end] * Op.BDₓᵀ * u[end]
    cache[end]        -= Op.τ * u[end]
end


function boundary(order,Δ)
    
    BDxT = BDₓᵀ(order,Δ)

    α,τ = SATpenalties()
    if forcing 
        α = -α; τ = -τ
    else
    end
    
    SAT_Dirichlet_internal!(SAT,node,u,c,Δ,α,τ,BDxT,order)
end









