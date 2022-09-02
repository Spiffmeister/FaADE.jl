





# struct SAT end
# function ()(uₓₓ,)




struct SAT_DirichletStruct{T} <: SimultanousApproximationTerm
    BDₓᵀ    :: Vector{T}
    RHS     :: Function
    type    :: BoundaryCondition
    side    :: NodeType
    order   :: Int
    solver  :: Symbol
    Δx      :: T
    α       :: T
    τ       :: T
end


function construct_Dirichlet()

    if solver == :cgie



    else
    end

end





function (Op::SAT_DirichletStruct)(::NodeType{:Left})
    cache[1:Op.order] .+= Op.α * c[1] * Op.BDₓᵀ * u[1]
    cache[1]        += Op.τ * u[1]
end
function (Op::SAT_DirichletStruct)(::NodeType{:Right})
    cache[end-Op.order+1:end] .+= Op.α * c[end] * Op.BDₓᵀ * u[end]
    cache[end]        += Op.τ * u[end]
end
function (Op::SAT_DirichletStruct)()
    cache[1:order] .-= Op.α * c[1] * Op.BDₓᵀ * u[1]
    cache[1]        -= Op.τ * u[1]
end
function (Op::SAT_DirichletStruct)()
    cache[end-Op.order+1:end] .-= Op.α * c[end] * Op.BDₓᵀ * u[end]
    cache[end]        .-= Op.τ * u[end]
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






function select_SAT_direction(axis::Int)
    if axis == 1
        return eachcol
    elseif axis == 2
        return eachrow
    else
        error("axis must be 1 or 2")
    end
end

