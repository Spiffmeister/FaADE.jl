





# struct SAT end
# function ()(uₓₓ,)




struct SAT_DirichletStruct{T} <: SimultanousApproximationTerm
    BDₓᵀ    :: Vector{T}
    RHS     :: Function
    type    :: BoundaryConditionType
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









abstract type Derivative{R,DIM} end

struct SecondDerivative2D <: Derivative{2,1}
    nx      ::Int
    Δx      ::Float64
    order   ::Int
end

function (D::Derivative{2,1})(uₓₓ,u,c)

end




abstract type Poly{R,T} end

struct Polynomial1D{R} <: Poly{R,1}
    coeffs::Vector{R}
end
struct Polynomial2D{R} <: Poly{R,2}
    coeffs::Vector{R}
end

function (p::Poly{R,1})(x) where R
    v = p.coeffs[end]
    for i = (length(p.coeffs)-1):-1:1
        v = v*x + p.coeffs[i]
    end
    return v
end


function (p::Poly{R,2})(x) where R
    v = p.coeffs[end]
    for i = (length(p.coeffs)-1):-1:1
        v = v*x + p.coeffs[i] + 1
    end
    return v
end

(p::Poly{Int,1})() = p(5)
(p::Poly{Int,2})() = p(1)



p = Polynomial1D{Int}([1,10,100])

p2 = Polynomial2D{Int}([1,10,100])

