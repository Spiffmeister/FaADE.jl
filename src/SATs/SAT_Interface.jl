
"""
SAT(type::BoundaryCondition,::NodeType{:Left},u::AbstractVector{Float64},Δx::Float64,g;
    order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],separate_forcing::Bool=false)
Returns the left SAT (or SAT+Forcing) term as a vector or SAT and Forcing term as separate vectors
Inputs: 
- type: Dirichlet, Neumann, Robin (see ? BoundaryCondition)
- ::NodeType: Left or Right (see ? NodeType)
- u: Solution vector
- Δx: step size
- g: boundary value
Optional inputs:
- order: 2 (default), 4, 6
- coefficient: c in ``\\frac{\\partial}{\\partial x}\\left(c\\frac{\\partial u}{\\partial x}\\right)``
- `αβ`: `[α,β]` Robin boundary coefficient ``\\alpha \\partial_x u + \\beta u = f``
- separate_forcing: true, false (default) -- controls the return type.
false; returns a single vector of the SAT
true; returns SAT and F(orcing) vectors individually
For periodic conditions call Periodic(u,Δx,c;order), for a split domain call Split_domain(u⁻,u⁺,Δx⁻,Δx⁺,c⁻,c⁺;order=2,order⁻=2,order⁺=2)
"""
function SAT end
function SAT(type::BoundaryCondition,node::NodeType,u::AbstractVector{Float64},Δ::Float64,RHS;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    SAT = zeros(Float64,order)
    if type == Dirichlet
        SAT = SAT_Dirichlet(node,u,Δ,RHS,c=c,order=order)
    elseif type == Neumann
        SAT = SAT_Neumann(node,u,Δ,RHS,c=c,order=order)
    elseif type == Robin
        SAT = SAT_Robin(node,u,Δ,RHS,a=αβ[1],b=αβ[2],order=order)
    end
    return SAT
end
function SAT(type::BoundaryCondition,node::NodeType,u::AbstractVector{Float64},Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    SAT = zeros(Float64,order)
    SAT!(SAT,type,node,u,Δ,order=order,c=c,αβ=αβ,forcing=forcing)
    return SAT
end

"""
    SAT!

Iterator for [`SAT`](@ref)
"""
function SAT!(SAT::AbstractVector,type::BoundaryCondition,node::NodeType,u::AbstractVector{Float64},Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractVector{Float64}}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    if type == Dirichlet
        SAT_Dirichlet!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
    elseif type == Neumann
        SAT_Neumann!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
    elseif type == Robin
        # SAT_Robin!(SAT,node,u,Δ,a=αβ[1],b=αβ[2],order=order,forcing=forcing)
    end
    # return SAT
end

