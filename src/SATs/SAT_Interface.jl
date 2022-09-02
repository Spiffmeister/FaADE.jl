



abstract type SimultanousApproximationTerm end




struct SimultanousApproximationTermContainer
    SATs        :: Vector{SimultanousApproximationTerm}
    axis        :: Vector{Int}
    edge        :: Vector{NodeType}
    type        :: BoundaryCondition

    function SimultanousApproximationTermContainer(boundaries...)
        
        ApproxTerm = []
        axis = []
        edge = []

        for term in boundaries
            push!(ApproxTerm,term)
            push!(axis,term.axis)
            push!(edge,term.side)
        end
        new(ApproxTerm,axis,edge)
    end

end





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
function SAT!(SAT::AbstractArray,type::BoundaryCondition,node::NodeType,u::AbstractArray,Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractArray}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)

    if type == Dirichlet
        SAT_Dirichlet!(SAT,node,u,Δ,α,τ,c=c,order=order,forcing=forcing)
    elseif type == Neumann
        SAT_Neumann!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
    elseif type == Robin
        # SAT_Robin!(SAT,node,u,Δ,a=αβ[1],b=αβ[2],order=order,forcing=forcing)
    end
    # return SAT
end

function SATAdd!(SAT::AbstractArray,::BoundaryCondition{:Dirichlet},node::NodeType,u::AbstractArray,Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractArray}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)

    α,τ = SATpenalties(Dirichlet,Δ,order)
    BD = BDₓᵀ(2,Δ)

    # FD(A,U,K) = SAT_Dirichlet_internal!(A,node,U,K,Δ,α,τ,BD,order)

    if !forcing
        map((A,U,K) -> SAT_Dirichlet!(A,node,U,Δ,α,τ,BD,c=K,order=order,forcing=forcing), eachcol(SAT),eachcol(u),eachcol(c))
        # map((A,U,K) -> SAT_Dirichlet_internal!(A,node,U,K,Δ,α,τ,BD,order), eachcol(SAT),eachcol(u),eachcol(c))
        # map((A,U,K) -> FD(A,U,K), eachcol(SAT),eachcol(u),eachcol(c))
        # SAT_Dirichlet!(SAT,node,u,Δ,α,τ,BD,c=c,order=order,forcing=forcing)
    else
        map((A,K) -> SAT_Dirichlet!(A,node,u,Δ,α,τ,BD,c=K,order=order,forcing=forcing), eachcol(SAT),eachcol(c))
    end
    SAT
end
function SATAdd!(SAT::AbstractArray,::BoundaryCondition{:Neumann},node::NodeType,u::AbstractArray,Δ::Float64;
        order=2::Int,c::Union{Float64,AbstractArray}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    if !forcing
        map((A,U) -> SAT_Neumann!(A,node,U,Δ,c=c,order=order,forcing=forcing), eachcol(SAT),eachcol(u))
    else
        map((A,U,K) -> SAT_Neumann!(A,node,U,Δ,c=K,order=order,forcing=forcing), eachcol(SAT),eachcol(u),eachcol(c))
    end
    SAT
end
# function SATAdd!(SAT::AbstractArray,::BoundaryCondition{:Periodic},node::NodeType,u::AbstractArray,Δ::Float64;
#         order=2::Int,c::Union{Float64,AbstractArray}=1.0,αβ::Vector{Float64}=[1.0,1.0],forcing=false)
    
#         map((A,U) -> SAT_Periodic!(A,node,U,Δ,c=c,order=order,forcing=forcing), eachcol(SAT),eachcol(u))
#     SAT
# end


# function SAT(u::AbstractArray,type::Dirichlet,node::NodeType,Δ::Float64,n::Int;
#         order::Int=2,c::Union{Float64,AbstractMatrix})
#     SAT = zeros((order,n))
#     for i = 1:n
#         SAT[1:order] = SAT_Dirichlet!(SAT[1:order],node,u,Δ,c=c,order=order,forcing=forcing)
#     end
# end



# function dimiter(u,:first)
#     return eachcol
# end
# function dimiter(u,:second)
#     return eachcol
# end







# struct Boundary 
#     nodes :: CartesianIndices
#     SAT :: Function
#     Forcing :: Function

#     function Boundary()

#         if condition == Dirichlet
#             BC(uₓₓ,u,c) = SAT_Dirichlet!()

#         elseif condition == Neumann
#         end

#     end

# end



function BDx(order,u)
    if order == 2
        return [-1.0,1.0]*u[1]
    elseif order == 4
        return [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]*u[1]
    elseif order == 6
        return [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]*u[1]
    end
end






