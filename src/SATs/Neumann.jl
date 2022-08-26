
"""
SAT_Neumann

1. SAT_Neumann(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
    c=1.0,order::Int64=2,forcing::Bool=false)
2. SAT_Neumann(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
    c=1.0,order::Int64=2,forcing::Bool=false)

Simulatenous approximation term for Neumann boundary conditions
``\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} = g(t)``
where ``i\\in\\{0,N\\}``
`NodeType` is either `Left` or `Right`
"""
function SAT_Neumann end
# function SAT_Neumann(::NodeType{:Left},u::Vector{Float64},Δx::Float64;
#         c=1.0,order::Int64=2,forcing::Bool=false)
#     # Penalties
#     τ = SATpenalties(Neumann,Δx,order)
#     SAT = zeros(Float64,order)

#     if !forcing
#         SAT[1] = τ*c[1]*BDₓ(u,Left,Δx,order)
#         return SAT
#     else
#         SAT[1] += -τ*u[1]
#         return SAT
#     end
# end
# function SAT_Neumann(::NodeType{:Right},u::Vector{Float64},Δx::Float64;
#         c=1.0,order::Int=2,forcing::Bool=false)
#     # Penalties
#     τ = SATpenalties(Neumann,Δx,order)
#     SAT = zeros(Float64,order)

#     if !forcing
#         SAT[end] = -τ*c[end]*BDₓ(u,Right,Δx,order)
#         return SAT
#     else
#         SAT[end] -= -τ*u[end] 
#         return SAT
#     end
# end
function SAT_Neumann(::NodeType{:Left},u::AbstractVector,Δx::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)

    SAT[1] = τ*(c[1]*BDₓ(u,Left,Δx,order) - RHS[1])
    return SAT
end
function SAT_Neumann(::NodeType{:Right},u::AbstractVector,Δx::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    τ = SATpenalties(Neumann,Δx,order)
    SAT = zeros(Float64,order)

    SAT[end] = -τ*(c[end]*BDₓ(u,Right,Δx,order) - RHS[end])
    return SAT
end

function SAT_Neumann(node,u::AbstractVector,Δ::Float64;
        c=1.0,order::Int=2,forcing::Bool=false)
    SAT = zeros(Float64,order)
    SAT = SAT_Neumann!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
    return SAT
end

"""
SAT_Neumann!

Iterator for [`SAT_Neumann`](@ref)
"""
function SAT_Neumann! end
function SAT_Neumann!(uₓₓ::AbstractVector,::NodeType{:Left},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int=2,forcing=false)

    τ = SATpenalties(Neumann,Δ,order)

    if !forcing
        uₓₓ[1] += τ*c[1]*BDₓ(u,Left,Δ,order)
        return uₓₓ
    else
        uₓₓ[1] -= τ*u[1]
        return uₓₓ
    end
end
function SAT_Neumann!(uₓₓ::AbstractVector,::NodeType{:Right},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int=2,forcing=false)

    τ = SATpenalties(Neumann,Δ,order)

    if !forcing
        uₓₓ[end] -= τ*c[end]*BDₓ(u,Right,Δ,order)
        return uₓₓ
    else
        uₓₓ[end] += τ*u[end]
        return uₓₓ
    end
end
