"""
    SAT_Dirichlet

1. SAT_Dirichlet(::NodeType{:Left},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)
2. SAT_Dirichlet(::NodeType{:Left,:Right},u::Vector{Float64},Δx::Float64,RHS;
        c=1.0,order::Int64=2,separate_forcing::Bool=false)
 

Simulatenous approximation term for Dirichlet boundary conditions 
    ``u(xᵢ)=g``
where ``i\\in \\{0,N\\}``.

`NodeType` is either `Left` or `Right`.

See [`NodeType`](@ref)
"""
function SAT_Dirichlet end
# Dirichlet for implicit integrators
function SAT_Dirichlet(node,u::AbstractArray,Δ::Float64; c=1.0,order::Int64=2,forcing=false)
    SAT = zeros(Float64,order)
    SAT = SAT_Dirichlet!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
end
function SAT_Dirichlet(node,u::AbstractVector,Δ::Float64,RHS; c=1.0,order::Int64=2,forcing=false)
    SAT = zeros(Float64,order)
end
# Dirichlet for explicit integrators





"""
    SAT_Dirichlet!


"""
function SAT_Dirichlet!(SAT::AbstractVector,::NodeType{:Left},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)

    SAT[1:order] += BDₓᵀ(u,Left,Δ,order) - BDₓᵀ(RHS,Left,Δ,order)

    SAT[1:order] .*= α * c[1]
    SAT[1]  += τ*(u[1] - RHS[1])
end
function SAT_Dirichlet(SAT::AbstractVector,::NodeType{:Right},u::AbstractVector,Δ::Float64,RHS;
        c=1.0,order::Int64=2,forcing::Bool=false)
    α,τ = SATpenalties(Dirichlet,Δ,order)

    SAT[end-order+1:end] .+= BDₓᵀ(u,Right,Δ,order) + BDₓᵀ(-RHS,Right,Δ,order)

    SAT[end-order+1:end] .*= α * c[end]
    SAT[end]  += τ*(u[end] - RHS[end])
end


#=== Implicit methods ===#
function SAT_Dirichlet!(uₓₓ::AbstractVector,::NodeType{:Left},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)

    α,τ = SATpenalties(Dirichlet,Δ,order)

    if !forcing
        uₓₓ[1:order] .+= α*c[1]*BDₓᵀ(u,Left,Δ,order)
        uₓₓ[1] += τ*u[1]
    else
        uₓₓ[1:order] .-= α*c[1]*BDₓᵀ(u,Left,Δ,order)
        uₓₓ[1] -= τ*u[1]
    end
end
function SAT_Dirichlet!(uₓₓ::AbstractVector,::NodeType{:Right},u::AbstractVector,Δ::Float64;
        c=1.0,order::Int64=2,forcing=false)

    α,τ = SATpenalties(Dirichlet,Δ,order)

    if !forcing
        uₓₓ[end-order+1:end] .+= α*c[end]*BDₓᵀ(u,Right,Δ,order)
        uₓₓ[end] += τ*u[end]
    else
        uₓₓ[end-order+1:end] .-= α*c[end]*BDₓᵀ(u,Right,Δ,order)
        uₓₓ[end] -= τ*u[end]
    end
end