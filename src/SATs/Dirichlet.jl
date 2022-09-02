

"""
    Boundary_Dirichlet
"""
struct Boundary_Dirichlet <: SimultanousApproximationTerm
    BDₓᵀ    :: Vector{Real}
    RHS     :: Function
    type    :: BoundaryCondition
    side    :: NodeType
    axis    :: Int
    order   :: Int
    Δx      :: Real
    penalties :: NamedTuple
    expression  :: String

    ### CONSTRUCTOR ###
    function Boundary_Dirichlet(RHS::Function,Δx::Real,side::NodeType,axis::Int,order::Int)

        order ∈ [2,4,6] ? nothing : error("Order must be 2,4 or 6 in position 5")
        side ∈ [Left,Right] ? nothing : error("Must be Left or Right in position 3")
        

        BD = BoundaryDerivativeTranspose(order,Δx)

        α,τ = SATpenalties(Dirichlet,Δx,order)
        penalties = (α=α,τ=τ)

        fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H Dₓᵀ) H⁻¹ E (u-f)"

        new(BD,RHS,Dirichlet,side,axis,order,Δx,penalties,fullsat)
    end
end

# Boundary_Dirichlet(RHS,Δx,side,axis,order)



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




# function (SAT::Boundary_Dirichlet)()




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



function generate_Dirichlet(::BoundaryCondition{:Dirichlet},side,RHS,order,Δx,axis,solver)

    loopdirection = select_SAT_direction(dim)
    α,τ = SATpenalties(Dirichlet,Δx,order)

    if solver == :cgie
        CGTerm(cache,u,c) = SAT_Dirichlet_internal!(cache,side,u,c,Δx,α,τ,BD,order,loopdirection)
        F(cache,c) = SAT_Dirichlet_internal_forcing!(cache,side,RHS,c,Δx,α,τ,BD,order,loopdirection)

        return (CGTerm, F)
    else
        Term(cache,u,c) = SAT_Dirichlet_internal!(cache,side,u,c,Δx,RHS,α,τ,BD,order,loopdirection)

        return Term
    end
end


#=== Explicit methods ===#
function SAT_Dirichlet_internal!(SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},RHS::T,Δ::T,α::T,τ::T,BD::AbstractVector{T},order::Int,loopaxis::Function) where T
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*(BD*U[1] .- BD*RHS)
        S[1] += τ*(U[1] - RHS)
    end
end
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,RHS::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*(BD*U[end] .- BD*RHS)
        S[end] += τ*(U[end] - RHS)
    end
end

#=== Implicit methods ===#
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*BD*U[1]
        S[1] += τ*U[1]
    end
end
function SAT_Dirichlet_internal!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*BD*U[end]
        S[end] += τ*U[end]
    end
end
function SAT_Dirichlet_internal_forcing!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[1:order] .+= α*C[1]*BD*u[1]
        S[1] += τ*u[1]
    end
end
function SAT_Dirichlet_internal_forcing!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[end-order+1:end] .+= α*C[end]*BD*u[end]
        S[end] += τ*u[end]
    end
end





function SAT_Dirichlet!(uₓₓ::AbstractArray,::NodeType{:Right},u::AbstractArray,Δ::Float64,α,τ,BD;
        c=1.0,order::Int64=2,forcing=false)

    # α,τ = SATpenalties(Dirichlet,Δ,order)

    if !forcing
        uₓₓ[end-order+1:end] .+= α*c[end]*[-u[1], u[1]]/Δ
        uₓₓ[end] += τ*u[end]
    else
        uₓₓ[end-order+1:end] .-= α*c[end]*[-u[1], u[1]]/Δ
        uₓₓ[end] -= τ*u[end]
    end
end














