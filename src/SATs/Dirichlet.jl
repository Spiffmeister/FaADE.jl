

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

WARNING: DEPRECATED

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
    # SAT_Dirichlet_implicit!(SAT,node,u,Δ,c=c,order=order,forcing=forcing)
end
function SAT_Dirichlet(node,u::AbstractVector,Δ::Float64,RHS; c=1.0,order::Int64=2,forcing=false)
    SAT = zeros(Float64,order)
end
# Dirichlet for explicit integrators






"""
    generate_Dirichlet

Generates mutating functions required for Dirichlet boundary conditions.

If `solver == :cgie` then two methods are generated, one for the boundary data and another for updating the solution.
If `solver ∈ [:euler]` then only one method is generated
"""
function generate_Dirichlet(SATD::Boundary_Dirichlet,solver)
    # Choose the axis to loop over
    loopdirection = select_SAT_direction(SATD.dim)

    α,τ = SATD.penalties
    BD = SATD.BDₓᵀ
    side = SATD.side
    Δx = SATD.Δx

    if solver == :cgie
        # Defines 2 methods
        CGTerm(cache,u,c) = SAT_Dirichlet_implicit!(cache,side,u,c,Δx,α,τ,BD,order,loopdirection)
        CGTerm(cache,c) = SAT_Dirichlet_implicit_forcing!(cache,side,SATD.RHS,c,Δx,α,τ,BD,order,loopdirection)

        return CGTerm
    elseif solver ∈ [:euler]
        Term(cache,u,c) = SAT_Dirichlet_explicit!(cache,side,u,c,Δx,SATD.RHS,α,τ,BD,order,loopdirection)

        return Term
    end
end


#=== Explicit methods ===#
function SAT_Dirichlet_explicit!(SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},RHS::T,Δ::T,α::T,τ::T,BD::AbstractVector{T},order::Int,loopaxis::Function) where T
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*(BD*U[1] .- BD*RHS)
        S[1] += τ*(U[1] - RHS)
    end
end
function SAT_Dirichlet_explicit!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,RHS::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*(BD*U[end] .- BD*RHS)
        S[end] += τ*(U[end] - RHS)
    end
end

#=== Implicit methods ===#
function SAT_Dirichlet_implicit!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*BD*U[1]
        S[1] += τ*U[1]
    end
end
function SAT_Dirichlet_implicit!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*BD*U[end]
        S[end] += τ*U[end]
    end
end
function SAT_Dirichlet_implicit_forcing!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[1:order] .+= α*C[1]*BD*u[1]
        S[1] += τ*u[1]
    end
end
function SAT_Dirichlet_implicit_forcing!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,Δ::Float64,α::Float64,τ::Float64,BD::AbstractVector,order::Int,loopaxis::Function)
    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[end-order+1:end] .+= α*C[end]*BD*u[end]
        S[end] += τ*u[end]
    end
end
















