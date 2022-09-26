

"""
    SAT_Dirichlet
"""
struct SAT_Dirichlet{T} <: SimultanousApproximationTerm
    type    :: BoundaryConditionType
    side    :: NodeType
    axis    :: Int
    order   :: Int
    EDₓᵀ    :: Vector{T}
    RHS     :: Function
    Δx      :: T
    α       :: T
    τ       :: T

    ### CONSTRUCTOR ###
    function SAT_Dirichlet(RHS::Function,Δx::T,side::NodeType,axis::Int,order::Int) where T

        check_boundary(side)

        ED = BoundaryDerivativeTranspose(order,Δx)
        α,τ = SATpenalties(Dirichlet,Δx,order)

        # fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H Dₓᵀ) H⁻¹ E (u-f)"

        new{T}(Dirichlet,side,axis,order,ED,RHS,Δx,α,τ)
    end
end







"""
    generate_Dirichlet

Generates mutating functions required for Dirichlet boundary conditions.

If `solver == :cgie` then two methods are generated, one for the boundary data and another for updating the solution.
If `solver ∈ [:euler]` then only one method is generated
"""
function generate_Dirichlet(SATD::SAT_Dirichlet,solver)
    # Choose the axis to loop over
    loopdirection = SelectLoopDirection(SATD.axis)

    let α = SATD.α, 
        τ = SATD.τ,
        BD = SATD.EDₓᵀ,
        side = SATD.side,
        order = SATD.order

        if solver == :cgie
            # Defines 2 methods
            CGTerm(cache::Array,u::Array,c::Array,::SATMode{:SolutionMode}) = 
                SAT_Dirichlet_implicit!(cache,side,u,c,α,τ,BD,order,loopdirection)
            CGTerm(cache::Array,RHS,c::Array,::SATMode{:DataMode}) = 
                    SAT_Dirichlet_implicit_data!(cache,side,RHS,c,α,τ,BD,order,loopdirection)

                return CGTerm
        elseif solver ∈ [:euler]
            Term(cache,u,c) = SAT_Dirichlet_explicit!(cache,side,u,c,SATD.RHS,α,τ,BD,order,loopdirection)

            return Term
        end
    end
end


#=== Explicit methods ===#
function SAT_Dirichlet_explicit! end
function SAT_Dirichlet_explicit!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,RHS,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)
        
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*(BD*U[1] .- BD*RHS)
        S[1] += τ*(U[1] - RHS)
    end
end
function SAT_Dirichlet_explicit!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,RHS,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*(BD*U[end] .- BD*RHS)
        S[end] += τ*(U[end] - RHS)
    end
end

#=== Implicit methods ===#
function SAT_Dirichlet_implicit! end
function SAT_Dirichlet_implicit!(SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*BD*U[1]
        S[1] += τ*U[1]
    end
end
function SAT_Dirichlet_implicit!(SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*BD*U[end]
        S[end] += τ*U[end]
    end
end
function SAT_Dirichlet_implicit_data! end
function SAT_Dirichlet_implicit_data!(SAT::AbstractArray,::NodeType{:Left},u,c::AbstractArray,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[1:order] .-= α*C[1]*BD*u[1]
        S[1] -= τ*u[1]
    end
end
function SAT_Dirichlet_implicit_data!(SAT::AbstractArray,::NodeType{:Right},u,c::AbstractArray,
        α::Float64,τ::Float64,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C) in zip(loopaxis(SAT),loopaxis(c))
        S[end-order+1:end] .-= α*C[end]*BD*u[end]
        S[end] -= τ*u[end]
    end
end
















