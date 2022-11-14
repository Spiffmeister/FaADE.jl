


"""
    SAT_Neumann

Simulatenous approximation term for Neumann boundary conditions
    ``\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} = g(t)``
where ``i\\in\\{0,N\\}``

`NodeType` is either `Left` or `Right`
"""
struct SAT_Neumann{T} <: SimultanousApproximationTerm
    type    :: BoundaryConditionType
    side    :: NodeType
    axis    :: Int
    order   :: Int
    EDₓ     :: Vector{T}
    RHS     :: Function
    Δx      :: T
    τ       :: T
    ### CONSTRUCTOR ###
    function SAT_Neumann(RHS::Function,Δx::T,side::NodeType,axis::Int,order::Int) where T
        check_boundary(side)

        ED = BoundaryDerivative(side,Δx,order)
        τ = SATpenalties(Neumann,Δx,order)

        new{T}(Neumann,side,axis,order,ED,RHS,Δx,τ)
    end
end


"""
    generate_Neumann
Generates SAT functions for Neumann boundary conditions.
"""
function generate_Neumann end
function generate_Neumann(SATN::SAT_Neumann,solver)
    loopdirection = SelectLoopDirection(SATN.axis)

    let τ = SATN.τ,
        BD = SATN.EDₓ,
        side = SATN.side,
        order = SATN.order

        if solver == :cgie
            CGTerm(cache::AbstractArray,u::AbstractArray,c::AbstractArray,::SATMode{:SolutionMode}) =
                SAT_Neumann_implicit!(cache,side,u,c,τ,BD,order,loopdirection)
            CGTerm(cache::AbstractArray,RHS,c::AbstractArray,::SATMode{:DataMode}) =
                SAT_Neumann_implicit_data!(cache,side,RHS,c,τ,BD,order,loopdirection)

            return CGTerm
        elseif solver ∈ [:euler]
            Term(cache,u,c) = SAT_Neumann_explicit!(cache,side,u,c,SATN.RHS,τ,BD,order,loopdirection)

            return Term
        end
    end
end






#=== Explicit methods ===#
"""
    SAT_Neumann_explicit!
Neumann boundary SAT for explicit solvers.
"""
function SAT_Neumann_explicit! end
function SAT_Neumann_explicit!(SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},RHS,
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[1] = τ*(C[1] * dot(BD,U[1:order]) - RHS[1])
    end
end
function SAT_Neumann_explicit!(SAT::AbstractArray{T},::NodeType{:Right},u::AbstractArray{T},c::AbstractArray{T},RHS,
        τ::T,BD::AbstractArray{T},
        order::Int,loopdirection::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[end] -= τ*(C[end] * dot(BD,U[end-order+1:end]) - RHS[end])
    end
end

#=== Implicit methods ===#
"""
    SAT_Neumann_implicit!

Solution term for the Neumann boundary conditions for SATs for implicit methods. See [`SAT_Neumann_implicit_data!`](@ref) for the data term.
"""
function SAT_Neumann_implicit! end
function SAT_Neumann_implicit!(SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},
    τ::T,BD::AbstractArray{T},
    order::Int,loopaxis::Function) where T
    
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1] += τ*(C[1] * dot(BD,U[1:order]))
    end
end
function SAT_Neumann_implicit!(SAT::AbstractArray{T},::NodeType{:Right},u::AbstractArray{T},c::AbstractArray{T},
    τ::T,BD::AbstractArray{T},
    order::Int,loopaxis::Function) where T
    
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end] -= τ*(C[end] * dot(BD,U[end-order+1:end]))
    end
end
"""
    SAT_Neumann_implicit_data!
Data term for the Neumann boundary conditions for SATs for implicit methods. See [`SAT_Neumann_implicit!`](@ref) for the solution term.
"""
function SAT_Neumann_implicit_data! end
function SAT_Neumann_implicit_data!(SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray,c::AbstractArray{T},
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T
    
    for (S,U) in zip(loopaxis(SAT),loopaxis(u))
        S[1] -= τ*U[1]
    end
end
function SAT_Neumann_implicit_data!(SAT::AbstractArray{T},::NodeType{:Right},u::AbstractArray,c::AbstractArray{T},
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T
    
    for (S,U) in zip(loopaxis(SAT),loopaxis(u))
        S[end] -= -τ*U[end]
    end
end
