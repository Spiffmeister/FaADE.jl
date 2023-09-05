


"""
    SAT_Neumann

Simulatenous approximation term for Neumann boundary conditions
    ``\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} = g(t)``
where ``i\\in\\{0,N\\}``

`NodeType` is either `Left` or `Right`
"""
struct SAT_Neumann{
    TN<:NodeType,
    TT<:Real,
    VT<:Vector{TT},
    F1<:Function} <: SimultanousApproximationTerm{:Neumann}
    
    type    :: BoundaryConditionType
    side    :: TN
    axis    :: Int
    order   :: Int
    ED₁     :: VT
    RHS     :: F1
    Δx      :: TT
    τ       :: TT
    ### CONSTRUCTOR ###
    function SAT_Neumann(RHS::F1,Δx::TT,side::TN,axis::Int,order::Int) where {TT,TN,F1}

        check_boundary(side)

        ED = BoundaryDerivative(side,Δx,order)
        τ = SATpenalties(Neumann,Δx,order)

        new{TN,TT,Vector{TT},F1}(Neumann,side,axis,order,ED,RHS,Δx,τ)
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
        BD = SATN.ED₁,
        side = SATN.side,
        order = SATN.order

        if solver == :cgie
            CGTerm(cache::AbstractArray,u::AbstractArray,c::AbstractArray,::SATMode{:SolutionMode}) =
                SAT_Neumann_implicit!(cache,side,u,c,τ,BD,order,loopdirection)
            CGTerm(cache::AbstractArray,RHS,c::AbstractArray,::SATMode{:DataMode}) =
                SAT_Neumann_implicit_data!(cache,side,RHS,c,τ,BD,order,loopdirection)

            return CGTerm
        elseif solver ∈ [:euler,:RK4]
            Term(cache,u,c,t) = SAT_Neumann_explicit!(SATN.RHS,cache,side,u,c,t,τ,BD,order,loopdirection)

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
function SAT_Neumann_explicit!(RHS::Function,SAT::AbstractArray{T},::NodeType{:Left},u::AbstractArray{T},c::AbstractArray{T},t::Float64,
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[1] = τ*(C[1] * dot(BD,U[1:order]) - RHS(t))
    end
end
function SAT_Neumann_explicit!(RHS::Function,SAT::AbstractArray{T},::NodeType{:Right},u::AbstractArray{T},c::AbstractArray{T},t::Float64,
        τ::T,BD::AbstractArray{T},
        order::Int,loopaxis::Function) where T

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(u),loopaxis(c))
        S[end] -= τ*(C[end] * dot(BD,U[end-order+1:end]) - RHS(t))
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


"""
Explicit Neumann SATs ----- NEEDS TESTING
"""
function (SN::SAT_Neumann{NodeType{:Left,DIM},TT})(cache::AT,u::AT,c::AT,t::TT) where {TT,DIM,AT}
    for (S,C,U) in zip(SN.loopaxis(cache),SN.loopaxis(c),SN.loopaxis(u))
        S[1] += SN.τ*(C[1]*dot(SN.ED₁,U[1:SN.order]) - SN.RHS(t))
    end
end
function (SN::SAT_Neumann{NodeType{:Right,DIM},TT})(cache::AT,u::AT,c::AT,t::TT) where {TT,DIM,AT}
    for (S,C,U) in zip(SN.loopaxis(cache),SN.loopaxis(c),SN.loopaxis(u))
        S[end] += -SN.τ*(C[end]*dot(SN.ED₁,U[end-SN.order+1:end]) - SN.RHS(t))
    end
end
"""
Implicit Neumann SATs ----- NEEDS TESTING
"""
function (SN::SAT_Neumann{NodeType{:Left,DIM},TT})(cache::AT,c::AT,u::AT,::SATMode{:SolutionMode}) where {TT,DIM,AT}
    for (S,C,U) in zip(SN.loopaxis(cache),SN.loopaxis(c),SN.loopaxis(u))
        S[1] += SN.τ * C[1] * dot(SN.ED₁,U[1:order])
    end
end
function (SN::SAT_Neumann{NodeType{:Right,DIM},TT})(cache::AT,c::AT,u::AT,::SATMode{:SolutionMode}) where {TT,DIM,AT}
    for (S,C,U) in zip(SN.loopaxis(cache),SN.loopaxis(c),SN.loopaxis(u))
        S[end] += -SN.τ * C[end] * dot(SN.ED₁,U[end-SN.order+1:end])
    end
end
function (SN::SAT_Neumann{NodeType{:Left,DIM},TT})(cache::AT,c::AT,u::AT,::SATMode{:DataMode}) where {TT,DIM,AT}
    for (S,U) in zip(SN.loopaxis(cache),SN.loopaxis(u))
        S[1] += -SN.τ*U[1]
    end
end
function (SN::SAT_Neumann{NodeType{:Right,DIM},TT})(cache::AT,c::AT,u::AT,::SATMode{:DataMode}) where {TT,DIM,AT}
    for (S,U) in zip(SN.loopaxis(cache),SN.loopaxis(u))
        S[end] += SN.τ*U[end]
    end
end
