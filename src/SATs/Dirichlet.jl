

"""
    SAT_Dirichlet
Storage of all objects needed for a Dirichlet SAT ``\\left. u\\right|_{x_i} = g(t)`` where ``i\\in\\{0,1\\}``.
"""
struct SAT_Dirichlet{
        TN<:NodeType,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function, F2<:Function, F3<:Function} <: SimultanousApproximationTerm

    type        :: BoundaryConditionType
    side        :: TN
    axis        :: Int
    order       :: Int
    ED₁ᵀ        :: VT
    RHS         :: F1
    Δx          :: TT
    α           :: TT
    τ           :: F2
    loopaxis    :: F3

    ### CONSTRUCTOR ###
    function SAT_Dirichlet(RHS::F1,Δx::TT,side::TN,axis::Int,order::Int) where {TT,TN,F1}

        check_boundary(side)

        ED = BoundaryDerivativeTranspose(side,order,Δx)
        α,τ = SATpenalties(Dirichlet,Δx,order)

        loopaxis = SelectLoopDirection(axis)

        # fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H D₁ᵀ) H⁻¹ E (u-f)"

        new{TN,TT,Vector{TT},F1,typeof(τ),typeof(loopaxis)}(
            Dirichlet,side,axis,order,ED,RHS,Δx,α,τ,loopaxis)
    end
end







"""
    generate_Dirichlet

Generates mutating functions required for Dirichlet boundary conditions.

If `solver == :cgie` then two methods are generated, one for the boundary data and another for updating the solution.
#TODO: If `solver ∈ [:euler]` then only one method is generated
"""
function generate_Dirichlet(SATD::SAT_Dirichlet,solver)
    # Choose the axis to loop over
    loopdirection = SelectLoopDirection(SATD.axis)

    let α = SATD.α, 
        τ = SATD.τ,
        BD = SATD.ED₁ᵀ,
        side = SATD.side,
        order = SATD.order

        if solver == :cgie
            # Defines 2 methods
            CGTerm!(cache::Array,u::Array,c::Array,::SATMode{:SolutionMode}) = 
                SAT_Dirichlet_implicit!(cache,side,u,c,α,τ,BD,order,loopdirection)
            CGTerm!(cache::Array,data::Array,c::Array,::SATMode{:DataMode}) = 
                    SAT_Dirichlet_implicit_data!(cache,side,data,c,α,τ,BD,order,loopdirection)

                return CGTerm!
        elseif solver ∈ [:euler,:RK4]
            Term!(cache::Array,u::Array,c::Array,t::Float64) = SAT_Dirichlet_explicit!(SATD.RHS,cache,side,u,c,t,α,τ,BD,order,loopdirection)

            return Term!
        end
    end
end


#=== Explicit methods ===#
"""
    SAT_Dirichlet_explicit!
Dirichlet boundary SAT for explicit solvers.
"""
function SAT_Dirichlet_explicit! end
function SAT_Dirichlet_explicit!(RHS::Function,SAT::AbstractArray,::NodeType{:Left},u::AbstractArray,c::AbstractArray,t::Float64,
        α::Float64,τ::Function,BD::AbstractVector,
        order::Int,loopaxis::Function)
        
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[1:order] .+= α*C[1]*(BD*U[1] .- BD*RHS(t))
        S[1] += τ(C[1])*(U[1] - RHS(t))
    end
end
function SAT_Dirichlet_explicit!(RHS::Function,SAT::AbstractArray,::NodeType{:Right},u::AbstractArray,c::AbstractArray,t::Float64,
        α::Float64,τ::Function,BD::AbstractVector,
        order::Int,loopaxis::Function)

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        S[end-order+1:end] .+= α*C[end]*(BD*U[end] .- BD*RHS(t))
        S[end] += τ(C[end])*(U[end] - RHS(t))
    end
end

#=== Implicit methods ===#
"""
    SAT_Dirichlet_implicit!
Solution term for the Dirichlet boundary conditions for SATs for implicit methods. See [`SAT_Dirichlet_implicit_data!`](@ref) for the data term.
"""
function SAT_Dirichlet_implicit! end
function SAT_Dirichlet_implicit!(SAT::AT,::NodeType{:Left},u::AT,c::AT,
        α::T,τ::Function,BD::AbstractVector,
        order::Int,loopaxis::Function) where {T,AT}

    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        for i = 1:order
            S[i] += -α*C[1]*BD[i]*U[1] #D₁ᵀE₀
        end
        S[1] += τ(C[1])*C[1]*U[1]
    end
    SAT
end
function SAT_Dirichlet_implicit!(SAT::AT,::NodeType{:Right},u::AT,c::AT,
        α::T,τ::Function,BD::AbstractVector,
        order::Int,loopaxis::Function) where {T,AT}
    for (S,C,U) in zip(loopaxis(SAT),loopaxis(c),loopaxis(u))
        for i = 1:order
            S[end-order+i] += α*C[end]*BD[i]*U[end] #D₁ᵀEₙ
        end
        # S[end-order+1:end] .+= α*C[end]*BD*U[end]
        S[end] += τ(C[end])*C[end]*U[end]
    end
    SAT
end
"""
    SAT_Dirichlet_implicit_data!
Data term for the Dirichlet boundary conditions for SATs for implicit methods. See [`SAT_Dirichlet_implicit!`](@ref) for the solution term.
"""
function SAT_Dirichlet_implicit_data! end
function SAT_Dirichlet_implicit_data!(SAT::AT,::NodeType{:Left},DATA::AT,c::AT,
        α::T,τ::Function,BD::AbstractVector,order::Int,loopaxis::Function) where {T,AT}

    for (S,U,C) in zip(loopaxis(SAT),loopaxis(DATA),loopaxis(c))
        for i = 1:order
            # S[i] += Δt * α*C[1]*BD[i]*RHS(t)
            S[i] += α*C[1]*BD[i]*U[1] #D₁ᵀE₀
        end
        # S[1] -= Δt* τ*C[1]*RHS(t)#U[1]
        S[1] -= τ(C[1])*C[1]*U[1]
    end
    # println("BD ",BD)
    # println("α ",α)
    # println("C ",τ(c[1]))
    SAT
end
function SAT_Dirichlet_implicit_data!(SAT::AT,::NodeType{:Right},DATA::AT,c::AT,
        α::T,τ::Function,BD::AbstractVector,order::Int,loopaxis::Function) where {T,AT}
    for (S,U,C) in zip(loopaxis(SAT),loopaxis(DATA),loopaxis(c))
        for i = 1:order
            # S[end-order+i] -= Δt* α*C[end]*BD[i]*RHS(t) #D₁ᵀEₙ
            S[end-order+i] -= α*C[end]*BD[i]*U[end] #D₁ᵀEₙ
        end
        # S[end] -= Δt* τ*C[end]*RHS(t)
        S[end] -= τ(C[end])*C[end]*U[end]
    end
    SAT
end














