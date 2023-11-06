


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
    F1<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Neumann}
    
    side    :: TN
    axis    :: Int
    order   :: Int
    RHS     :: F1
    H⁻¹E    :: TT
    D₁      :: VT
    Δx      :: TT
    τ       :: TT
    loopaxis:: LAT
    ### CONSTRUCTOR ###
    function SAT_Neumann(RHS::F1,Δx::TT,side::TN,axis::Int,order::Int) where {TT,TN,F1}

        check_boundary(side)

        LA = SelectLoopDirection(axis)


        Hinv = _InverseMassMatrix(order,Δx,side)
        E = _BoundaryOperator(TT,side)
        D₁ = _BoundaryDerivative(Δx,order,side)


        # τ * H⁻¹ * E
        if typeof(side) <: NodeType{:Left}
            n = -1 # normal vector
            H⁻¹E    = n*Hinv[1]*E
        elseif typeof(side) <: NodeType{:Right}
            n = -1 # normal vector
            H⁻¹E    = n*Hinv[end]*E
        end

        τ = TT(1)

        new{TN,TT,Vector{TT},F1,typeof(LA)}(side,axis,order,RHS,H⁻¹E,D₁,Δx,τ,LA)
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

#== NEW ==#
function SAT_Neumann_implicit_solution!(dest::AT,u::AT,c::AT,SN::SAT_Neumann{TN}) where {AT,TN<:Union{NodeType{:Left},NodeType{:Up}}}
    for (S,U,C) in zip(SN.loopaxis(dest),SN.loopaxis(u),SN.loopaxis(c))
        for i = 1:SN.order
            S[1] += SN.τ*(C[1] * SN.ED₁[i]*U[i])
        end
    end
end
function SAT_Neumann_implicit_solution!(dest::AT,u::AT,c::AT,SN::SAT_Neumann{TN}) where {AT,TN<:Union{NodeType{:Right},NodeType{:Down}}}
    for (S,U,C) in zip(SN.loopaxis(dest),SN.loopaxis(u),SN.loopaxis(c))
        for i = 1:SN.order
            S[end] -= SN.τ*(C[end] * SN.ED₁[i]*U[end-SN.order+i])
        end
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

# #== NEW ==#
# function SAT_Neumann_data!(dest::AT,u::AT,SN::SAT_Neumann{TN}) where {AT,TN<:NodeType{SIDE}} where SIDE
#     SIDE == :Left ? j = 1 : j = lastindex(dest)
#     SIDE == :Left ? b =-1 : b = -1
#     dest[j] -= b*SN.τ*SN.H⁻¹E*u[1]
#     # dest[j] -= SN.τ*SN.H⁻¹E*u[1] #-1 for right
# end
# #== NEW ==#
# function SAT_Neumann_solution!(dest::VT,u::VT,c::VT,SN::SAT_Neumann{TN}) where {VT,TN<:NodeType{SIDE}} where SIDE
#     SIDE == :Left ? j = 1 : j = lastindex(dest)
#     SIDE == :Left ? m = 0 : m = j-SN.order
#     SIDE == :Left ? b = -1 : b = 1
#     for i = 1:SN.order
#         dest[j] += b*SN.τ * SN.H⁻¹E * c[1] * SN.D₁[i]*u[m+i]
#         # dest[j] += SN.τ * SN.H⁻¹E * c[j] * SN.D₁[i]*u[m+i] #-1 for left
#     end
# end



#== NEW ==#
function SAT_Neumann_data!(dest::AT,u::AT,SN::SAT_Neumann{TN}) where {AT,TN<:NodeType{:Left}}
    dest[1] -= SN.τ*SN.H⁻¹E*u[1]
end
function SAT_Neumann_solution!(dest::VT,u::VT,c::VT,SN::SAT_Neumann{TN}) where {VT,TN<:NodeType{:Left}}
    for i = 1:SN.order
        dest[1] += SN.τ * SN.H⁻¹E * c[1] * SN.D₁[i]*u[i]
    end
end
function SAT_Neumann_data!(dest::AT,u::AT,SN::SAT_Neumann{TN}) where {AT,TN<:NodeType{:Right}}
    dest[end] -= SN.τ*SN.H⁻¹E*u[end]
end
function SAT_Neumann_solution!(dest::VT,u::VT,c::VT,SN::SAT_Neumann{TN}) where {VT,TN<:NodeType{:Right}}
    for i = 1:SN.order
        dest[end] += SN.τ * SN.H⁻¹E * c[end] * SN.D₁[i]*u[end-SN.order+i]
    end
end