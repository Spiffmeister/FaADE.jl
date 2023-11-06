

"""
    SAT_Dirichlet
Storage of all objects needed for a Dirichlet SAT ``\\left. u\\right|_{x_i} = g(t)`` where ``i\\in\\{0,1\\}``.
"""
struct SAT_Dirichlet{
        TN<:NodeType,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function, PT<:Function, LAT<:Function} <: SimultanousApproximationTerm{:Dirichlet}

    side        :: TN
    axis        :: Int
    order       :: Int
    RHS         :: F1
    H⁻¹EH⁻¹E    :: TT
    H⁻¹D₁ᵀE     :: VT
    Δx          :: TT
    α           :: TT
    τ           :: PT
    loopaxis    :: LAT

    ### CONSTRUCTOR ###
    function SAT_Dirichlet(RHS::F1,Δx::TT,side::TN,axis::Int,order::Int;α=nothing,τ=nothing) where {TT,TN,F1}
        # fullsat = "τH⁻¹ E H⁻¹E(u-f) + α H⁻¹ (K H D₁ᵀ) H⁻¹ E (u-f)"

        check_boundary(side)

        loopaxis = SelectLoopDirection(axis)

        if α === nothing
            α = TT(1)
        end
        if τ === nothing
            τ = t -> -TT(1 + 1/max(t,eps(TT)))
        end

        Hinv = _InverseMassMatrix(order,Δx,side)
        D₁ᵀ = _DerivativeTranspose(order,Δx,side)
        E = _BoundaryOperator(TT,side)
        
        # τ * H⁻¹ * E * SD.H⁻¹ * E * (u-f)
        if typeof(side) <: NodeType{:Left}
            H⁻¹EH⁻¹E    = Hinv[1]*E*Hinv[1]*E
        elseif typeof(side) <: NodeType{:Right}
            H⁻¹EH⁻¹E    = Hinv[end]*E*Hinv[end]*E
        end
        
        # α H⁻¹ * (K H D₁ᵀ) * H⁻¹ * E * (u-f)
        H⁻¹D₁ᵀE     = Hinv.*D₁ᵀ.*E

        new{TN,TT,Vector{TT},typeof(RHS),typeof(τ),typeof(loopaxis)}(
            side,axis,order,RHS,H⁻¹EH⁻¹E,H⁻¹D₁ᵀE,Δx,α,τ,loopaxis)
    end
end
SAT_Dirichlet(RHS,Δx,side::NodeType{SIDE,AX},order) where {SIDE,AX} = SAT_Dirichlet(RHS,Δx,side,AX,order)





###
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Left},NodeType{:Up}}}
    for i = 1:SD.order
        dest[i] += SD.α*c[1]*SD.H⁻¹D₁ᵀE[i]*(U[1] - RHS[1])
    end
    S[1] += SD.τ(c[1])*SD.H⁻¹EH⁻¹E*(u[1] - RHS[1])
end
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Right},NodeType{:Down}}}
    j = lastindex(dest)
    for i = 1:SD.order
        dest[j-SD.order+i] += SD.α*C[j]*SD.H⁻¹D₁ᵀE[i]*(U[j] - RHS[j]) #D₁ᵀE₀
    end
    S[j] += SD.τ(c[j])*c[j]SD.H⁻¹EH⁻¹E*(U[j] - RHS[j])
end




# Operation along vector
function SAT_Dirichlet_data!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? m = 0 : m = j-SD.order
    SIDE == :Left ? o = 1 : o = lastindex(data)

    dest[j] -= SD.τ(c[j])*SD.H⁻¹EH⁻¹E*c[j]*data[o]
    for i = 1:SD.order #nodes SD.nodes
        dest[m+i] -= SD.α*c[j]*SD.H⁻¹D₁ᵀE[i]*data[o] #u[Left]
    end
    dest
end

# function SAT_Dirichlet_data!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{:Right}}
#     for i = 1:SD.order #nodes SD.nodes
#         dest[end-SD.order+i] = dest[i] + SD.α*c[end]*SD.ED₁ᵀ[i]*data[end] #u[Left]
#     end
#     dest[end] = dest[1] - SD.τ(c[end])*c[end]*data[end]
# end


function SAT_Dirichlet_solution!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? m = 0 : m = j-SD.order
    # SIDE == :Left ? o = 1 : o = lastindex(data)

    dest[j] += SD.τ(c[j])*SD.H⁻¹EH⁻¹E*c[j]*data[j]
    for i = 1:SD.order #nodes SD.nodes
        dest[m+i] += SD.α*c[j]*SD.H⁻¹D₁ᵀE[i]*data[j] #u[Right]
    end
    dest
end

# function SAT_Dirichlet_solution!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{:Right}}
#     τ = SD.τ(c[end])
#     for i = 1:SD.order #nodes SD.nodes
#         dest[end-SD.order+i] = dest[end-SD.order+i] + SD.α*c[end]*SD.ED₁ᵀ[i]*data[end] #u[Right] dest + SD.α * SD.h[i]/SD.Δx * SD.ED₁ᵀ[i] * data[end]
#     end
#     dest[end] = dest[end] + τ*c[end]*data[end]
# end


