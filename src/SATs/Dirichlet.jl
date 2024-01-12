

struct Injection_Dirichlet{
        TN<:NodeType,
        F1<:Function}
    side    :: TN
    RHS     :: F1
end



"""
    SAT_Dirichlet
Storage of all objects needed for a Dirichlet SAT ``\\left. u\\right|_{x_i} = g(t)`` where ``i\\in\\{0,1\\}``.

In Cartesian coordinates the SAT reads
``\\tau H^{-1} E H^{-1} E (u - f) + \\alpha H^{-1} (K H D_1^T) H^{-1} E (u - f)``

In Curivlinear coordinates cross derivatives appear in the first term
``\\tau H^{-1} E H^{-1} E (u - f)  ``
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
    Δy          :: TT
    coordinates :: Symbol

    ### CONSTRUCTOR ###
    function SAT_Dirichlet(RHS::F1,Δx::TT,side::TN,axis::Int,order::Int;α=nothing,τ=nothing,Δy=0.0,coord=:Cartesian) where {TT,TN,F1}
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
            side,axis,order,RHS,H⁻¹EH⁻¹E,H⁻¹D₁ᵀE,Δx,α,τ,loopaxis,Δy,coord)
    end
end
SAT_Dirichlet(RHS,Δx,side::NodeType{SIDE,AX},order,Δy=0.0,coord=:Cartesian) where {SIDE,AX} = SAT_Dirichlet(RHS,Δx,side,AX,order,Δy=Δy,coord=coord)





###
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Left},NodeType{:Up}}}
    for i = 1:SD.order
        dest[i] += SD.α*c[1]*SD.H⁻¹D₁ᵀE[i]*(u[1] - RHS[1])
    end
    S[1] += SD.τ(c[1])*SD.H⁻¹EH⁻¹E*(u[1] - RHS[1])
end
function SAT_Dirichlet_explicit!(dest::VT,u::VT,RHS::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:Union{NodeType{:Right},NodeType{:Down}}}
    j = lastindex(dest)
    for i = 1:SD.order
        dest[j-SD.order+i] += SD.α*C[j]*SD.H⁻¹D₁ᵀE[i]*(u[j] - RHS[j]) #D₁ᵀE₀
    end
    S[j] += SD.τ(c[j])*c[j]SD.H⁻¹EH⁻¹E*(u[j] - RHS[j])
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
function SAT_Dirichlet_data!(dest::AT,data::AT,c::AT,SD::SAT_Dirichlet{TN,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(c))
        SAT_Dirichlet_data!(DEST,DATA,C,SD)
    end

    
    if SD.coordinates == :Curvilinear
        n = size(dest,SD.axis)
        m = size(dest,mod1(SD.axis+1,2))
        
        # println("===")
        # @show n, m
        # @show SD.side, SD.axis
        # @show size(dest), size(data)

        if SD.side == Left
            # println("hi")
            DEST = view(dest,   1,1:m)
            SRC = view(data,    1,1:m)
        elseif SD.side == Right
            DEST = view(dest,   n,1:m)
            SRC = view(data,    1,1:m)
        elseif SD.side == Up
            DEST = view(dest,   1:m,1)
            SRC = view(data,    1:m,1)
        else
            DEST = view(dest,   1:m,n)
            SRC = view(data,    1:m,1)
        end
        # @show size(DEST),size(SRC)
        # println("===")

        FirstDerivativeTranspose!(DEST,SRC,m,SD.Δy,SD.order,TT(1))
    end
    dest
end


function SAT_Dirichlet_solution!(dest::VT,data::VT,c::VT,SD::SAT_Dirichlet{TN}) where {VT<:AbstractVector,TN<:NodeType{SIDE}} where SIDE

    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? m = 0 : m = j-SD.order

    dest[j] += SD.τ(c[j])*SD.H⁻¹EH⁻¹E*c[j]*data[j]
    for i = 1:SD.order #nodes SD.nodes
        dest[m+i] += SD.α*c[j]*SD.H⁻¹D₁ᵀE[i]*data[j] #u[Right]
    end
    dest
end
function SAT_Dirichlet_solution!(dest::AT,data::AT,c::AT,SD::SAT_Dirichlet{TN,TT}) where {AT<:AbstractMatrix,TN<:NodeType,TT}
    for (DEST,DATA,C) in zip(SD.loopaxis(dest),SD.loopaxis(data),SD.loopaxis(c))
        SAT_Dirichlet_solution!(DEST,DATA,C,SD)
    end

    if SD.coordinates == :Curvilinear
        n = size(dest,SD.axis)
        m = size(dest,mod1(SD.axis+1,2))
        
        # @show n, m
        # @show SD.side, SD.axis
        # @show size(dest)

        if SD.side == Left
            # println("hi")
            DEST = view(dest,   1:n,1)
            SRC = view(data,    1:n,1)
        elseif SD.side == Right
            DEST = view(dest,   1:n,m)
            SRC = view(data,    1:n,m)
        elseif SD.side == Up
            DEST = view(dest,   1,1:m)
            SRC = view(data,    1,1:m)
        else
            DEST = view(dest,   n,1:m)
            SRC = view(data,    n,1:m)
        end
        # @show size(DEST),size(SRC)

        FirstDerivativeTranspose!(DEST,SRC,n,SD.Δy,SD.order,TT(1))
    end
    dest
end


