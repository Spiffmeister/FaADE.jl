
"""
    SAT_Periodic{TN<:NodeType,COORD,TT<:Real,VT<:Vector{TT},F1<:Function,F2<:Function} <: SimultanousApproximationTerm{:Periodic}

Storage of all objects needed for a Periodic SAT ``u(x_0) = u(x_N)`` and ``\\left.\\partial_x u\\right|_{x_0} = \\left.\\partial_x u\\right|_{x_N}``.

In Cartesian the SAT is
``\\tau_0 H^{-1}L_0u + \\tau_1 H^{-1}KD_x^TL_1u \\alpha_0 H^{-1}L_1KD_xu``

In curvilinear coordinates the SAT is
``\\tau_0 H^{-1}L_0u + \\tau_1 H^{-1}K_qD_q^TL_1u + \\tau_2H^{-1}K_{qr}D_r^TL_1 u + \\alpha_1 H^{-1}L_1K_qD_qu + \\alpha_2 H^{-1}L_1K_{qr}D_ru``
"""
struct SAT_Periodic{
        TN<:NodeType,
        COORD,
        TT<:Real,
        VT<:Vector{TT},
        F1<:Function,F2<:Function} <: SimultanousApproximationTerm{:Periodic}
    type        :: BoundaryConditionType
    side        :: TN
    axis        :: Int
    order       :: Int
    D₁ᵀE₀       :: VT
    D₁ᵀEₙ       :: VT
    E₀D₁        :: VT
    EₙD₁        :: VT
    Δx          :: TT
    α₀          :: TT
    τ₁          :: TT
    τ₀          :: F1
    loopaxis    :: F2

    Δy          :: TT
    coordinates :: Symbol

    @doc """
        SAT_Periodic(Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian) where {TT,SIDE,AX}
    
    Constructor for periodic SAT

    Required arguments:
    - `Δx`: Grid spacing
    - `side`: Boundary side, either `Left`, `Right`, `Up`, `Down`
    - `order`: Order of the ``D_x`` operators, either `2` or `4`
    """
    function SAT_Periodic(Δx::TT,side::NodeType{SIDE,AX},order::Int,Δy=TT(0),coord=:Cartesian) where {TT,SIDE,AX}
        D₁ᵀE₀ = BoundaryDerivativeTranspose(Left,order,Δx)
        D₁ᵀEₙ = BoundaryDerivativeTranspose(Right,order,Δx)
        E₀D₁ = BoundaryDerivative(Left,Δx,order)
        EₙD₁ = BoundaryDerivative(Right,Δx,order)
    
        α₀, τ₁, τ₀ = SATpenalties(Periodic,Δx,order)
    
        loopaxis = SelectLoopDirection(AX)
    
        new{typeof(side),coord,TT,Vector{TT},typeof(τ₀),typeof(loopaxis)}(
            Periodic,side,AX,order,D₁ᵀE₀,D₁ᵀEₙ,E₀D₁,EₙD₁,Δx,α₀,τ₁,τ₀,loopaxis,Δy,coord)
    end
end








"""
Periodic SAT.

1D methods for `Left` and `Down` boundaries and `Right` and `Up` boundaries respectively are

    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Left}},::SATMode{:SolutionMode}) where {AT}
    SAT_Periodic!(cache::AT,u::AT,c::AT,SP::SAT_Periodic{NodeType{:Right}},::SATMode{:SolutionMode}) where {AT}

The 2D method loops over the 1D methods by `loopaxis`

    SAT_Periodic!(dest::AT,u::AT,c::AT,SP::SAT_Periodic{TN,:Cartesian,TT}) where {TT,AT<:AbstractMatrix{TT},TN<:NodeType}

The curvilinear method wraps the 2D Cartesian method and calls the [`FirstDerivativeTranspose`](@ref) function.

    SAT_Periodic!(dest::AT,u::AT,cx::AT,cxy::AT,SP::SAT_Periodic{TN,:Curvilinear,TT}) where {TT,AT<:AbstractMatrix{TT},TN<:NodeType}

"""
function SAT_Periodic! end
function SAT_Periodic!(dest::VT,u::VT,c::VT,SP::SAT_Periodic{TN}) where {VT<:AbstractVector,TN<:NodeType{:Left}}
    # for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        dest[1]  += SP.τ₀(c) * (u[1] - u[end])
        # Symmeteriser
        L₁u = (u[1] - u[end])
        for i = 1:SP.order
            dest[i]    += SP.τ₁*c[1]*SP.D₁ᵀE₀[i]*L₁u
            #Neumann terms
            dest[1]    += SP.α₀ * (c[1]*SP.E₀D₁[i]*u[i] - c[end]*SP.EₙD₁[i]*u[end-SP.order+i])
        end
        # println(dest[1])
    # end
    dest
end
function SAT_Periodic!(dest::VT,u::VT,c::VT,SP::SAT_Periodic{TN}) where {VT<:AbstractVector,TN<:NodeType{:Right}}
    # for (S,U,K) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        # Dirichlet terms
        dest[end]+= SP.τ₀(c) * (u[end] - u[1])
        # Symmeteriser
        L₁u = (u[1] - u[end])
        for i = 1:SP.order
            dest[end-SP.order+i]  += SP.τ₁*c[end]*SP.D₁ᵀEₙ[i]*L₁u
            #Neumann terms
            dest[end]   += SP.α₀ * (c[1]*SP.E₀D₁[i]*u[i] - c[end]*SP.EₙD₁[i]*u[end-SP.order+i])
        end
    # end
    dest
end
function SAT_Periodic!(dest::AT,u::AT,c::AT,SP::SAT_Periodic{TN,:Cartesian,TT}) where {TT,AT<:AbstractMatrix{TT},TN<:NodeType}
    for (DEST,U,C) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(c))
        SAT_Periodic!(DEST,U,C,SP)
    end
    dest
end
function SAT_Periodic!(dest::AT,u::AT,cx::AT,cxy::AT,SP::SAT_Periodic{TN,:Curvilinear,TT}) where {TT,AT<:AbstractMatrix{TT},TN<:NodeType}
    for (DEST,U,C) in zip(SP.loopaxis(dest),SP.loopaxis(u),SP.loopaxis(cx))
        SAT_Periodic!(DEST,U,C,SP)
    end

    n = size(dest,SP.axis)
    m = size(dest,mod1(SP.axis+1,2))

    if SP.side == Left
        DEST = view(dest,   1,1:m)
        SRC = view(u[1,1:m]-u[n,1:m], :,:)
        C = view(cxy[1,1:m]-cxy[n,1:m],:,:)
    elseif SP.side == Right
        DEST = view(dest,   n,1:m)
        SRC = view(u[n,1:m]-u[1,1:m],:,:)
        C = view(cxy[n,1:m]-cxy[1,1:m],:,:)
    elseif SP.side == Up
        DEST = view(dest,   1:m,1)
        SRC = view(u[:,1]-u[:,n], :,1)
        C = view(cxy[:,1]-cxy[:,n],:,1)
    else
        DEST = view(dest,   1:m,n)
        SRC = view(-u[:,1]+u[:,n],:,1)
        C = view(-cxy[:,1]+cxy[:,n],:,1)
    end
    FirstDerivativeTranspose!(DEST,SRC,m,SP.Δy,SP.order,TT(1))

    dest
end











