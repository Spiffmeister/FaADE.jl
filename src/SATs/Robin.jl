
"""
    SAT_Robin{
            TN<:NodeType,
            COORD,
            TT<:Real,
            F1<:Function, LAT<:Function
        } <: SimultanousApproximationTerm{:Robin}

Simulatenous approximation term for Robin boundary conditions
    ``\\left.\\alpha u(x_i) + \\beta\\left.\\frac{\\partial u}{\\partial x}\\right|_{x_i} - g(t) = 0``

In Cartesian coordinates the SAT reads

``\\operatorname{SAT}_R = \\tau H^{-1} E (\\alpha u + \\kappa D_x u - g)``

In curvilinear coordinates the SAT is

``\\operatorname{SAT}_R = \\tau H^{-1}_q E (\\alpha u + K_q D_q u + K_{qr} D_r u - g)``

TODO: Testing
"""
struct SAT_Robin{
            TN<:NodeType,
            COORD,
            TT<:Real,
            F1<:Function, LAT<:Function
        } <: SimultanousApproximationTerm{:Robin}
    
    side        :: TN
    axis        :: Int
    order       :: Int
    RHS         :: F1
    H⁻¹ED₁      :: Vector{TT}
    H⁻¹E        :: TT
    Δx          :: TT
    τ           :: TT
    α           :: TT
    β           :: TT
    loopaxis    :: LAT
    Δy          :: TT
    coordinates :: Symbol

    @doc """
        SAT_Robin(RHS::F1,Δx::TT,side::TN,order::Int,Δy=TT(0),coord=:Cartesian;α=TT(1),β=TT(1),τ=nothing) where {TT,TN<:NodeType{SIDE,AX},F1} where {SIDE,AX}

    Constructor for Robin SAT
        
    Required arguments:
    - `RHS`: Right hand side function
    - `Δx`: Grid spacing
    - `side`: Boundary side, either `Left`, `Right`, `Up`, `Down`
    - `order`: Order of the ``D_x`` operators, either `2` or `4`
    """
    function SAT_Robin(RHS::F1,Δx::TT,side::TN,order::Int,Δy=TT(0),coord=:Cartesian;α=TT(1),β=TT(1),τ=nothing) where {TT,TN<:NodeType{SIDE,AX},F1} where {SIDE,AX}
    
        # H⁻¹B (α u + n β D₁ - g), n = ∓ 1

        check_boundary(side)

        loopaxis = SelectLoopDirection(AX)


        Hinv = _InverseMassMatrix(order,Δx,side)
        D₁ = _BoundaryDerivative(order,Δx,side)

        if SIDE == :Left
            H⁻¹E = Hinv[1] # B
            H⁻¹ED₁ = -Hinv[1]*D₁ # B (-D₁)
        elseif SIDE == :Right
            H⁻¹E = Hinv[end] # 
            H⁻¹ED₁ = Hinv[end]*D₁
        end


        if isnothing(τ)
            τ = -TT(1)
        end


        new{TN,coord,TT,F1,typeof(loopaxis)}(side,AX,order,RHS,H⁻¹ED₁,H⁻¹E,Δx,τ,α,β,loopaxis,Δy,coord)
    end
end









"""
Robin SAT method for implicit solvers. Applys portion of SAT related to the solution.

2D method which loops over 2D domains using the `loopaxis` value from [`SAT_Robin`](@ref).

    SAT_Robin_solution!(dest::AT,u::AT,c::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE,AX}} where {SIDE,AX}

The curvilinear call wraps around the Cartesian method. Calls the [`D₁!`](@ref) operator from [`FaADE.Derivatives`](@ref)

    SAT_Robin_solution!(dest::AT,u::AT,c::AT,cxy::AT,SR::SAT_Robin{TN,:Curvilinear,TT}) where {TT,AT,TN}

"""
function SAT_Robin_solution! end
function SAT_Robin_solution!(dest::AT,u::AT,c::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE,AX}} where {SIDE,AX}
    SIDE == :Left ? j = 1 : j = size(dest)[AX]
    SIDE == :Left ? m = 0 : m = j-SR.order
    SIDE == :Left ? n = 1 : n = -1

    for (DEST,U,C) in zip(SR.loopaxis(dest),SR.loopaxis(u),SR.loopaxis(c))
        DEST[j] += SR.τ * SR.α * SR.H⁻¹E * U[j]
        for i = 1:SR.order
            DEST[j] += SR.τ * SR.β * C[j] * SR.H⁻¹ED₁[i] * U[m+i]
        end
    end
    dest
end
function SAT_Robin_solution!(dest::AT,u::AT,c::AT,cxy::AT,SR::SAT_Robin{TN,:Curvilinear,TT}) where {TT,AT,TN}
    SAT_Robin_solution!(dest,u,c,SR)

    n = size(dest,SR.axis)
    m = size(dest,mod1(SR.axis+1,2))

    if SR.side == Left
        DEST =  view(dest,              1,1:m)
        SRC =   view(u,                 1,1:m)
        C =     view(SR.τ*SR.H⁻¹E*cxy,  1,1:m)
    elseif SR.side == Right
        DEST =  view(dest,              n,1:m)
        SRC =   view(u,                 1,1:m)
        C =     view(SR.τ*SR.H⁻¹E*cxy,  n,1:m)
    elseif SR.side == Down
        DEST =  view(dest,              1:m,1)
        SRC =   view(u,                 1:m,1)
        C =     view(SR.τ*SR.H⁻¹E*cxy,  1:m,1)
    elseif SR.side == Up
        DEST =  view(dest,              1:m,n)
        SRC =   view(u,                 1:m,1)
        C =     view(SR.τ*SR.H⁻¹E*cxy,  1:m,n)
    end

    D₁!(DEST,C,SRC,m,SR.Δy,Val(SR.order),TT(1))
    dest
end





"""
Robin SAT method for implicit solvers. Applies boundary data i.e. the portion of the SAT which operates on ``g`` only.

Only one method is used for both Cartesian and curvilinear coordinates.

    SAT_Robin_data!(dest::AT,data::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE,AX}} where {SIDE,AX}

"""
function SAT_Robin_data!(dest::AT,data::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE,AX}} where {SIDE,AX}
    SIDE == :Left ? j = 1 : j = size(dest)[AX]
    # SIDE == :Left ? n = 1 : n = -1

    for (DEST,DATA) in zip(SR.loopaxis(dest),SR.loopaxis(data))
        DEST[j] += -SR.τ * SR.H⁻¹E * DATA[1]
    end
    dest
end


