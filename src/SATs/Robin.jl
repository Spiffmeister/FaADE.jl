
"""
    SAT_Robin

Simulatenous approximation term for Robin boundary conditions
`NodeType` is either `Left` or `Right`

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

    function SAT_Robin(RHS::F1,Δx::TT,side::TN,order::Int;α=TT(1),β=TT(1),τ=nothing,Δy=TT(0),coord=:Cartesian) where {TT,TN<:NodeType{SIDE,AX},F1} where {SIDE,AX}
    
        # H⁻¹B (α u + n β D₁ - g), n = ∓ 1

        check_boundary(side)

        loopaxis = SelectLoopDirection(AX)


        Hinv = _InverseMassMatrix(order,Δx,side)
        D₁ = _BoundaryDerivative(order,Δx,side)

        if SIDE == :Left
            H⁻¹E = -Hinv[1] # B
            H⁻¹ED₁ = Hinv[1]*D₁ # B (-D₁)
        elseif SIDE == :Right
            H⁻¹E = -Hinv[end] # 
            H⁻¹ED₁ = -Hinv[end]*D₁
        end


        if isnothing(τ)
            τ = TT(1)
        end


        new{TN,coord,TT,F1,typeof(loopaxis)}(side,AX,order,RHS,H⁻¹ED₁,H⁻¹E,Δx,τ,α,β,loopaxis,Δy,coord)
    end
end











function SAT_Robin_solution!(dest::AT,u::AT,c::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE}} where SIDE
    SIDE == :Left ? j = 1 : j = lastindex(dest)
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

function SAT_Robin_data!(dest::AT,data::AT,SR::SAT_Robin{TN}) where {AT,TN<:NodeType{SIDE}} where SIDE
    SIDE == :Left ? j = 1 : j = lastindex(dest)
    SIDE == :Left ? n = 1 : n = -1

    for (DEST,DATA) in zip(SR.loopaxis(dest),SR.loopaxis(data))
        DEST[j] -= SR.τ * SR.H⁻¹E * DATA[1]
    end
    dest
end
