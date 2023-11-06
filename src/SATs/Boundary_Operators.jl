"""
    BoundaryDerivativeTranspose

Used to construct `EH⁻¹D₁ᵀ` where `E=E₀` or `Eₙ`.
"""
function BoundaryDerivativeTranspose end
function BoundaryDerivativeTranspose(::NodeType{:Left},order::Int,Δx::Float64)
    if order == 2
        return [-2.0,1.0]/Δx
    elseif order == 4
        #[-24.0/17.0 * 48.0/17.0, 59.0/34.0 * 48.0/59.0, -4.0/17.0 * 48.0/43.0, -3.0/34.0 * 48.0/49.0]/Δx
        return [-1152.0/289.0, 24.0/17.0, -192.0/731.0, -72.0/833.0]/Δx #H⁻¹Dᵀ
    elseif order == 6
        return [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]/Δx
    end
end
function BoundaryDerivativeTranspose(::NodeType{:Right},order::Int,Δx::Float64)
    if order == 2
        return [-1.0,2.0]/Δx
    elseif order == 4
        # return -[-3.0/34.0 * 48.0/49.0, -4.0/17.0 * 48.0/43.0, 59.0/34.0 * 48.0/59.0, -24.0/17.0 * 48.0/17.0]/Δx
        return [72.0/833.0, 192.0/731.0, -24.0/17.0, 1152.0/289.0]/Δx #H⁻¹Dᵀ
    elseif order == 6
        return [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]/Δx
    end
end


function _DerivativeTranspose(order::Int,Δx::TT,::NodeType{TN}) where {TT,TN}
    if order == 2
        return [-1,1]/Δx
    elseif order == 4
        DT = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]/Δx
        return DT
    elseif order == 6
        error("order 6 not set")
    end

    # if TN == :Left
    #     return DT
    # else
    #     return reverse(DT)
    # end
end
function _InverseMassMatrix(order::Int,Δx::TT,::NodeType{TN}) where {TT,TN}
    if order == 2
        Hinv = inv.([1/2, 1]*Δx) 
    elseif order == 4
        Hinv = inv.([17/48, 59/48, 43/48, 49/48]*Δx)
    elseif order == 6
        error("order 6 not set")
    end
    
    if TN == :Left
        return Hinv
    else
        return reverse(Hinv)
    end
end
_BoundaryOperator(TT,::NodeType{:Left}) = TT(-1)
_BoundaryOperator(TT,::NodeType{:Right}) = TT(1)


function _BoundaryDerivative(order::Int,Δx::TT,::NodeType{TN}) where {TT,TN}
    if order == 2
        return [-1.0,1.0]/Δx
    elseif order == 4
        D = [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]/Δx
    elseif order == 6
        D = [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]/Δx
    end

    if TN == :Left
        return D
    elseif TN == :Right
        return reverse(-D)
    end
end

"""
    BoundaryDerivative

Used to construct `E₀D₁` or `EₙD₁`.
"""
function BoundaryDerivative end
function BoundaryDerivative(::NodeType{:Left},Δx::Real,order::Int)
    if order == 2
        return [-1.0,1.0]/Δx
    elseif order == 4
        return [-24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0]/Δx
    elseif order == 6
        return [-1.582533518939116, 2.033378678700676, -0.141512858744873, -0.450398306578272, 0.104488069284042, 0.036577936277544]/Δx
    end
end
function BoundaryDerivative(::NodeType{:Right},Δx::Real,order::Int)
    if order == 2
        return [-1.0,1.0]/Δx
    elseif order == 4
        return [3.0/34.0, 4.0/17.0, -59.0/34.0, 24.0/17.0]/Δx
    elseif order == 6
        return [-0.036577936277544, -0.104488069284042, 0.450398306578272, 0.141512858744873, -2.033378678700676, 1.582533518939116]/Δx
    end
end





"""
    SATpenalties

Determines the penatly parameters for the given boundary conditions.

TODO: Remove inlineing
"""
function SATpenalties end
@inline function SATpenalties(::BoundaryConditionType{:Robin},a,Δx::Float64,order::Int64)
    h = hval(order)

    τ = 1.0/(a * h * Δx) # τ=1/a H^{-1}
    return τ
end
@inline function SATpenalties(::BoundaryConditionType{:Periodic},Δx::Real,order::Int64)
    h = hval(order)

    α₀ = 0.5/(h*Δx)
    τ₁ = -0.5/Δx #1/h accounted for in BoundaryDerivativeTranspose
    τ₀(c) = -(1.0+1.0/max(c[1],c[end]))*max(c[1],c[end])/(2*(h*Δx)^2)

    return α₀,τ₁,τ₀
end
@inline function SATpenalties(::BoundaryConditionType{:Interface},Δx⁺,Δx⁻,order⁺,order⁻)

    h⁻ = hval(order⁻)
    h⁺ = hval(order⁺)

    α₀ = -0.5/(h⁻*Δx⁻)
    τ₁ = 0.5/(Δx⁺)
    τ₀(c) = max(c/2(h⁻*Δx⁻),c/2(h⁺*Δx⁺))

    return α₀,τ₁,τ₀
end






"""
    hval(order::Int64)

Returns the value of ``h^{-1}`` for the penalties

TODO: Remove inlineing
"""
@inline function hval(order::Int64)
    if order == 2
        return 0.5
    elseif order == 4
        return 17.0/48.0
    elseif order == 6
        return 13649.0/43200.0
    end
end


