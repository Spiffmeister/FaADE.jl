


"""
    mul!
"""
function mul! end



"""
    mul!(dest::VT,u::VT,K::VT,D::DiffusionOperator{TT,DO,:Constant},α) where {TT<:Real,VT<:AbstractVector{TT},DO}
"""
function mul!(dest::VT,u::VT,K::VT,D::DiffusionOperator{TT,DO,:Constant},α) where {TT<:Real,VT<:AbstractVector{TT},DO}
    if !D.periodic
        SecondDerivativeInternal!(dest,u,K,D.Δx,D.n,    Val(DO),α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,   Val(DO),α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,  Val(DO),α)
    elseif D.periodic
        SecondDerivativePeriodic!(dest,u,K,D.Δx,        Val(DO),D.n,α)
    end
    dest
end
function mul!(dest::VT,u::VT,K::VT,KDu::VT,D::DiffusionOperator{TT,DO,:Variable},α) where {TT<:Real,VT<:AbstractVector{TT},DO}
    order = Val(DO)
    if !D.periodic
        # D_{qq}^{K_q} u = D_q^T K_q D_q u
        SecondDerivativeInternal!(dest,u,K,D.Δx,D.n,    order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,   order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,  order,α)
        # D_q (K_{qr} D_r u)
        FirstDerivativeInternal!(dest,KDu,D.Δx,D.n,     order,TT(1))
        FirstDerivativeBoundary!(dest,KDu,D.Δx,Left,    order,TT(1))
        FirstDerivativeBoundary!(dest,KDu,D.Δx,Right,   order,TT(1))
    elseif D.periodic
        SecondDerivativePeriodic!(dest,u,K,D.Δx,        order,D.n,α)
        FirstDerivativePeriodic!(dest,KDu,u,D.Δx,       order,D.n,TT(1))
    end
    dest
end

"""
    mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Constant},α::TT) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO}
Multidimensional version of constant coefficient
"""
function mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Constant},α::TT) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO}
    cx = c[1]
    cy = c[2]

    for (DEST,U,K) in zip(eachcol(dest),eachcol(u),eachcol(cx))
        mul!(DEST,U,K,D.DO[1],TT(0))
    end
    for (DEST,U,K) in zip(eachrow(dest),eachrow(u),eachrow(cy))
        mul!(DEST,U,K,D.DO[2],TT(1))
    end
    dest
end
"""
    mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Variable},α::TT) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO}
Multidimensional version of variable coefficient
"""
function mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Variable},α::TT) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO}
    cx = c[1]
    cy = c[2]
    cxy = c[3]
    cache = D.cache

    for (TMP,U) in zip(eachrow(cache),eachrow(u))
        # D_r u
        FirstDerivativeInternal!(TMP,U,D.DO[2].Δx,D.DO[2].n,Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[2].Δx,Left,     Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[2].Δx,Right,    Val(DO),TT(0))
    end
    @. cache = cxy*cache # K_{qr} (D_r u)
    for (DEST,U,Kx,TMP) in zip(eachcol(dest),eachcol(u),eachcol(cx),eachcol(cache))
        # K_qr D_r u + D_{qq}^{K_q} u 
        mul!(DEST,U,Kx,TMP,D.DO[1],TT(0))
    end

    for (TMP,U) in zip(eachcol(cache),eachcol(u))
        # D_q u
        FirstDerivativeInternal!(TMP,U,D.DO[1].Δx,D.DO[1].n,Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[1].Δx,Left,     Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[1].Δx,Right,    Val(DO),TT(0))
    end
    @. cache = cxy*cache # K_{qr} (D_q u)
    for (DEST,U,Ky,TMP) in zip(eachrow(dest),eachrow(u),eachrow(cy),eachrow(cache))
        # K_{qr} D_q u + D_{rr}^{K_r} u
        mul!(DEST,U,Ky,TMP,D.DO[2],TT(1))
    end
    dest
end

