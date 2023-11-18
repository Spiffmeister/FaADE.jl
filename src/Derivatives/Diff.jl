



function mul! end



function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,1,DO}) where {AT,KT,TT,DO}
    # D₂!(dest,u,K,D.nx,D.Δx,D.order,TT(0))
    SecondDerivativeInternal!(dest,u,K,D.Δx,D.nx,D.order,TT(0))
    if !D.xperiodic
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,D.order,TT(0))
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,D.order,TT(0))
    else
        SecondDerivativeBoundaryPeriodic!(dest,u,K,D.Δx,Left,D.order,TT(0))
        SecondDerivativeBoundaryPeriodic!(dest,u,K,D.Δx,Right,D.order,TT(0))
    end
    dest
end


function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,2,DO,:Constant}) where {TT,AT<:AbstractArray{TT},KT,DO<:DerivativeOrder}

    for (A,B,C) in zip(eachcol(dest),eachcol(u),eachcol(K[1]))
        SecondDerivativeInternal!(A,B,C,D.Δx,D.nx,D.order,TT(0))
        if !D.xperiodic
            SecondDerivativeBoundary!(A,B,C,D.Δx,Left,D.order,TT(0))
            SecondDerivativeBoundary!(A,B,C,D.Δx,Right,D.order,TT(0))
        else
            SecondDerivativeBoundaryPeriodic!(A,B,C,D.Δx,Left,D.order,TT(0))
            SecondDerivativeBoundaryPeriodic!(A,B,C,D.Δx,Right,D.order,TT(0))
        end
    end
    for (A,B,C) in zip(eachrow(dest),eachrow(u),eachrow(K[2]))
        D₂!(A,B,C,D.ny,D.Δy,D.order,TT(1))
    end

    dest
end



"""
    mul! for Diffusion problems
"""
function mul!(dest::VT,u::VT,K::VT,D::DiffusionOperator{TT,DO,:Constant},α) where {TT<:Real,VT<:AbstractVector{TT},DO<:DerivativeOrder}
    if !D.periodic
        SecondDerivativeInternal!(dest,u,K,D.Δx,D.n,D.order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,D.order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,D.order,α)
    elseif D.periodic
        SecondDerivativePeriodic!(dest,u,K,D.Δx,D.order,D.n,α)
    end
    dest
end
function mul!(dest::VT,u::VT,K::VT,Kxy::VT,D::DiffusionOperator{TT,DO,:Variable},α) where {TT<:Real,VT<:AbstractVector{TT},DO<:DerivativeOrder}
    if !D.periodic
        SecondDerivativeInternal!(dest,u,K,D.Δx,D.n,D.order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,D.order,α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,D.order,α)
        FirstDerivativeInternal!(dest,Kxy,u,D.Δx,D.n,D.order,TT(1))
        FirstDerivativeBoundary!(dest,Kxy,u,D.Δx,Left,D.order,TT(1))
        FirstDerivativeBoundary!(dest,Kxy,u,D.Δx,Right,D.order,TT(1))
    elseif D.periodic
        SecondDerivativePeriodic!(dest,u,K,D.Δx,D.order,D.n,α)
        FirstDerivativePeriodic!(dest,Kxy,u,D.Δx,D.order,D.n,TT(1))
    end
    dest
end

"""
    Multidimensional version of constant coefficient
"""
function mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Constant}) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO<:DerivativeOrder}
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
    Multidimensional version of variable coefficient
"""
function mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Variable},α) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO<:DerivativeOrder}
    cx = c[1]
    cy = c[2]
    cxy = c[3]

    for (DEST,U,Kx,Kxy) in zip(eachcol(dest),eachcol(u),eachcol(cx),eachcol(cxy))
        mul!(DEST,U,Kx,Kxy,D.DO[1],α)
    end
    for (DEST,U,Ky,Kxy) in zip(eachrow(dest),eachrow(u),eachrow(cy),eachrow(cxy))
        mul!(DEST,U,Ky,Kxy,D.DO[2],TT(1))
    end
    dest
end

