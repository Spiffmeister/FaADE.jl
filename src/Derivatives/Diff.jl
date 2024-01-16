



function mul! end



# function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,1,DO}) where {AT,KT,TT,DO}
#     # D₂!(dest,u,K,D.nx,D.Δx,D.order,TT(0))
#     SecondDerivativeInternal!(dest,u,K,D.Δx,D.nx,D.order,TT(0))
#     if !D.xperiodic
#         SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,D.order,TT(0))
#         SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,D.order,TT(0))
#     else
#         SecondDerivativeBoundaryPeriodic!(dest,u,K,D.Δx,Left,D.order,TT(0))
#         SecondDerivativeBoundaryPeriodic!(dest,u,K,D.Δx,Right,D.order,TT(0))
#     end
#     dest
# end


# function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,2,DO,:Constant}) where {TT,AT<:AbstractArray{TT},KT,DO<:Val}

#     for (A,B,C) in zip(eachcol(dest),eachcol(u),eachcol(K[1]))
#         SecondDerivativeInternal!(A,B,C,D.Δx,D.nx,D.order,TT(0))
#         if !D.xperiodic
#             SecondDerivativeBoundary!(A,B,C,D.Δx,Left,D.order,TT(0))
#             SecondDerivativeBoundary!(A,B,C,D.Δx,Right,D.order,TT(0))
#         else
#             SecondDerivativeBoundaryPeriodic!(A,B,C,D.Δx,Left,D.order,TT(0))
#             SecondDerivativeBoundaryPeriodic!(A,B,C,D.Δx,Right,D.order,TT(0))
#         end
#     end
#     for (A,B,C) in zip(eachrow(dest),eachrow(u),eachrow(K[2]))
#         D₂!(A,B,C,D.ny,D.Δy,D.order,TT(1))
#     end

#     dest
# end



"""
    mul! for Diffusion problems
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
    if !D.periodic
        SecondDerivativeInternal!(dest,u,K,D.Δx,D.n,    Val(DO),α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Left,   Val(DO),α)
        SecondDerivativeBoundary!(dest,u,K,D.Δx,Right,  Val(DO),α)
        FirstDerivativeInternal!(dest,KDu,D.Δx,D.n,     Val(DO),TT(1))
        FirstDerivativeBoundary!(dest,KDu,D.Δx,Left,    Val(DO),TT(1))
        FirstDerivativeBoundary!(dest,KDu,D.Δx,Right,   Val(DO),TT(1))
    elseif D.periodic
        SecondDerivativePeriodic!(dest,u,K,D.Δx,        Val(DO),D.n,α)
        FirstDerivativePeriodic!(dest,KDu,u,D.Δx,       Val(DO),D.n,TT(1))
    end
    dest
end

"""
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
    Multidimensional version of variable coefficient
"""
function mul!(dest::AT,u::AT,c::KT,D::DiffusionOperatorND{TT,2,DO,:Variable},α::TT) where {TT,AT<:AbstractMatrix{TT},KT<:AbstractVector{AT},DO}
    cx = c[1]
    cy = c[2]
    cxy = c[3]
    cache = D.cache

    for (TMP,U) in zip(eachrow(cache),eachrow(u))
        FirstDerivativeInternal!(TMP,U,D.DO[2].Δx,D.DO[2].n,Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[2].Δx,Left,     Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[2].Δx,Right,    Val(DO),TT(0))
    end
    @. cache = cxy*cache
    for (DEST,U,Kx,TMP) in zip(eachcol(dest),eachcol(u),eachcol(cx),eachcol(cache))
        mul!(DEST,U,Kx,TMP,D.DO[1],TT(0))
    end

    for (TMP,U) in zip(eachcol(cache),eachcol(u))
        FirstDerivativeInternal!(TMP,U,D.DO[1].Δx,D.DO[1].n,Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[1].Δx,Left,     Val(DO),TT(0))
        FirstDerivativeBoundary!(TMP,U,D.DO[1].Δx,Right,    Val(DO),TT(0))
    end
    @. cache = cxy*cache
    for (DEST,U,Ky,TMP) in zip(eachrow(dest),eachrow(u),eachrow(cy),eachrow(cache))
        mul!(DEST,U,Ky,TMP,D.DO[2],TT(1))
    end
    dest
end

