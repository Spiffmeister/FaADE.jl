

"""
    generate_SecondDerivative

Returns a function `Diff(uₓₓ,u,c...)` which computes the derivative in 1 or 2 dimensions
"""
function generate_SecondDerivative end
function generate_SecondDerivative(n::Int64,Δx::Float64,order::DerivativeOrder)
    # DO = DerivativeOrder{order}()
    let n = n, Δx=Δx, order=order
        Diff!(uₓₓ,u,c) = D₂!(uₓₓ,u,c,n,Δx,order,0.0)
        return Diff!
    end
end
function generate_SecondDerivative(nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,order::Integer)
    DO = DerivativeOrder{order}()
    let nx=nx, ny=ny,
            Δx=Δx, Δy=Δy,
            DO=DO
        
        Diff!(uₓₓ,u,cx,cy) = D₂!(uₓₓ,u,cx,cy, nx,ny,Δx,Δy,DO,DO,0.0)
        return Diff!
    end
end



"""
    Diff!(cache,u,K,D::DerivativeOperator)
Apply the derivative operator to `u` given the diffusion coefficient `K` and store result in `cache`
"""
function Diff!(cache::AT,u::AT,K::KT,D::DerivativeOperator) where {AT,KT}
    for (C,U) in zip(cache,u,K)
        D(C,U,K)
    end
    cache
end


# function Diff!(cache::AT,u::AT,K::KT,D::DerivativeOperator{TT,1,ORD}) where {AT,KT,TT,ORD}
#      (cache,,U,K)
#     cache
# end

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
        # D₂!(A,B,C,D.nx,D.Δx,D.order,TT(0))
    end
    for (A,B,C) in zip(eachrow(dest),eachrow(u),eachrow(K[2]))
        D₂!(A,B,C,D.ny,D.Δy,D.order,TT(1))
    end

    dest
end

function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,2,DO,:Variable}) where {TT,AT<:AbstractArray{TT},KT,DO<:DerivativeOrder}

    for (A,B,Kx,DKy) in zip(eachcol(dest),eachcol(u),eachcol(K[1]),eachcol(K[4]))
        D₂!(A,B,Kx,D.nx,D.Δx,D.order,TT(0))
        D₁!(A,DKy,DKy,D.nx,D.Δx,D.order,TT(1))
    end
    for (A,B,Ky,DKx) in zip(eachrow(dest),eachrow(u),eachrow(K[2]),eachrow(K[3]))
        D₂!(A,B,Ky,D.ny,D.Δy,D.order,TT(1))
        D₁!(A,DKx,DKx,D.ny,D.Δy,D.order,TT(1))
    end

    dest
end



# function MixedOp(dest,u,cx,Dcy,Δx,Δy,DO,α)
#     local cache :: TT
#     for i = m:nx-m+1
#         cache = SecondDerivativeInternal(u,cx,Δx,DO,i)
#         cache += Dcy[i]*FirstDerivativeInternal(u,Δy,DO,i)
#         dest[i] = α*dest[i] + cache
#     end
# end