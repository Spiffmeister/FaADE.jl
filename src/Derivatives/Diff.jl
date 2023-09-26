

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
        
        Diff!(uₓₓ,u,cx,cy) = D₂!(uₓₓ,u,cx,cy, nx,ny,Δx,Δy,DO,DO)
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
    D₂!(dest,u,K,D.nx,D.Δx,D.order,TT(0))
    dest
end
function mul!(dest::AT,u::AT,K::KT,D::DerivativeOperator{TT,2,DO}) where {AT,KT,TT,DO}
    D₂!(dest,u,K[1],K[2],D.nx,D.ny,D.Δx,D.Δy,D.order,D.order)
    dest
end
