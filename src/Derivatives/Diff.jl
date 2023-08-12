

"""
    generate_SecondDerivative

Returns a function `Diff(uₓₓ,u,c...)` which computes the derivative in 1 or 2 dimensions
"""
function generate_SecondDerivative end
function generate_SecondDerivative(n::Int64,Δx::Float64,order::DerivativeOrder)
    DO = DerivativeOrder{order}()
    let n = n, Δx=Δx, order=DO
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



function Diff(cache,u,op::DerivativeOperator)
    
end


