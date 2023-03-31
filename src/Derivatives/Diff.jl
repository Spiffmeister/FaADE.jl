

"""
    generate_SecondDerivative

Returns a function `Diff(uₓₓ,u,c...)` which computes the derivative in 1 or 2 dimensions
"""
function generate_SecondDerivative end
function generate_SecondDerivative(n::Int64,Δx::Float64,order::Int64)
    let n = n, Δx=Δx, order=order
        Diff(uₓₓ,u,c) = Dₓₓ!(uₓₓ,u,c,n,Δx,order)
        return Diff
    end
end
function generate_SecondDerivative(nx::Int64,ny::Int64,Δx::Float64,Δy::Float64,order::Int64)
    let nx=nx, ny=ny,
            Δx=Δx, Δy=Δy,
            order=order
        Diff(uₓₓ,u,cx,cy) = Dₓₓ!(uₓₓ,u,cx,cy, nx,ny,Δx,Δy,order)
        return Diff
    end
end
