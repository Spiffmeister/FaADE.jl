using LinearAlgebra
using Revise
using FaADE





u0(x,y) = x
soln(x,y) = 1.0


order = 2

error = []

for n in [21,41,81,161]
    
    Dom = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
        v->[0.0, v/2 - 0.25],
        v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
            (1-v)*[0.0,0.25] + v*[0.0,0.25],
        u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
        n,n)
    
    
    
    Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx, Dom.Δx, order, false, :variable)
    Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny, Dom.Δy, order, false, :variable)


    u = [u0(Dom[i,j]...) for i in 1:Dom.nx, j in 1:Dom.ny]
    
    exact_du = [soln(Dom[i,j]...) for i in 1:Dom.nx, j in 1:Dom.ny]
    
    diff_u = zeros(size(Dom))

    for (UX, U, C) in zip(eachcol(diff_u), eachcol(u), eachcol(Dom.qx))
        D₁!(UX, C, U, Dx.n, Dx.Δx, Dx.order, 0.0)
    end
    for (UX, U, C) in zip(eachrow(diff_u), eachrow(u), eachrow(Dom.rx))
        D₁!(UX, C, U, Dy.n, Dx.Δx, Dy.order, 1.0)
    end
    
    @show tmp = norm(diff_u .- exact_du)/norm(exact_du)
    push!(error,tmp)
    
end


error[1:end-1] ./ error[2:end]