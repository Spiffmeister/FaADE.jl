using LinearAlgebra
using Printf
using Plots
pyplot()

using Pkg
Pkg.activate(".")
using SBP_operators


###
function rate(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,k;order_x=2,order_y=2)
    uₓₓ = Dₓₓ(u,k,nx,ny,Δx,ky,order=order_x) + Dₓₓ(u,nx,ny,Δy,ky,dim=2,order=order_y)
    return uₓₓ
end


###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
n = 51

Δx = 𝒟x[2]/(n-1)
x = collect(range(𝒟x[1],𝒟x[2],step=Δx))




kx = zeros(Float64,n) .+ 1.0
ky = zeros(Float64,n) .+ 1.0

Δt = 1.0 * Δx^2
t_f = 1000Δt
N = ceil(Int64,t_f/Δt)

u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

g(t) = [0.0, 1.0]

order = 2
method = :cgie

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
soln = SBP_operators.time_solver(rate,u₀,n,x,Δx,t_f,Δt,k,g,:Dirichlet,method=method,order=order)


###
anim = @animate for i=1:N
    plot(soln.x,soln.u[:,i],label="t=$(@sprintf("%.5f",i*Δt))")
end

gif(anim,"yes.gif",fps=50)