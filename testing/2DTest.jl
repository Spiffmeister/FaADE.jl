using LinearAlgebra
using Printf
using Plots
# pyplot()

using BenchmarkTools

using Pkg
Pkg.activate(".")
using SBP_operators


###
function rate(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky;order_x=2,order_y=2)
    # uₓₓ = Dₓₓ(u,nx,ny,Δx,kx,dim=1,order=order_x) + Dₓₓ(u,nx,ny,Δy,ky,dim=2,order=order_y)
    uₓₓ = Dₓₓ!(uₓₓ,u,nx,ny,Δx,Δy,kx,ky,order_x=order_x,order_y=order_y)
    return uₓₓ
end


###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 41
ny = 41

Δx = 𝒟x[2]/(nx-1)
Δy = 𝒟y[2]/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],length=nx))
y = collect(range(𝒟y[1],𝒟y[2],length=ny))


kx = zeros(Float64,nx,ny) .+ 1.0
ky = zeros(Float64,nx,ny) .+ 1.0

Δt = 0.01*1.00 * min(Δx^2,Δy^2)
t_f = 1000Δt
N = ceil(Int64,t_f/Δt)

u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)
# u₀(x,y) = (x-0.5)/(0.68-0.5)

gx(t) = [0.0, 0.0]
gy(t) = [0.0, 0.0]

order = 2
method = :cgie

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
soln, = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,Periodic,
    method=method,order_x=order,order_y=order,rtol=1e-14,adaptive=false)


# ky = zeros(Float64,nx,ny)
@benchmark SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,Periodic,
method=method,order_x=order,order_y=order,rtol=1e-14,adaptive=false) seconds = 30




###

# skip = 5
# fps = 50

# anim = @animate for i = 1:length(soln.u)
#     l = @layout [a{0.7w} [b; c]]
#     p = surface(soln.u[i],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0),xlabel="y",ylabel="x")
#     plot!(p[2],soln.u[i][25,:],ylims=(-0.05,1.05))
#     plot!(p[3],soln.u[i][:,25],ylims=(-0.05,1.05))
# end

# gif(anim,"yes.gif",fps=fps)


# @benchmark SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Neumann,SBP_operators.Periodic,
#     method=method,order_x=order,order_y=order,samplefactor=Inf,tol=1e-5,rtol=1e-10,adaptive=false)