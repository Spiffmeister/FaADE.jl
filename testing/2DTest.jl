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
    uₓₓ = zeros(Float64,nx,ny)
    uₓₓ = Dₓₓ(u,nx,ny,Δy,ky,dim=1,order=order_y) + Dₓₓ(u,nx,ny,Δx,kx,dim=2,order=order_x)
    return uₓₓ
end


###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]
nx = 51
ny = 51

Δx = 𝒟x[2]/(nx-1)
Δy = 𝒟y[2]/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],step=Δx))
y = collect(range(𝒟y[1],𝒟y[2],step=Δy))


kx = zeros(Float64,nx,ny) .+ 1.0
ky = zeros(Float64,nx,ny) .+ 1.0e-14

Δt = 1.00 * min(Δx^2,Δy^2)
t_f = 500Δt
N = ceil(Int64,t_f/Δt)

u₀(x,y) = exp(-((x-0.5)^2 + (y-0.5)^2) / 0.02)

gx(t) = [0.0, 0.0]
gy(t) = [0.0, 0.0]

order = 2
method = :cgie

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
@time u = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,:Periodic,:Dirichlet,method=method,order_x=order,order_y=order,samplefactor=5)

###

skip = 1
fps = 10
#=
@time anim = @animate for i=1:skip:size(u)[3]
    surface(u[:,:,i]',label="t=$(@sprintf("%.5f",i*Δt))",zlims=(-0.3,1.3),xlabel="x",ylabel="y",camera=(20+5*cos(100*π*i*Δt),50))
end

gif(anim,"yes.gif",fps=fps)


###
@time anim = @animate for i=1:skip:size(u)[3]
    plot(u[:,25,i],label="t=$(@sprintf("%.5f",i*Δt))",ylims=(0.0,1.0))
end

gif(anim,"yes2.gif",fps=fps)

anim = @animate for i = 1:1:N
    l = @layout [a b; c d]
    p = plot(u[1,:,i],layout=l,label="t=$(@sprintf("%.5f",i*Δt))")
    plot!(p[2],u[end,:,i])
    plot!(p[3],u[:,1,i])
    plot!(p[4],u[:,end,i])
end

gif(anim,"yes2.gif",fps=5)
=#

anim = @animate for i = 1:skip:size(u)[3]
    l = @layout [a{0.7w} [b; c]]
    p = surface(u[:,:,i],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0))
    plot!(p[2],u[25,:,i])
    plot!(p[3],u[:,25,i])
end

gif(anim,"yes.gif",fps=fps)
