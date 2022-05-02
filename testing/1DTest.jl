using LinearAlgebra
using Printf
using Plots

# cd("..")
using Pkg
Pkg.activate(".")
using SBP_operators


###
function rate(uₓₓ,u,n,x,Δx,t,Δt,k;order=2)
    uₓₓ = Dₓₓ!(uₓₓ,u,k,n,Δx,order=order)
    return uₓₓ
end

###
𝒟 = [0.0,1.0]
n = 51
Δx = 𝒟[2]/(n-1)
x = collect(range(𝒟[1],𝒟[2],step=Δx))

k = zeros(Float64,n) .+ 1.0

Δt = 0.1 * Δx^2
t_f = 1000Δt
N = ceil(Int64,t_f/Δt)

u₀(x) = exp.(-(x.-0.5).^2 ./ 0.02)

g(t) = [0.0, 0.0]

order = 2
method = :euler

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


###
soln = SBP_operators.time_solver(rate,u₀,n,x,Δx,t_f,Δt,k,g,:Neumann,method=method,order=order)


###
anim = @animate for i=1:N
    plot(soln.x,soln.u[:,i],label="t=$(@sprintf("%.5f",i*Δt))")
end

gif(anim,"yes.gif",fps=50)