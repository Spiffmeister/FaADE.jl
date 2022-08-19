using LinearAlgebra
using Printf
using Plots
# pyplot()
# using GLMakie
using Interpolations
using JLD2

using BenchmarkTools
# using ProfileView
using Profile
using PProf

cd("..")
using Distributed
# addprocs(1)
@everywhere push!(LOAD_PATH,"./plas_diff")
@everywhere push!(LOAD_PATH,"./SBP_operators")
# @everywhere push!(LOAD_PATH,".")
@everywhere using SBP_operators
@everywhere using plas_diff
using SharedArrays

###
function rate(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky;order_x=2,order_y=2)
    # uₓₓ = Dₓₓ(u,nx,ny,Δx,kx,dim=1,order=order_x) + Dₓₓ(u,nx,ny,Δy,ky,dim=2,order=order_y)
    uₓₓ = Dₓₓ!(uₓₓ,u,nx,ny,Δx,Δy,kx,ky,order_x=order_x,order_y=order_y)
    return uₓₓ
end




###
𝒟x = [0.5,0.68]
𝒟y = [-π,π]
nx = 41
ny = 41

Δx = (𝒟x[2]-𝒟x[1])/(nx-1)
Δy = (𝒟y[2]-𝒟y[1])/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],length=nx))
y = collect(range(𝒟y[1],𝒟y[2],length=ny))

kx = zeros(Float64,nx,ny) .+ 1.0e-8
ky = zeros(Float64,nx,ny) .+ 1.0e-8


Δt = 1.0 * min(Δx^2,Δy^2)
# t_f = 200Δt
t_f = 100.0
N = ceil(Int64,t_f/Δt)

# u₀(x,y) = exp(-(x-0.5)^2/0.02 - (y-π)^2/0.5)
# u₀(x,y) = 0.5sin(4*2π*x) + 0.5sin(4*y)
u₀(x,y) = (x-0.5)/(0.68-0.5)

gx(t) = [0.0, 1.0] #Dirichlet
gy(t) = [0.0, 0.0] #Periodic

order = order_x = order_y = 2
method = :cgie

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)




# params = plas_diff.SampleFields.H_params([0.],[0.],[0.])
χₘₙ = 2.1e-3
params = plas_diff.SampleFields.H_params([χₘₙ/2., χₘₙ/3.],[2.0, 3.0],[1.0, 2.0])

@everywhere function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[1] = x[2] #p_1            qdot        θ
    χ[2] = -sum(p.ϵₘₙ .*(sin.(p.m*x[1] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)

H_x = SBP_operators.build_H(ny,order_x)
H_x = 1.0 ./H_x.^2
# H_x = diagm(H_x)

H_y = SBP_operators.build_H(nx,order_y)
H_y = 1.0 ./H_y.^2
# H_y = diagm(H_y)

# Hinv = diag(kron(I(nx),H_y) + kron(H_x,I(ny)))

κ_para = 1.0
τ_para = -1.0




### Parallel Penalty ###
function penalty_fn(u,uₒ,Δt)
    uₚ = SharedArray(zeros(Float64,nx,ny))
    umw = zeros(Float64,nx,ny)

    interp = LinearInterpolation((x,y),u)

    @sync @distributed for j = 1:ny
        for i = 1:nx

            umw[i,j] = 2u[i,j] - (u[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]] + u[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]])


            if 𝒟x[1] ≥ gdata.z_planes[1].x[i,j]
                w_f  = 0.0
            elseif 𝒟x[2] ≤ gdata.z_planes[1].x[i,j]
                w_f = 1.0
            else
                w_f = interp(gdata.z_planes[1].x[i,j],gdata.z_planes[1].y[i,j])
            end

            if 𝒟x[1] ≥ gdata.z_planes[2].x[i,j]
                w_b  = 0.0
            elseif 𝒟x[2] ≤ gdata.z_planes[2].x[i,j]
                w_b = 1.0
            else
                w_b = interp(gdata.z_planes[2].x[i,j],gdata.z_planes[2].y[i,j])
            end

            uₚ[i,j] = 1.0/(1.0 - κ_para * τ_para/2.0 * Δt * (H_y[i] + H_x[j])) *
                (u[i,j] - Δt*τ_para/4.0 *(H_y[i] + H_x[j])*(w_f + w_b))


        end
    end
    return uₚ, norm(umw)
end

"""
@time SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,0.1,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=1.0,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

###
@time soln,umw = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=1.0,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

###
"""

SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,2Δt,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=1.0,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

###
@benchmark SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=1000.0,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

###

# soln,_ = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
#     method=method,order_x=order,order_y=order,samplefactor=1.0,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)


# println("plotting")

# pdata = plas_diff.poincare(plas_diff.SampleFields.χ_h!,params,N_trajs=1000,N_orbs=100,x=𝒟x,y=𝒟y)


# plas_diff.plot_grid(gdata)

# N = length(soln.u)
# skip = 1
# fps = 1

# energy = zeros(N)
# maxerry = zeros(N)
# maxerrx = zeros(N)
# for i = 1:N
#     energy[i] = norm(soln.u[i][:,:],2)
#     maxerry[i] = norm(soln.u[i][:,1]-soln.u[i][:,end],Inf)
#     maxerrx[i] = norm(soln.u[i][1,:]-soln.u[i][end,:],Inf)
# end


# anim = @animate for i = 1:skip:N
#     l = @layout [a{0.7w} [b; c]]
#     p = surface(soln.u[i][:,:],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0),clims=(0.0,1.0),xlabel="y",ylabel="x",camera=(30,30))
#     plot!(p[2],soln.t[1:i],maxerry[1:i],ylims=(0.0,max(maximum(maxerrx),maximum(maxerry))),label="y_0 - y_N")
#     plot!(p[2],soln.t[1:i],maxerrx[1:i],label="x_0 - x_N")
#     # plot!(p[2],u[15,:,i],ylabel="u(x=0.5)")
#     plot!(p[3],soln.t[1:i],energy[1:i],ylabel="||u||_2")
# end
# gif(anim,"yes.gif",fps=fps)


# Slice info
# anim = @animate for i = 1:skip:N
#     l = @layout [a b; c d]
#     # p = surface(u[:,:,i],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0),clims=(0.0,1.0),xlabel="y",ylabel="x")
#     # plot!(p[2],maxerr[1:i],ylims=(0.0,maximum(maxerr)),ylabel="y_0 - y_N")
#     p = plot(u[13,:,i],layout=l,ylabel="u(x_13)")
#     plot!(p[2],u[14,:,i],ylabel="u(x_14)")
#     plot!(p[3],u[15,:,i],ylabel="u(x_15)")
#     plot!(p[4],u[16,:,i],ylabel="u(x_16)")
# end
# gif(anim,"yes2.gif",fps=fps)


# Mid region animation
# anim = @animate for i = 1:skip:N
#     # l = @layout [a b; c d]
#     p = surface(u[10:20,:,i],title="u[10-20,:,i]",zlims=(0.0,1.0),clims=(0.0,1.0),xlabel="y",ylabel="x",camera=(30,30))

#     # plot!(p[2],maxerr[1:i],ylims=(0.0,maximum(maxerr)),ylabel="y_0 - y_N")
# end
# gif(anim,"yes3.gif",fps=fps)


# Boundaries matching
# anim = @animate for i = 1:skip:N
#     l = @layout [a b ; c d]
#     p = plot(u[1,:,i],layout=l,ylims=(0.0,1.0),label="x_0")
#     plot!(p[2],u[end,:,i],ylims=(0.0,1.0),label="x_N")
#     plot!(p[3],u[:,1,i],ylims=(0.0,1.0),label="y_0")
#     plot!(p[4],u[:,end,i],ylims=(0.0,1.0),label="y_N")
# end
# gif(anim,"yes4.gif",fps=fps)


# ψ = repeat(gdata.x,1,gdata.ny)
# θ = repeat(gdata.y',gdata.nx,1)
# anim = @animate for i = 1:skip:N
#     # l = @layout [a b; c d]
#     p = contour(x,y,u[:,:,i],xlabel="y",ylabel="x",clims=(0.0,1.0))
#     scatter!(ψ[:],θ[:])
#     # plot!(p[2],maxerr[1:i],ylims=(0.0,maximum(maxerr)),ylabel="y_0 - y_N")
# end
# gif(anim,"yes5.gif",fps=fps)



# println("saving")

# if !(typeof(pdata.poincare) <: Nothing)
#     pdata.poincare = nothing
# end
# save_object("testrun.jld2",(soln,gdata,umw,pdata))

