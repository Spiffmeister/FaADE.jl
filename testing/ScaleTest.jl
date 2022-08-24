using LinearAlgebra
using Printf
using Plots
# pyplot()
# using GLMakie
using JLD2

using BenchmarkTools
# using ProfileView
using Profile
using PProf

cd("..")
using Distributed
# addprocs(1)
@everywhere using Interpolations
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
nx = 21
ny = 21

Δx = (𝒟x[2]-𝒟x[1])/(nx-1)
Δy = (𝒟y[2]-𝒟y[1])/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],length=nx))
y = collect(range(𝒟y[1],𝒟y[2],length=ny))

kx = zeros(Float64,nx,ny) .+ 1.0e-8
ky = zeros(Float64,nx,ny) .+ 1.0e-8


Δt = 1.0 * min(Δx^2,Δy^2)
# t_f = 1000Δt
t_f = 50.0
N = ceil(Int64,t_f/Δt)

u₀(x,y) = (x-0.5)/(0.68-0.5)

gx(t) = [0.0, 1.0] #Dirichlet
gy(t) = [0.0, 0.0] #Periodic

order = order_x = order_y = 2
method = :cgie

κ_para = 1.0
τ_para = -1.0
χₘₙ = 2.1e-3
params = plas_diff.SampleFields.H_params([χₘₙ/2., χₘₙ/3.],[2.0, 3.0],[1.0, 2.0])
@everywhere function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[1] = x[2] #p_1            qdot        θ
    χ[2] = -sum(p.ϵₘₙ .*(sin.(p.m*x[1] - p.n*t) .* p.m)) #q_1        pdot        ψ
end


println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)


gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)
H_x = SBP_operators.build_H(ny,order_x)
H_x = 1.0 ./H_x.^2
H_y = SBP_operators.build_H(nx,order_y)
H_y = 1.0 ./H_y.^2
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
SBP_operators.time_solver(rate,u₀,21,21,Δx,Δy,x,y,2Δt,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=Inf,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)










nx = 41
ny = 41

Δx = (𝒟x[2]-𝒟x[1])/(nx-1)
Δy = (𝒟y[2]-𝒟y[1])/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],length=nx))
y = collect(range(𝒟y[1],𝒟y[2],length=ny))

kx = zeros(Float64,nx,ny) .+ 1.0e-8
ky = zeros(Float64,nx,ny) .+ 1.0e-8

Δt = 1.0 * min(Δx^2,Δy^2)
# t_f = 1000Δt
t_f = 50.0
N = ceil(Int64,t_f/Δt)

gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)
H_x = SBP_operators.build_H(ny,order_x)
H_x = 1.0 ./H_x.^2
H_y = SBP_operators.build_H(nx,order_y)
H_y = 1.0 ./H_y.^2


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


println("(nx,ny)=(",nx,",",ny,")    Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)
@time SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=Inf,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

@time SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,Dirichlet,SBP_operators.Periodic,
    method=method,order_x=order,order_y=order,samplefactor=Inf,tol=1e-5,rtol=1e-10,penalty_fn=penalty_fn,adaptive=true)

