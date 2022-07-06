using LinearAlgebra
using Printf
using Plots
pyplot()

# using BenchmarkTools

cd("..")
push!(LOAD_PATH,"./SBP_operators")
push!(LOAD_PATH,"./plas_diff")
using SBP_operators
using plas_diff


###
function rate(uₓₓ,u,nx,ny,x,y,Δx,Δy,t,Δt,kx,ky;order_x=2,order_y=2)
    uₓₓ = Dₓₓ(u,nx,ny,Δx,kx,dim=1,order=order_x) + Dₓₓ(u,nx,ny,Δy,ky,dim=2,order=order_y)
    return uₓₓ
end




###
𝒟x = [0.0,1.0]
𝒟y = [0.0,2π]
nx = 31
ny = 31

Δx = 𝒟x[2]/(nx-1)
Δy = 𝒟y[2]/(ny-1)
x = collect(range(𝒟x[1],𝒟x[2],step=Δx))
y = collect(range(𝒟y[1],𝒟y[2],step=Δy))
if y[end] != 2π
    y = vcat(y,2π)
end

kx = zeros(Float64,nx,ny) .+ 1.0e-2
ky = zeros(Float64,nx,ny) .+ 1.0e-2

Δt = 1.0 * min(Δx^2,Δy^2)
t_f = 400Δt
N = ceil(Int64,t_f/Δt)

u₀(x,y) = exp(-(x-0.5)^2/0.02 - (y-π)^2/0.5)

gx(t) = [0.0, 0.0]
gy(t) = [0.0, 0.0]

order = 2
method = :cgie

println("Δx=",Δx,"      ","Δt=",Δt,"        ","final time=",t_f)




params = plas_diff.SampleFields.H_params([0.],[0.],[0.])
gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,plas_diff.SampleFields.χ_h!,params)

H_x = SBP_operators.build_H(nx,2)
H_x = 1.0 ./H_x.^2
H_x = diagm(H_x)

H_y = SBP_operators.build_H(nx,2)
H_y = 1.0 ./H_y.^2
H_y = diagm(H_y)

Hinv = diag(kron(I(nx),H_y) + kron(H_x,I(ny)))

κ_para = 1.0
τ_para = -1.0




### Parallel Penalty ###
function penalty_fn(u,uₒ)
    uₚ = zeros(Float64,nx,ny)
    for i = 1:nx
        for j = 1:ny
            uₚ[i,j] = κ_para * τ_para/2.0 * (Hinv[(i-1)*j+1]) * (u[i,j] - uₒ[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]])
            uₚ[i,j] += κ_para * τ_para/2.0 * (Hinv[(i-1)*j+1]) * (u[i,j] - uₒ[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]])
        end
    end
    return uₚ
end

function penalty_cg_fn_inner(u)
    uₚ = zeros(Float64,nx,ny)
    for i = 1:nx
        for j = 1:ny
            uₚ[i,j] = κ_para * τ_para/2.0 * Hinv[(i-1)*j+1] * u[i,j]
            uₚ[i,j] += κ_para * τ_para/2.0 * Hinv[(i-1)*j+1] * u[i,j]
        end
    end
    return uₚ
end
function penalty_cg_fn_outer(u)
    uₚ = zeros(Float64,nx,ny)
    for i = 1:nx
        for j = 1:ny
            uₚ[i,j] = κ_para * τ_para/2.0 * Hinv[(i-1)*j+1] * -u[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]]
            uₚ[i,j] += κ_para * τ_para/2.0 * Hinv[(i-1)*j+1] * -u[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]]
        end
    end
    return uₚ
end




###
@time u = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,:Dirichlet,:Periodic,
    # method=method,order_x=order,order_y=order,samplefactor=1,tol=1e-14,penalty_fn=penalty_cg_fn_outer,penalty_fn_outer=penalty_cg_fn_outer)
    method=method,order_x=order,order_y=order,samplefactor=1,tol=1e-14,penalty_fn=penalty_fn)

###

# plas_diff.plot_grid(gdata)
# savefig("yes2.png")

skip = 5
fps = 25

energy = zeros(N)
for i = 1:N
    energy[i] = norm(u[:,:,i],2)
end

anim = @animate for i = 1:skip:size(u)[3]
    l = @layout [a{0.7w} [b; c]]
    p = surface(u[:,:,i],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0),xlabel="y",ylabel="x")
    plot!(p[2],u[floor(Int64,nx/2),:,i])
    # plot!(p[3],u[:,floor(Int64,ny/2),i])
    plot!(p[3],energy[1:i])
end

gif(anim,"yes.gif",fps=fps)
