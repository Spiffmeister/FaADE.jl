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

kx = zeros(Float64,nx,ny) .+ 1.0e-8
ky = zeros(Float64,nx,ny) .+ 1.0e-8

Δt = 1.0 * min(Δx^2,Δy^2)
t_f = 200Δt
# t_f = 10.0
N = ceil(Int64,t_f/Δt)

# u₀(x,y) = exp(-(x-0.5)^2/0.02 - (y-π)^2/0.5)
u₀(x,y) = x

gx(t) = [0.0, 1.0] #Dirichlet
gy(t) = [0.0, 0.0] #Periodic

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
function penalty_fn(u,uₒ,Δt)
    uₚ = zeros(Float64,nx,ny)
    for i = 1:nx
        for j = 1:ny
            # uₚ[i,j] = κ_para * τ_para/2.0 * (Hinv[(i-1)*j+1]) * (u[i,j] - uₒ[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]])
            # uₚ[i,j] += κ_para * τ_para/2.0 * (Hinv[(i-1)*j+1]) * (u[i,j] - uₒ[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]])

            # uₚ[i,j] = u[i,j] + Δt * κ_para * τ_para/2.0 * (Hinv[(i-1)*j+1]) * 
                # (2.0*u[i,j] - (uₒ[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]] + uₒ[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]]))

            uₚ[i,j] = 1.0/(1.0 - 2.0 * κ_para * τ_para/2.0 * Δt * Hinv[(i-1)*j+1]) *
                (u[i,j] - Δt*τ_para/2.0 *Hinv[(i-1)*j+1]*(u[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]] + u[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]]))

        end
    end
    return uₚ
end




###
@time soln = SBP_operators.time_solver(rate,u₀,nx,ny,Δx,Δy,x,y,t_f,Δt,kx,ky,gx,gy,:Dirichlet,:Periodic,
    method=method,order_x=order,order_y=order,samplefactor=1,tol=1e-14,penalty_fn=penalty_fn,adaptive=false)

###

# plas_diff.plot_grid(gdata)
# savefig("yes2.png")

# u = soln.u
u = soln.u

println("plotting")

N = length(u)
skip = 1
fps = 10

energy = zeros(N)
maxerry = zeros(N)
maxerrx = zeros(N)
for i = 1:N
    energy[i] = norm(u[i][:,:],2)
    maxerry[i] = norm(u[i][:,1]-u[i][:,end],Inf)
    maxerrx[i] = norm(u[i][1,:]-u[i][end,:],Inf)
end


anim = @animate for i = 1:skip:N
    l = @layout [a{0.7w} [b; c]]
    p = surface(u[i][:,:],layout=l,label="t=$(@sprintf("%.5f",i*Δt))",zlims=(0.0,1.0),clims=(0.0,1.0),xlabel="y",ylabel="x",camera=(30,30))
    plot!(p[2],maxerry[1:i],ylims=(0.0,max(maximum(maxerrx),maximum(maxerry))),label="y_0 - y_N")
    plot!(p[2],maxerrx[1:i],label="x_0 - x_N")
    # plot!(p[2],u[15,:,i],ylabel="u(x=0.5)")
    plot!(p[3],energy[1:i],ylabel="||u||_2")
end
gif(anim,"yes.gif",fps=fps)


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



