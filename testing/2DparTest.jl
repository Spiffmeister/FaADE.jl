using LinearAlgebra
using Printf
using Plots
using JLD2

using BenchmarkTools
using Profile
# using PProf

cd("..")
using Interpolations
push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./SBP_operators")
using SBP_operators
using plas_diff




###
𝒟x = [0.5,0.68]
𝒟y = [-π,π]
nx = 41
ny = 41
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx(x,y) = 1.0e-8;
ky(x,y) = 1.0e-8;


Δt = 1.0 * min(Dom.Δx^2,Dom.Δy^2)
t_f = 1000Δt

u₀(x,y) = (x-0.5)/(0.68-0.5)


BoundaryLeft = Boundary(Dirichlet,(y,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(y,t) -> 1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)


order = 2
method = :cgie


println(method)
println("(Δx,Δy)=(",Dom.Δx,",",Dom.Δy,")      ","Δt=",Δt,"        ","final time=",t_f)

P = VariableCoefficientPDE2D(u₀,kx,ky,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)



# params = plas_diff.SampleFields.H_params([0.],[0.],[0.])
χₘₙ = 2.1e-3
params = plas_diff.SampleFields.H_params([χₘₙ/2., χₘₙ/3.],[2.0, 3.0],[1.0, 2.0])

function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[1] = x[2] #p_1            qdot        θ
    χ[2] = -sum(p.ϵₘₙ .*(sin.(p.m*x[1] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)

H_x = SBP_operators.build_H(order,ny)
H_x = 1.0 ./H_x.^2

H_y = SBP_operators.build_H(order,nx)
H_y = 1.0 ./H_y.^2

# Hinv = diag(kron(I(nx),H_y) + kron(H_x,I(ny)))

κ_para = 1.0
τ_para = -1.0




### Parallel Penalty ###
function penalty_fn(u,uₒ,Δt)
    # umw = zeros(Float64,nx,ny)

    interp = LinearInterpolation((Dom.gridx,Dom.gridy),uₒ)

    for j = 1:ny
        for i = 1:nx

            # umw[i,j] = 2u[i,j] - (uₒ[gdata.z_planes[1].xproj[i,j],gdata.z_planes[1].yproj[i,j]] + uₒ[gdata.z_planes[2].xproj[i,j],gdata.z_planes[2].yproj[i,j]])


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

            u[i,j] = 1.0/(1.0 - κ_para * τ_para/2.0 * Δt * (H_y[i] + H_x[j])) *
                (uₒ[i,j] - Δt*κ_para*τ_para/4.0 *(H_y[i] + H_x[j])*(w_f + w_b))

        end
    end
    # return u#, norm(umw)
end



# println("Benchmarking")
# @benchmark solve($P,$Dom,$Δt,$t_f,:cgie,penalty_func=$penalty_fn)

@time soln = solve(P,Dom,Δt,5.1Δt,:cgie,adaptive=true,penalty_func=penalty_fn)
@time soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=true,penalty_func=penalty_fn)
println("Plotting")
using Plots
surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
    xlabel="y",ylabel="x",zlabel="Temp")

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)

