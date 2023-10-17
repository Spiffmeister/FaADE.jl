using LinearAlgebra
using Printf
using Plots
using JLD2

using BenchmarkTools
using Profile
# using PProf

using Interpolations
# push!(LOAD_PATH,"..")
using FaADE




###
𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 81
ny = 81
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

kx(x,y) = 1.0
ky(x,y) = 1.0


Δt = 1.0 * min(Dom.Δx^2,Dom.Δy^2)
# t_f = 100Δt
t_f = 10.0

u₀(x,y) = x


BoundaryLeft = Boundary(Dirichlet,(y,t) -> 0.0,Left,1)
BoundaryRight = Boundary(Dirichlet,(y,t) -> 1.0,Right,1)
BoundaryUpDown = PeriodicBoundary(2)


order = 2
method = :cgie

P = VariableCoefficientPDE2D(u₀,ky,kx,order,BoundaryLeft,BoundaryRight,BoundaryUpDown)



# params = plas_diff.SampleFields.H_params([0.],[0.],[0.])
χₘₙ = 2.1e-3 + 5.0e-3
params = (ϵₘₙ = [χₘₙ/2., χₘₙ/3.],m=[2.0, 3.0],n=[1.0, 2.0])
function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[2] = x[1] #p_1            qdot        θ
    χ[1] = -sum(p.ϵₘₙ .*(sin.(p.m*x[2] - p.n*t) .* p.m)) #q_1        pdot        ψ
end

dH(X,x,p,t) = χ_h!(X,x,params,t)
PGrid = FaADE.construct_grid(dH,Dom,[-2π,2π])
# Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2,κ=1e8)
PData = FaADE.ParallelData(PGrid,Dom,κ=1e12)


println("(Δx,Δy)=(",Dom.Δx,",",Dom.Δy,")      ","Δt=",Δt,"        ","final time=",t_f)


# using Profile

Pfn = FaADE.generate_parallel_penalty(PGrid,Dom,2)
# @time soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=true,penalty_func=Pfn)
@time soln = solve(P,Dom,Δt,t_f,:cgie,adaptive=true,Pgrid=PData)

surface(soln.u[2])

include("../../paper_JCP2023/FieldLines.jl")

poindata = FieldLines.construct_poincare(dH,[0.0,1.0],[-π,π])
scatter(poindata.θ,poindata.ψ,markercolor=:black,markersize=0.7,ylims=𝒟x,xlims=𝒟y)
contour!(Dom.gridy,Dom.gridx,soln.u[2],linewidth=3)

#=
println("Plotting")
using Plots
surface(soln.grid.gridy,soln.grid.gridx,soln.u[2],
    xlabel="y",ylabel="x",zlabel="Temp")

# @time solve(P,Dom,Δt,t_f,:cgie)
# Profile.clear_malloc_data()
# @time solve(P,Dom,Δt,t_f,:cgie)



pdata = plas_diff.poincare(χ_h!,params,x=[0.0,1.0],y=[-π,π])


plas_diff.plot_grid(gdata)

=#


#=
p1 = scatter(pdata.θ,pdata.ψ,markercolor=:black,markersize=0.7,ylims=𝒟x,xlims=𝒟y,ylabel="ψ",xlabel="θ",legend=false,dpi=600,fmt=:png)
savefig(p1,"FaADE//figures//CTAC_Poincare")


p2 = contour(soln.grid.gridy,soln.grid.gridx,soln.u[2],dpi=600,fmt=:png,linewidth=2)
scatter!(pdata.θ,pdata.ψ,markercolor=:black,markersize=0.7,ylims=𝒟x,ylabel="ψ",xlabel="θ",legend=false)
savefig(p2,"FaADE//figures/CTAC_Contour")




###
𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 6
ny = 6
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)

pdata = plas_diff.poincare(χ_h!,params,x=[0.0,1.0],y=[-π,π])




p3 = GLMakie.Figure()
ax3 = GLMakie.Axis3(p3[1,1])

ψ = repeat(gdata.x,1,gdata.ny)
θ = repeat(gdata.y',gdata.nx,1)

GLMakie.scatter!(ax3,zeros(size(ψ))[:],θ[:],ψ[:],color=:black,label="ζ=0 plane")
GLMakie.wireframe!(ax3,zeros(size(ψ)),θ,ψ,color=:grey)

t = collect(range(0.0,2π,length=100));

pickapoint = CartesianIndex((4,4))
GLMakie.scatter!(ax3,zeros(2), [θ[pickapoint],gdata.z_planes[1].y[pickapoint]], [ψ[pickapoint],gdata.z_planes[1].x[pickapoint]])
GLMakie.lines!(ax3, π*sin.(t), -(π*cos.(t) .+ (-(θ[pickapoint]*(1.0.-t/2π) + gdata.z_planes[1].y[pickapoint]*t/2π) .- π)), collect(range(ψ[pickapoint],gdata.z_planes[1].x[pickapoint],length=100)) )


GLMakie.save("FaADE//figures//AustMS_2022.png", p3, resolution=(1600,1200), transparency=true)





p3 = GLMakie.Figure()
ax3 = GLMakie.Axis3(p3[1,1])

ψ = repeat(gdata.x,1,gdata.ny)
θ = repeat(gdata.y',gdata.nx,1)

obsZ = Observable([zeros(size(ψ))[:] θ[:] ψ[:]])
GLMakie.scatter!(ax3,obsZ,color=:black,label="FD Mesh")


flt = Observable(Point3f[])
flt_ends = Observable(Point3f[])

GLMakie.scatter!(ax3,flt_ends,marker=:cross)

GLMakie.Legend()

on(events(ax3).mousebutton, priority=2) do event
    if event.button == Mouse.left
        if event.action == Mouse.press
            plt, pickapoint = pick(ax3)
            if plt != nothing
                GLMakie.lines!(ax3, π*sin.(t), -(π*cos.(t) .+ (-(θ[pickapoint]*(1.0.-t/2π) + gdata.z_planes[1].y[pickapoint]*t/2π) .- π)), collect(range(ψ[pickapoint],gdata.z_planes[1].x[pickapoint],length=100)) )

                push!(flt_ends[], Point3f(0.0, gdata.z_planes[1].y[pickapoint], gdata.z_planes[1].x[pickapoint]))
                notify(flt_ends)

                return Consume(true)
            end
        end
    # elseif event.button == Mouse.right
    #     if event.action == Mouse.press
    #         plt, pickapoint = pick(ax3)
    #         if plt != nothing
    #             deleteat!(flt_ends[],pickapoint)
    #         end
    #     end
    end
end






GLMakie.scatter(ax3,soln.grid.gridy,soln.grid.gridx)



GLMakie.surface!(ax3,soln.grid.gridy,soln.grid.gridx,soln.u[2]',colormap=(:viridis, 0.5),transparency=true,alpha=0.5)
GLMakie.wireframe!(ax3,soln.grid.gridy,soln.grid.gridx,soln.u[2]',color=(:black,0.2),transparency=true,linewidth=1.0)
GLMakie.scatter!(ax3,pdata.θ[0.0 .≤ pdata.ψ .≤ 1.0],pdata.ψ[0.0 .≤ pdata.ψ .≤ 1.0],zeros(length(pdata.ψ[0.0 .≤ pdata.ψ .≤ 1.0])),color=:black,markersize=2.0)


GLMakie.scale!(ax3,(1.0,2.0,1.0))


using GLMakie
# GLMakie.wireframe(soln.grid.gridy,soln.grid.gridx,soln.u[2])
# GLMakie.scatter(pdata.θ[0.0 .≤ pdata.ψ .≤ 1.0],pdata.ψ[0.0 .≤ pdata.ψ .≤ 1.0],markersize=1.0,color=:black)
# GLMakie.contour!(soln.grid.gridy,soln.grid.gridx,soln.u[2]',linewidth=2,levels=10)
# GLMakie.Colorbar!()


pdata = plas_diff.poincare(χ_h!,params)


scatter(pdata.θ,pdata.ψ,ylims=(0.0,1.0),xlims=(0,2π))

plas_diff.plot_grid(gdata)

scatter!(gdata.z_planes[1].y[:],gdata.z_planes[1].x[:])
scatter!(gdata.z_planes[2].y[:],gdata.z_planes[2].x[:])

=#
