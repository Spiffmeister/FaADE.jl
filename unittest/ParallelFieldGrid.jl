using LinearAlgebra

cd("..")
using Interpolations
push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./SBP_operators")
using SBP_operators
using plas_diff




###
𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 21
ny = 21
Dom = Grid2D(𝒟x,𝒟y,nx,ny)



# params = plas_diff.SampleFields.H_params([0.],[0.],[0.])
χₘₙ = 2.1e-3 + 5.0e-3
params = plas_diff.SampleFields.H_params([χₘₙ/2., χₘₙ/3.],[2.0, 3.0],[1.0, 2.0])

function χ_h!(χ,x::Array{Float64},p,t)
    # Hamiltons equations for the field-line Hamiltonian
    # H = ψ²/2 - ∑ₘₙ ϵₘₙ(cos(mθ - nζ))
    χ[1] = x[2] #p_1            qdot        θ
    χ[2] = -sum(p.ϵₘₙ .*(sin.(p.m*x[1] - p.n*t) .* p.m)) #q_1        pdot        ψ
end



###
𝒟x = [0.0,1.0]
𝒟y = [-π,π]
nx = 6
ny = 6
Dom = Grid2D(𝒟x,𝒟y,nx,ny)

gdata = plas_diff.construct_grid(𝒟x,𝒟y,nx,ny,χ_h!,params)

pdata = plas_diff.poincare(χ_h!,params,x=[0.0,1.0],y=[-π,π])





# using GLMakie
using CairoMakie




# p3 = Figure(resolution=(1600,1200),fontsize=20)
p3 = Figure()
ax3 = Axis3(p3[1,1],xlabel="",ylabel="y",zlabel="x",
    ylabelsize=30,
    zlabelsize=30)

ψ = repeat(gdata.x,1,gdata.ny);
θ = repeat(gdata.y',gdata.nx,1);

scatter!(ax3,zeros(size(ψ))[:],θ[:],ψ[:],color=:black,label="ζ=0 plane")
wireframe!(ax3,zeros(size(ψ)),θ,ψ,color=:grey)

pickapoint = CartesianIndex((4,4))


f(x,y,t) = Point3f(
    π*sin.(t),
    -(π*cos.(t) .+ (-(θ[pickapoint]*(1.0.-t/2π) + y*t/2π) .- π)),
    (ψ[pickapoint]*(2π.-t) + x*t)/2π
    )

#=
F(x) = Point3f(
    x[3],
    x[2],
    -sum(χₘₙ*[1/2, 1/3].*sin.([2,3]*x[1] - [1,2]*x[3]).*[2,3]))
    
streamplot(F,0.0..1.0,-π..π,0.0..2π,arrowsize=1e-5,density=0.1,arrow_size=0.1)
=#

t = collect(range(0.0,2π,length=100));
zplane = gdata.z_planes[1]

pts = f.(zplane.x[pickapoint],zplane.y[pickapoint],t)
ax3_forwardln = lines!(ax3, pts)
ax3_forwardpt = scatter!(ax3,[pts[end]],markercolor=ax3_forwardln.color)
ax3_forwardarrow = arrows!(ax3,[pts[50]],[Point3f([-1e-2,0.0,0.0])],linecolor=ax3_forwardln.color,arrowcolor=ax3_forwardln.color,arrowsize=[1e-1,1e-1,1e-2])

# Backward

zplane = gdata.z_planes[2]

f(x,y,t) = Point3f(
    -π*sin.(t),
    -(π*cos.(t) .+ (-(θ[pickapoint]*(1.0.-t/2π) + y*t/2π) .- π)),
    (ψ[pickapoint]*(2π.-t) + x*t)/2π
    )
    
t = collect(range(0.0,2π,length=100));
pts = f.(zplane.x[pickapoint],zplane.y[pickapoint],t)

ax3_backwardln = lines!(ax3, pts)
ax3_backwardpt = scatter!(ax3,[pts[end]],markercolor=ax3_backwardln.color)
ax3_backwardarrow = arrows!(ax3,[pts[50]],[Point3f([1e-2,0.0,0.0])],linecolor=ax3_backwardln.color,arrowcolor=ax3_backwardln.color,arrowsize=[1e-1,1e-1,1e-2])

# Pretty bit

ax3.yticks = (-π:π/2:π,["π", "-π/2", "0", "π/2", "π"])
ax3.xticks = ([0],["0"])

axislegend(ax3,[[ax3_forwardpt,ax3_forwardln],[ax3_backwardpt,ax3_backwardln]],["Foward","Backward"],"Field line tracing direction",position=:rt,labelsize=30)



save("FieldLineTracing.pdf", p3)#, resolution=(1600,1200), transparency=true)



#=
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
