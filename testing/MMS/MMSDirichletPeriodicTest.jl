using LinearAlgebra

using Plots

using FaADE






###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

###=== MMS ===###

# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::Grid2D,t::Float64)
    u_MMS = zeros(grid.nx,grid.ny)
    for j = 1:grid.ny
        for i = 1:grid.nx
            u_MMS[i,j] = MMS(grid.gridx[i,j],grid.gridy[i,j],t)
        end
    end
    return u_MMS
end





function comp_MMS(Dx,Dy,npts,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        BoundaryY0,BY0Type,BoundaryYL,BYLType,
        F,ũ,ũ₀,order;
        dt_scale=0.1,t_f=0.1,kx=1.0,ky=kx,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    relerr = []

    # Loop
    for n in npts
        Dom = Grid2D(Dx,Dy,n,n)

        # X boundaries
        Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left,  1,order)
        BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right, 1,order)
        By0 = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Up)
        ByL = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Down)

        BD = FaADE.Inputs.SATBoundaries(Bx0,BxL,By0,ByL)


        Δt = dt_scale*Dom.Δx
        nt = round(t_f/Δt)
        Δt = t_f/nt

        P = newProblem2D(order,ũ₀,kx,ky,Dom,BD,F,nothing)

        println("Solving n=",Dom.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=θ)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2])

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        push!(relerr, norm(u_MMS .- soln.u[2])/norm(u_MMS))
    end

    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end




###=== MMS TESTS ===###
npts = [21,31,41,51,61,71,81,91,101]

nameappend = "spatial"
θ = 0.5
cx=1.0
cy=0.0
ωx=7.5
ωy=5.0
ωt=1.0

# Solution
ũ(x,y,t;
    ωt=1.0,
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)

# Initial condition
ũ₀(x,y;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = sin(2π*ωx*x + cx) * sin(2π*ωy*y + cy)


K = 1.0

F(x,y,t;
    ωt=1.0,
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0,
    K = 1.0) = 
        -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
            K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
            K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
            
    
println("=== K=",K," ===")

println("=====")
println("Dirichlet")


println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
FD(x,y,t) = F(x,y,t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

BxLũ(y,t)           = cos(2π*ωt*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x=0
BxRũ(y,t;Lx=1.0)    = cos(2π*ωt*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=Lx

order = 2
println("order=",order)
O2_DirichletPeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order,
    kx=K,ky=K,θ=θ)

order = 4
println("order=",order)
O4_DirichletPeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order,
    kx=K,ky=K,θ=θ)

println("Order 2 Dirichlet convergence rates=",O2_DirichletPeriodicMMS.conv_rate)
println("Order 4 Dirichlet convergence rates=",O4_DirichletPeriodicMMS.conv_rate)

# pD = plot(axis=:log,minorgrid=true)
# plot!(pD,  O2_DirichletPeriodicMMS.npts,   O2_DirichletPeriodicMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
# plot!(pD,  O2_DirichletPeriodicMMS.npts,   O2_DirichletPeriodicMMS.npts.^2,     label=L"$\mathcal{O}(h^2)$", markershape=:circle)
# plot!(pD,  O4_DirichletPeriodicMMS.npts,   O4_DirichletPeriodicMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:circle)
# plot!(pD,  O4_DirichletPeriodicMMS.npts,   O4_DirichletPeriodicMMS.npts.^4,     label=L"$\mathcal{O}(h^4)$", markershape=:circle)    
# savefig(pD,"2DMMSDirichlet.png")


println("=====")





O2Conv = (n=npts,
    conv_DP = O2_DirichletPeriodicMMS.conv_rate,
    relerr_DP = O2_DirichletPeriodicMMS.relerr,
    )

O4Conv = (n=npts,
    conv_DP = O4_DirichletPeriodicMMS.conv_rate,
    relerr_DP = O4_DirichletPeriodicMMS.relerr,
    )

using JLD2
jldsave("testing/MMS/DirichletPeriodicMMS.jld2";O2Conv,O4Conv)


using DelimitedFiles

open(string("testing/MMS/MMS_DP_",nameappend,"_Tests_theta",θ,".csv"),"w") do io
    writedlm(io,[npts O2_DirichletPeriodicMMS.relerr O4_DirichletPeriodicMMS.relerr])
end
open(string("testing/MMS/MMS_DP_",nameappend,"_Rates_theta",θ,".csv"),"w") do io
    writedlm(io,[O2_DirichletPeriodicMMS.conv_rate O4_DirichletPeriodicMMS.conv_rate])
end



