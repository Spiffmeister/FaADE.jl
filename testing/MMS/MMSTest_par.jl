using Distributed
addprocs(4)

@everywhere using LinearAlgebra

# using Plots
# using LaTeXStrings

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using SPADE






###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

###=== MMS ===###







# Generates the exact MMS solution
@everywhere function generate_MMS(MMS::Function,grid::SPADE.Helpers.Grid2D,t::Float64)
    u_MMS = zeros(grid.nx,grid.ny)
    for j = 1:grid.ny
        for i = 1:grid.nx
            u_MMS[i,j] = MMS(grid.gridx[i],grid.gridy[j],t)
        end
    end
    return u_MMS
end





function comp_MMS(Dx,Dy,npts,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        BoundaryY0,BY0Type,BoundaryYL,BYLType,
        F,ũ,ũ₀,order;
        dt_scale=0.01,t_f=0.01,kx=1.0,ky=kx)

    comp_soln   = []
    MMS_soln    = []
    grids       = []
    relerr      = []

    # Loop
    @everywhere function subtest(n,Dx,Dy,kx,ky,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        BoundaryY0,BY0Type,BoundaryYL,BYLType,
        F,ũ,ũ₀,order,
        dt_scale,t_f)

        Kx = zeros(Float64,n,n) .+ kx
        Ky = zeros(Float64,n,n) .+ ky

        # X boundaries
        if BX0Type != Periodic
            Bx0 = Boundary(BX0Type,BoundaryX0,Left,1)
            BxL = Boundary(BXLType,BoundaryXL,Right,1)
        else
            Bx0L = PeriodicBoundary(1)
        end
        # Y boundaries
        if BY0Type != Periodic
            By0 = Boundary(BY0Type,BoundaryY0,Up,2)
            ByL = Boundary(BYLType,BoundaryYL,Down,2)
        else
            By0L = PeriodicBoundary(2)
        end


        if (BX0Type != Periodic) & (BY0Type != Periodic)
            P = VariableCoefficientPDE2D(ũ₀,Kx,Ky,order,Bx0,BxL,By0,ByL)
        elseif (BX0Type != Periodic) & (BY0Type == Periodic) 
            P = VariableCoefficientPDE2D(ũ₀,Kx,Ky,order,Bx0,BxL,By0L)
        elseif (BX0Type == Periodic) & (BY0Type != Periodic)
            P = VariableCoefficientPDE2D(ũ₀,Kx,Ky,order,Bx0L,By0,ByL)
        else
            P = VariableCoefficientPDE2D(ũ₀,Kx,Ky,order,Bx0L,By0L)
        end

        Dom = Grid2D(Dx,Dy,n,n)
        
        Δt = dt_scale*Dom.Δx^2

        println("Solving n=",Dom.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,:cgie,source=F)

        # u_MMS = generate_MMS(ũ,Dom,t_f)
        u_MMS = zeros(Dom.nx,Dom.ny)
        for j = 1:Dom.ny
            for i = 1:Dom.nx
                u_MMS[i,j] = ũ(Dom.gridx[i],Dom.gridy[j],t_f)
            end
        end

        test = (soln=soln,grid=Dom,MMS=u_MMS)
        return test
    end


    test = pmap(n -> subtest(n,
    Dx,Dy,kx,ky,BoundaryX0,BX0Type,BoundaryXL,BXLType,BoundaryY0,BY0Type,BoundaryYL,BYLType,F,ũ,ũ₀,order,dt_scale,t_f),
        npts)

    for i in 1:length(npts)
        push!(comp_soln,test[i].soln)
        push!(grids,test[i].grid)
        push!(MMS_soln,test[i].MMS)

        push!(relerr, norm(test[i].MMS .- test[i].soln.u[2])/norm(test[i].MMS))
    end
    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
    # return test
end




###=== MMS TESTS ===###

npts = [21,31,41,51,61,71,81,91,101]


# Solution
@everywhere ũ(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)

# Initial condition
@everywhere ũ₀(x,y;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = sin(2π*ωx*x + cx) * sin(2π*ωy*y + cy)


@everywhere K = 1.0

@everywhere F(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0,
    K = 1.0) = 
        -2π*sin(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
            K * 4π^2 * (ωx^2 + ωy^2) * cos(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) 
            
    
println("=== K=",K," ===")

# Dirichlet
println("=====")
println("Dirichlet")
@everywhere cx=1.0
@everywhere cy=0.0
@everywhere ωx=9.0
@everywhere ωy=7.5

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

@everywhere analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

@everywhere BxLũ(y,t)           = cos(2π*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x=0
@everywhere BxRũ(y,t;Lx=1.0)    = cos(2π*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=Lx
@everywhere ByLũ(x,t)           = cos(2π*t) * sin(2π*x*ωx + cx) * sin(cy) #Boundary condition y=0
@everywhere ByRũ(x,t;Ly=1.0)    = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*Ly*ωy + cy) #Boundary condition y=Ly

order = 2
println("order=",order)
O2_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    ByLũ,Dirichlet,ByRũ,Dirichlet,
    FD,analytic,IC,order,
    kx=K,ky=K)

order = 4
println("order=",order)
O4_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    ByLũ,Dirichlet,ByRũ,Dirichlet,
    FD,analytic,IC,order,
    kx=K,ky=K)

println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

println("=====")



# Neumann
println("=====")
println("Neumann")

cx=1.0
cy=0.0
ωx=9.0
ωy=7.5

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

@everywhere analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K=K)

@everywhere BxLũ(y,t) =         2π*ωx * K * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0
@everywhere BxRũ(y,t;Lx=1.0) =  2π*ωx * K * cos(2π*t) * cos(2π*Lx*ωx + cx)  * sin(2π*y*ωy + cy) #Boundary condition x=Lx
@everywhere ByLũ(x,t) =         2π*ωy * K * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(cy) #Boundary condition y=0
@everywhere ByRũ(x,t;Ly=1.0) =  2π*ωy * K * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(2π*Ly*ωy + cy) #Boundary condition y=Ly

order = 2
println("order=",order)
O2_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    ByLũ,Neumann,ByRũ,Neumann,
    FD,analytic,IC,order,
    kx=K, ky=K)

order = 4
println("order=",order)
O4_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    ByLũ,Neumann,ByRũ,Neumann,
    FD,analytic,IC,order,
    kx=K, ky=K)

println("Order 2 Neumann convergence rates=",O2_NeumannMMS.conv_rate)
println("Order 4 Neumann convergence rates=",O4_NeumannMMS.conv_rate)

println("=====")



# Periodic
println("=====")
println("Dirichlet/Periodic")

cx=1.0
cy=0.0
ωx=2.0
ωy=2.0

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

@everywhere analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
@everywhere FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K=K)

@everywhere BxLũ(y,t)           = cos(2π*t) * sin(cx)               * sin(2π*y*ωy + cy) #Boundary condition x=0
@everywhere BxRũ(y,t;Lx=1.0)    = cos(2π*t) * sin(2π*Lx*ωx + cx)    * sin(2π*y*ωy + cy) #Boundary condition x=Lx

order = 2
O2_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order,
    kx=K, ky=K)

order = 4
O4_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order,
    kx=K, ky=K)

println("Order 2 Dirichlet/Periodic convergence rates=",O2_PeriodicMMS.conv_rate)
println("Order 4 Dirichlet/Periodic convergence rates=",O4_PeriodicMMS.conv_rate)

println("=====")






O2Conv = (n=npts,
    conv_D = O2_DirichletMMS.conv_rate,
    conv_N = O2_NeumannMMS.conv_rate,
    conv_P = O2_PeriodicMMS.conv_rate,
    relerr_D = O2_DirichletMMS.relerr,
    relerr_N = O2_NeumannMMS.relerr,
    relerr_P = O2_PeriodicMMS.relerr
    )

O4Conv = (n=npts,
    conv_D = O4_DirichletMMS.conv_rate,
    conv_N = O4_NeumannMMS.conv_rate,
    conv_P = O4_PeriodicMMS.conv_rate,
    relerr_D = O4_DirichletMMS.relerr,
    relerr_N = O4_NeumannMMS.relerr,
    relerr_P = O4_PeriodicMMS.relerr
    )

using JLD2
jldsave("testing/MMS/FullMMS.jld2";O2Conv,O4Conv)




#=

pO2 = plot(axis=:log,minorgrid=true)
plot!(pO2,    (O2_DirichletMMS.npts),     (O2_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
plot!(pO2,    (O2_NeumannMMS.npts),       (O2_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^2)$", markershape=:square)
plot!(pO2,    (O2_PeriodicMMS.npts),      (O2_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$", markershape=:x)

plot!(pO2,    log.(O2_DirichletMMS.npts),     log.(O2_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
plot!(pO2,    log.(O2_NeumannMMS.npts),       log.(O2_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^2)$"), markershape=:circle
plot!(pO2,    log.(O2_PeriodicMMS.npts),      log.(O2_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$", markershape=:circle)

plot!(pO2, log.([npts[2],npts[end-1]]), -log.([npts[2],npts[end-1]].^2) .+ log(npts[2]^2) .+ log(sum(O2_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dash, linecolor=:black,label=L"$\mathcal{O}(h^2)$")

savefig(pO2,".//testing//MMS//MMSTests_order2.eps")
savefig(pO2,".//testing//MMS//MMSTests_order2.png")


pO4 = plot()
plot!(pO4,    log.(O4_DirichletMMS.npts),     log.(O4_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    log.(O4_NeumannMMS.npts),       log.(O4_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    log.(O4_PeriodicMMS.npts),      log.(O4_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^4)$", markershape=:x)

plot!(pO4, log.([npts[2],npts[end-1]]), 
    -log.([npts[2],npts[end-1]].^2) .+ log(npts[2]^2) .+ log(sum(O4_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dash, linecolor=:black, label=L"$\mathcal{O}(h^2)$")
plot!(pO4, log.([npts[2],npts[end-1]]), 
    -log.([npts[2],npts[end-1]].^4) .+ log(npts[2]^4) .+ log(sum(O4_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dashdot, linecolor=:black, label=L"$\mathcal{O}(h^4)$")

savefig(pO4,".//testing//MMS//MMSTests_order4.eps")
savefig(pO4,".//testing//MMS//MMSTests_order4.png")

surface(O4_PeriodicMMS.comp_soln[end].u[2] .- O4_PeriodicMMS.MMS_soln[end])


surface(O4_PeriodicMMS.grids[end].gridx,O4_PeriodicMMS.grids[end].gridy,O4_PeriodicMMS.comp_soln[end].u[2] .- O4_PeriodicMMS.MMS_soln[end],xlabel="x",ylabel="y")
surface(O4_PeriodicMMS.grids[end].gridx,O4_PeriodicMMS.grids[end].gridy,O4_PeriodicMMS.comp_soln[end].u[2] .- O4_PeriodicMMS.MMS_soln[end],xlabel="x",ylabel="y")

surface(O2_PeriodicMMS.grids[4].gridx,O2_PeriodicMMS.grids[4].gridy,O2_PeriodicMMS.comp_soln[4].u[2] .- O2_PeriodicMMS.MMS_soln[4],xlabel="x",ylabel="y")
surface(O2_PeriodicMMS.grids[end].gridx,O2_PeriodicMMS.grids[end].gridy,O2_PeriodicMMS.comp_soln[end].u[2] .- O2_PeriodicMMS.MMS_soln[end],xlabel="x",ylabel="y")



surface(O4_DirichletMMS.grids[end].gridx,O4_DirichletMMS.grids[end].gridy,O4_DirichletMMS.comp_soln[end].u[2] .- O4_DirichletMMS.MMS_soln[end],xlabel="x",ylabel="y")
surface(O2_DirichletMMS.grids[end].gridx,O2_DirichletMMS.grids[end].gridy,O2_DirichletMMS.comp_soln[end].u[2] .- O2_DirichletMMS.MMS_soln[end],xlabel="x",ylabel="y")

=#




using DelimitedFiles

nameappend=string("K=",K)

open(string("testing/MMS/MMS_Tests_O2",nameappend,".csv"),"w") do io
    writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr O2_PeriodicMMS.relerr])
end
open(string("testing/MMS/MMS_Rates_O2",nameappend,".csv"),"w") do io
    writedlm(io,[O2_DirichletMMS.conv_rate O2_NeumannMMS.conv_rate O2_PeriodicMMS.conv_rate])
end

open(string("testing/MMS/MMS_Tests_O4",nameappend,".csv"),"w") do io
    writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr O4_PeriodicMMS.relerr])
end
open(string("testing/MMS/MMS_Rates_O4",nameappend,".csv"),"w") do io
    writedlm(io,[O4_DirichletMMS.conv_rate O4_NeumannMMS.conv_rate O4_PeriodicMMS.conv_rate])
end


