using LinearAlgebra

using Plots
using LaTeXStrings

using Pkg
Pkg.activate(".")
using SBP_operators






###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

method = :cgie
# Number of grid points in each solution
# npts = [21,31,41,51,61,71,81,91,101]

###=== MMS ===###



# Neumann boundaries
Nx0_Lũ(y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) =        2π*ωx * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy)
NxL_Rũ(y,t;
    ωx=1.0,cx=0.0,Lx=1.0,
    ωy=1.0,cy=0.0) =        2π*ωx * cos(2π*t) * cos(2π*Lx*ωx + cx)  * sin(2π*y*ωy + cy) 

Ny0_Lũ(x,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) =        2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(cy)
NyL_Rũ(x,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=1.0,Ly=1.0) = 2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(2π*Ly*ωy + cy)






# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::SBP_operators.Helpers.Grid2D,t::Float64)
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

    comp_soln = []
    MMS_soln = []
    grids = []
    relerr = []
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

    # Construct the correct problem
    function MakeProb(kx,ky)
        if (BX0Type != Periodic) & (BY0Type != Periodic)
            return VariableCoefficientPDE2D(ũ₀,kx,ky,order,Bx0,BxL,By0,ByL)
        elseif (BX0Type != Periodic) & (BY0Type == Periodic) 
            return VariableCoefficientPDE2D(ũ₀,kx,ky,order,Bx0,BxL,By0L)
        elseif (BX0Type == Periodic) & (BY0Type != Periodic)
            return VariableCoefficientPDE2D(ũ₀,kx,ky,order,Bx0L,By0,ByL)
        else
            return VariableCoefficientPDE2D(ũ₀,kx,ky,order,Bx0L,By0L)
        end
    end

    # Loop
    for n in npts
        Dom = Grid2D(Dx,Dy,n,n)
        
        Δt = dt_scale*Dom.Δx^2

        Kx = zeros(Float64,n,n) .+ kx
        Ky = zeros(Float64,n,n) .+ ky

        P = MakeProb(Kx,Ky)

        println("Solving n=",Dom.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,:cgie,source=F)

        u_MMS = generate_MMS(ũ,Dom,t_f)

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        push!(relerr, norm(u_MMS .- soln.u[2])/norm(MMS_soln))
    end

    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end




###=== MMS TESTS ===###
p = plot()

npts = [21,31,41,51,61,71,81,91,101]


# Solution
ũ(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)

# Initial condition
ũ₀(x,y;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = sin(2π*ωx*x + cx) * sin(2π*ωy*y + cy)


K = 1.0
F(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = 
        -2π*sin(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
            K * 4π^2 * (ωx^2 + ωy^2) * cos(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) 
            
    


# Dirichlet
println("=====")
println("Dirichlet")
cx=1.0
cy=0.0
ωx=9.0
ωy=7.5

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)

BxLũ(y,t) = cos(2π*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x=0
BxRũ(y,t;Lx=1.0) = cos(2π*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=Lx
ByLũ(x,t) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(cy) #Boundary condition y=0
ByRũ(x,t;Ly=1.0) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*Ly*ωy + cy) #Boundary condition y=Ly

order = 2
println("order=",order)
O2_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    ByLũ,Dirichlet,ByRũ,Dirichlet,
    FD,analytic,IC,order)

order = 4
println("order=",order)
O4_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    ByLũ,Dirichlet,ByRũ,Dirichlet,
    FD,analytic,IC,order)

println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

plot!(p,    O2_DirichletMMS.npts,     O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$")
plot!(p,    O4_DirichletMMS.npts,     O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$")

println("=====")



# Neumann
println("=====")
println("Neumann")

cx=1.0
cy=0.0
ωx=9.0
ωy=7.5

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)

BxLũ(y,t) =         2π*ωx * cos(2π*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0
BxRũ(y,t;Lx=1.0) =  2π*ωx * cos(2π*t) * cos(2π*Lx*ωx + cx)  * sin(2π*y*ωy + cy) #Boundary condition x=Lx
ByLũ(x,t) =         2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(cy) #Boundary condition y=0
ByRũ(x,t;Ly=1.0) =  2π*ωy * cos(2π*t) * sin(2π*x*ωx + cx)   * cos(2π*Ly*ωy + cy) #Boundary condition y=Ly

order = 2
println("order=",order)
O2_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    ByLũ,Neumann,ByRũ,Neumann,
    FD,analytic,IC,order)

order = 4
println("order=",order)
O4_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    ByLũ,Neumann,ByRũ,Neumann,
    FD,analytic,IC,order)

println("Order 2 Neumann convergence rates=",O2_NeumannMMS.conv_rate)
println("Order 4 Neumann convergence rates=",O4_NeumannMMS.conv_rate)

plot!(p,    O2_NeumannMMS.npts,       O2_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^2)$")
plot!(p,    O4_NeumannMMS.npts,       O4_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^4)$")

println("=====")



# Periodic
println("=====")
println("Dirichlet/Periodic")

cx=1.0
cy=0.0
ωx=7.5
ωy=6.0

println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

analytic(x,y,t) = ũ(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
FD(x,y,t) = F(x,y,t, ωx=ωx, cx=cx, ωy=ωy, cy=cy)

BxLũ(y,t) = cos(2π*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x=0
BxRũ(y,t;Lx=1.0) = cos(2π*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=Lx

order = 2
O2_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order)

order = 4
O4_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order)

println("Order 2 Dirichlet/Periodic convergence rates=",O2_PeriodicMMS.conv_rate)
println("Order 4 Dirichlet/Periodic convergence rates=",O4_PeriodicMMS.conv_rate)

plot!(p,    O2_PeriodicMMS.npts,      O2_PeriodicMMS.relerr,      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$")
plot!(p,    O4_PeriodicMMS.npts,      O4_PeriodicMMS.relerr,      label=L"Dirichlet/Periodic $\mathcal{O}(h^4)$")

println("=====")


plot!(p,xaxis=:log,yaxis=:log)


order2rate = npts


#=


savefig(p,".//testing//MMS//MMSTests.eps")
=#



#=

pO2 = plot()

plot!(pO2,    log.(O2_DirichletMMS.npts),     log.(O2_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
plot!(pO2,    log.(O2_NeumannMMS.npts),       log.(O2_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^2)$")
plot!(pO2,    log.(O2_PeriodicMMS.npts),      log.(O2_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$")

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

=#