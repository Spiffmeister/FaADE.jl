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




# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::SBP_operators.Helpers.Grid1D,t::Float64)
    u_MMS = zeros(grid.n)
    for i = 1:grid.n
        u_MMS[i] = MMS(grid.grid[i],t)
    end
    return u_MMS
end





function comp_MMS(Dx,npts,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        F,ũ,ũ₀,order;
        dt_scale=0.01,t_f=0.01,k=1.0)

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

    # Construct the correct problem
    function MakeProb(k)
        if (BX0Type != Periodic)
            return VariableCoefficientPDE1D(ũ₀,k,order,Bx0,BxL)
        elseif (BX0Type == Periodic)
            return VariableCoefficientPDE1D(ũ₀,k,order,Bx0L)
        end
    end

    # Loop
    for n in npts
        Dom = Grid1D(Dx,n)
        
        Δt = dt_scale*Dom.Δx^2

        K = zeros(Float64,n) .+ k

        P = MakeProb(K)

        println("Solving n=",Dom.n," case with Δt=",Δt)
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
ũ(x,t;ωx=1.0,cx=0.0) = cos(2π*t) * sin(2π*x*ωx + cx)

# Initial condition
ũ₀(x;ωx=1.0,cx=0.0) = sin(2π*ωx*x + cx)


K = 1.0
F(x,t;ωx=1.0,cx=0.0) = 
        -2π*sin(2π*t)*sin(2π*x*ωx + cx) + 
            K * 4π^2 * ωx^2 * cos(2π*t)*sin(2π*x*ωx + cx)
            
    


# Dirichlet
println("=====")
println("Dirichlet")
cx=1.0
ωx=8.5

println("ωx=",ωx,",  cx=",cx)

analytic(x,t) = ũ(x,t, ωx=ωx, cx=cx)
IC(x) = ũ₀(x, ωx=ωx, cx=cx)
FD(x,t) = F(x,t, ωx=ωx, cx=cx)

BxLũ(t) = cos(2π*t) * sin(cx) #Boundary condition x=0
BxRũ(t;Lx=1.0) = cos(2π*t) * sin(2π*Lx*ωx + cx) #Boundary condition x=Lx

order = 2
println("order=",order)
O2_DirichletMMS = comp_MMS(𝒟x,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
    FD,analytic,IC,order)

order = 4
println("order=",order)
O4_DirichletMMS = comp_MMS(𝒟x,npts,
    BxLũ,Dirichlet,BxRũ,Dirichlet,
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
ωx=8.5

println("ωx=",ωx,"  cx=",cx)

analytic(x,t) = ũ(x,t, ωx=ωx, cx=cx)
IC(x) = ũ₀(x, ωx=ωx, cx=cx)
FD(x,t) = F(x,t, ωx=ωx, cx=cx)

BxLũ(t) =         2π*ωx * cos(2π*t) * cos(cx) #Boundary condition x=0
BxRũ(t;Lx=1.0) =  2π*ωx * cos(2π*t) * cos(2π*Lx*ωx + cx) #Boundary condition x=Lx

order = 2
println("order=",order)
O2_NeumannMMS = comp_MMS(𝒟x,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    FD,analytic,IC,order)

order = 4
println("order=",order)
O4_NeumannMMS = comp_MMS(𝒟x,npts,
    BxLũ,Neumann,BxRũ,Neumann,
    FD,analytic,IC,order)

println("Order 2 Neumann convergence rates=",O2_NeumannMMS.conv_rate)
println("Order 4 Neumann convergence rates=",O4_NeumannMMS.conv_rate)

plot!(p,    O2_NeumannMMS.npts,       O2_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^2)$")
plot!(p,    O4_NeumannMMS.npts,       O4_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^4)$")

println("=====")



# Periodic
println("=====")
println("Periodic")

cx=1.0
ωx=8.0

println("ωx=",ωx,",  cx=",cx)

analytic(x,t) = ũ(x,t, ωx=ωx, cx=cx)
IC(x) = ũ₀(x, ωx=ωx, cx=cx)
FD(x,t) = F(x,t, ωx=ωx, cx=cx)

order = 2
O2_PeriodicMMS = comp_MMS(𝒟x,npts,
    nothing,Periodic,nothing,Periodic,
    FD,analytic,IC,order)

order = 4
O4_PeriodicMMS = comp_MMS(𝒟x,npts,
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
jldsave("testing/MMS/FullMMS1D.jld2";O2Conv,O4Conv)




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

=#


#=

using DelimitedFiles


open("testing/MMS/1DMMSTests_O2.csv","w") do io
    writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr O2_PeriodicMMS.relerr])
end


open("testing/MMS/1DMMSTests_O4.csv","w") do io
    writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr O4_PeriodicMMS.relerr])
end


=#