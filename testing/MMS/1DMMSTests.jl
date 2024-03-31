using LinearAlgebra

# using Plots
# using LaTeXStrings

using FaADE



rundirichlet = true
runneumann = true
runperiodic = true

###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]


###=== MMS ===###

# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::Grid1D,t::Float64)
    u_MMS = zeros(grid.n)
    for i = 1:grid.n
        u_MMS[i] = MMS(grid.grid[i],t)
    end
    return u_MMS
end





function comp_MMS(Dx,npts,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        F,ũ,ũ₀,order;
        dt_scale=0.10,t_f=5.0,k=1.0,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    relerr = []
    interr = []



    # Loop
    for n in npts
        Dom = Grid1D(Dx,n)
        
        # X boundaries
        if BX0Type == Periodic
            Bx0 = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
            BxL = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
        elseif BX0Type == Dirichlet
            Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left, order)
            BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right,order)
        elseif BX0Type == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,1,order)
        end
        BD = FaADE.Inputs.SATBoundaries(Bx0,BxL)

        Δt = dt_scale*Dom.Δx

        K = k

        P = Problem1D(order,ũ₀,K,Dom,BD,F,nothing)

        println("Solving n=",Dom.n," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=θ)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2])

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        
        push!(relerr, norm(u_MMS .- soln.u[2])/norm(u_MMS))
        push!(interr, norm(u_MMS[order:end-order+1] .- soln.u[2][order:end-order+1])/norm(u_MMS[order:end-order+1]))
    end

    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,interr=interr,conv_rate=conv_rate,npts=npts)
end




###=== MMS TESTS ===###
npts = [21,31,41,51,61,71,81,91,101,111,121,131,141,151,161,171,181,191,201]
# npts = [25,51,101,201,401,801]
# npts = [21,31,41,51,61,71,81,91,101,111,121,131,141,151]


θ = 0.5


# Solution
ũ(x,t;ωt=1.0, ωx=1.0,cx=0.0) = cos(2π*ωt*t) * sin(2π*x*ωx + cx)

# Initial condition
ũ₀(x;ωt=1.0, ωx=1.0,cx=0.0) = sin(2π*ωx*x + cx)


K = 1.0
F(x,t;ωt=1.0, ωx=1.0,cx=0.0,K=1.0) = 
        -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx) + 
            K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)
            
    
println("K=",K)


if rundirichlet
    # Dirichlet
    println("=====")
    println("Dirichlet")
    cx=1.0
    ωx=9.0
    ωt=1.0


    analytic(x,t) = ũ(x,t, ωt=ωt, ωx=ωx, cx=cx)
    IC(x) = ũ₀(x, ωt=ωt, ωx=ωx, cx=cx)
    FD(x,t) = F(x,t, ωt=ωt, ωx=ωx, cx=cx, K=K)

    BxLũ(t) = cos(2π*ωt*t) * sin(cx) #Boundary condition x=0
    BxRũ(t;Lx=1.0) = cos(2π*ωt*t) * sin(2π*Lx*ωx + cx) #Boundary condition x=Lx

    order = 2
    println("order=",order)
    O2_DirichletMMS = comp_MMS(𝒟x,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        FD,analytic,IC,order,
        k=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS = comp_MMS(𝒟x,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        FD,analytic,IC,order,
        k=K,θ=θ)

    println("ωx=",ωx,",  cx=",cx,", ωt=",ωt,", θ=",θ)
    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

    println("=====")
end

if runneumann
    # Neumann
    println("=====")
    println("Neumann")

    cx=1.0
    ωx=9.0
    ωt=1.0

    println("ωx=",ωx,"  cx=",cx)

    analytic(x,t) = ũ(x,t, ωt=ωt, ωx=ωx, cx=cx)
    IC(x) = ũ₀(x, ωt=ωt, ωx=ωx, cx=cx)
    FD(x,t) = F(x,t, ωt=ωt, ωx=ωx, cx=cx, K=K)

    BxLũ(t) =         2π*ωx * K * cos(2π*ωt*t) * cos(cx) #Boundary condition x=0
    BxRũ(t;Lx=1.0) =  2π*ωx * K * cos(2π*ωt*t) * cos(2π*Lx*ωx + cx) #Boundary condition x=Lx

    order = 2
    println("order=",order)
    O2_NeumannMMS = comp_MMS(𝒟x,npts,
        BxLũ,Neumann,BxRũ,Neumann,
        FD,analytic,IC,order,
        k=K,θ=θ)

    order = 4
    println("order=",order)
    O4_NeumannMMS = comp_MMS(𝒟x,npts,
        BxLũ,Neumann,BxRũ,Neumann,
        FD,analytic,IC,order,
        k=K,θ=θ)

    println("Order 2 Neumann convergence rates=",O2_NeumannMMS.conv_rate)
    println("Order 4 Neumann convergence rates=",O4_NeumannMMS.conv_rate)

    println("=====")
end

if runperiodic
    # Periodic
    println("=====")
    println("Periodic")

    cx=1.0
    ωx=8.0
    ωt=1.0

    println("ωx=",ωx,",  cx=",cx)

    analytic(x,t) = ũ(x,t,ωt=ωt, ωx=ωx, cx=cx)
    IC(x) = ũ₀(x, ωt=ωt, ωx=ωx, cx=cx)
    FD(x,t) = F(x,t, ωt=ωt, ωx=ωx, cx=cx, K=K)

    order = 2
    O2_PeriodicMMS = comp_MMS(𝒟x,npts,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        k=K)

    order = 4
    O4_PeriodicMMS = comp_MMS(𝒟x,npts,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        k=K)

    println("Order 2 Periodic convergence rates=",O2_PeriodicMMS.conv_rate)
    println("Order 4 Periodic convergence rates=",O4_PeriodicMMS.conv_rate)

    println("=====")
end


#=


savefig(p,".//testing//MMS//MMSTests.eps")
=#



#=
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
=#





# using Plots
# using LaTeXStrings


# pO2 = plot(O2_DirichletMMS.npts,     O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle, axis=:log)
# plot!(pO2,    O2_NeumannMMS.npts,       O2_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^2)$", markershape=:square)
# plot!(pO2,    O2_PeriodicMMS.npts,      O2_PeriodicMMS.relerr,      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$", markershape=:x)

# pO4 = plot(O4_DirichletMMS.npts,     O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:circle, axis=:log)
# plot!(pO4,    O4_NeumannMMS.npts,       O4_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^4)$", markershape=:square)
# plot!(pO4,    O4_PeriodicMMS.npts,      O4_PeriodicMMS.relerr,      label=L"Dirichlet/Periodic $\mathcal{O}(h^4)$", markershape=:x)

# plot!(pO2, [npts[2],npts[end-1]], -log.([npts[2],npts[end-1]].^2) .+ log(npts[2]^2),    linestyle=:dash, linecolor=:black,label=L"$\mathcal{O}(h^2)$")

# plot!(pO2, npts, npts.^2 , linestyle=:dash, linecolor=:black,label=L"$\mathcal{O}(h^2)$")

# p = plot(O2_DirichletMMS.npts,     O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle, axis=:log)
# plot!(p, O4_DirichletMMS.npts,     O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:circle)



#=
pO2 = plot(log.(O2_DirichletMMS.npts),     log.(O2_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
plot!(pO2,    log.(O2_NeumannMMS.npts),       log.(O2_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^2)$", markershape=:circle)
plot!(pO2,    log.(O2_PeriodicMMS.npts),      log.(O2_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^2)$", markershape=:circle)

plot!(pO2, log.([npts[2],npts[end-1]]), -log.([npts[2],npts[end-1]].^2) .+ log(npts[2]^2) .+ log(sum(O2_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dash, linecolor=:black,label=L"$\mathcal{O}(h^2)$")
=#
# savefig(pO2,".//testing//MMS//MMSTests_order2.eps")
# savefig(pO2,".//testing//MMS//MMSTests_order2.png")

#=
pO4 = plot(axis=:log,minorgrid=true)
#=
plot!(pO4,    O4_DirichletMMS.npts,     O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    O4_NeumannMMS.npts,       O4_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    O4_PeriodicMMS.npts,      O4_PeriodicMMS.relerr,      label=L"Dirichlet/Periodic $\mathcal{O}(h^4)$", markershape=:x)
=#
plot!(pO4,    log.(O4_DirichletMMS.npts),     log.(O4_DirichletMMS.relerr),     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    log.(O4_NeumannMMS.npts),       log.(O4_NeumannMMS.relerr),       label=L"Neumann $\mathcal{O}(h^4)$", markershape=:x)
plot!(pO4,    log.(O4_PeriodicMMS.npts),      log.(O4_PeriodicMMS.relerr),      label=L"Dirichlet/Periodic $\mathcal{O}(h^4)$", markershape=:x)

plot!(pO4, log.([npts[2],npts[end-1]]), 
    -log.([npts[2],npts[end-1]].^2) .+ log(npts[2]^2) .+ log(sum(O4_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dash, linecolor=:black, label=L"$\mathcal{O}(h^2)$")
plot!(pO4, log.([npts[2],npts[end-1]]), 
    -log.([npts[2],npts[end-1]].^4) .+ log(npts[2]^4) .+ log(sum(O4_DirichletMMS.relerr[1:2]/2)),
    linestyle=:dashdot, linecolor=:black, label=L"$\mathcal{O}(h^4)$")
=#
# savefig(pO4,".//testing//MMS//MMSTests_order4.eps")
# savefig(pO4,".//testing//MMS//MMSTests_order4.png")








# plot(O4_DirichletMMS.comp_soln[8].u[2],label="comp")
# plot!(O4_DirichletMMS.MMS_soln[8],label="exact")

# plot(O4_DirichletMMS.comp_soln[8].u[2] .- O4_DirichletMMS.MMS_soln[8],label="err")




#=
using DelimitedFiles

nameappend=string("K=",K)

open(string("testing/MMS/1DMMS_Tests_O2",nameappend,".csv"),"w") do io
    writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr O2_PeriodicMMS.relerr])
end
open(string("testing/MMS/1DMMS_Rates_O2",nameappend,".csv"),"w") do io
    writedlm(io,[O2_DirichletMMS.conv_rate O2_NeumannMMS.conv_rate O2_PeriodicMMS.conv_rate])
end

open(string("testing/MMS/1DMMS_Tests_O4",nameappend,".csv"),"w") do io
    writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr O4_PeriodicMMS.relerr])
end
open(string("testing/MMS/1DMMS_Rates_O4",nameappend,".csv"),"w") do io
    writedlm(io,[O4_DirichletMMS.conv_rate O4_NeumannMMS.conv_rate O4_PeriodicMMS.conv_rate])
end
=#