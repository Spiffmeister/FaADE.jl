using LinearAlgebra

# using Plots
# using LaTeXStrings
using Revise
using FaADE



rundirichlet    = false
runneumann      = false
runperiodic     = false
runrobin        = true

saverates = false

# plots = false

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
        dt_scale=0.01,t_f=5.0,k=1.0,θ=1.0)

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
            Bx0 = FaADE.SATs.SAT_Periodic(Dom.Δx,order,Left)
            BxL = FaADE.SATs.SAT_Periodic(Dom.Δx,order,Right)
        elseif BX0Type == Dirichlet
            Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left, order)
            BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right,order)
        elseif BX0Type == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,1,order)
        elseif BX0Type == Robin
            Bx0 = SAT_Robin(BoundaryX0,Dom.Δx,Left, order)
            BxL = SAT_Robin(BoundaryXL,Dom.Δx,Right,order)
        end
        BD = (Bx0,BxL)

        Δt = dt_scale*Dom.Δx

        K = k

        P = Problem1D(order,ũ₀,K,Dom,BD,F,nothing)

        println("Solving n=",Dom.n," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=θ)
        # soln = solve(P,Dom,Δt,2.1Δt,solver=:theta,θ=θ)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2]+Δt)

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        
        push!(relerr, norm(u_MMS .- soln.u[2])/norm(u_MMS))
        push!(interr, norm(u_MMS[order:end-order+1] .- soln.u[2][order:end-order+1])/norm(u_MMS[order:end-order+1]))
    end

    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,interr=interr,conv_rate=conv_rate,npts=npts)
end


function time_MMS(Dx,npts,
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
            Bx0 = FaADE.SATs.SAT_Periodic(Dom.Δx,order,Left)
            BxL = FaADE.SATs.SAT_Periodic(Dom.Δx,order,Right)
        elseif BX0Type == Dirichlet
            Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left, order)
            BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right,order)
        elseif BX0Type == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,1,order)
        end
        BD = (Bx0,BxL)

        Δt = dt_scale*Dom.Δx

        K = k

        P = Problem1D(order,ũ₀,K,Dom,BD,F,nothing)

        println("Solving n=",Dom.n," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f,solver=:theta,θ=θ)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2]+Δt)

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
# npts = [51,61,71,81,91,101,111,121,131,141,151,161,171,181,191,201]
npts = [51,61,71,81,91,101]

@show θ = 0.5


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
    # cx=0.7
    # ωx=12.5
    # ωt=1.0
    cx=1.0
    ωx=16.5
    ωt=21.0


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

    cx=0.7
    ωx=12.5
    ωt=1.0
    # cx=0.0
    # ωx=1.0
    # ωt=5.0

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
    ωx=12.0
    ωt=1.0
    # cx=1.0
    # ωx=1.0
    # ωt=12.0
    
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





if runrobin
    println("=====")
    println("Robin")

    cx=1.0
    ωx=6.5
    ωt=1.0

    analytic(x,t) = ũ(x,t,ωt=ωt, ωx=ωx, cx=cx)
    IC(x) = ũ₀(x, ωt=ωt, ωx=ωx, cx=cx)
    FD(x,t) = F(x,t, ωt=ωt, ωx=ωx, cx=cx, K=K)

    α = 1.0
    BxLũ(t) = α*cos(2π*ωt*t) * sin(cx)          - 2π*ωx*K * cos(2π*ωt*t) * cos(cx) #Boundary condition x=0
    BxRũ(t) = α*cos(2π*ωt*t) * sin(2π*ωx + cx)  + 2π*ωx*K * cos(2π*ωt*t) * cos(2π*ωx + cx) #Boundary condition x=L

    order = 2
    O2_RobinMMS = comp_MMS(𝒟x,npts,
        BxLũ,Robin,BxRũ,Robin,
        FD,analytic,IC,order,
        k=K,θ=θ)

    order = 4
    O4_RobinMMS = comp_MMS(𝒟x,npts,
        BxLũ,Robin,BxRũ,Robin,
        FD,analytic,IC,order,
        k=K,θ=θ)

    println("Order 2 Robin convergence rates=",O2_RobinMMS.conv_rate)
    println("Order 2 Robin convergence rates=",O4_RobinMMS.conv_rate)
end







if saverates
    using DelimitedFiles

    nameappend=string("spatial")

    open(string("testing/MMS/1DMMS_Tests_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr O2_PeriodicMMS.relerr O2_RobinMMS.relerr])
    end

    open(string("testing/MMS/1DMMS_Tests_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr O4_PeriodicMMS.relerr O4_RobinMMS.relerr])
    end
end




# if plots
#     using GLMakie

#     f = Figure()

#     Ax = Axis(f[1,1])
#     menu1 = Menu(f[2,1], options = ["$N" for N in npts], default="MMS")

#     on(menu1.selection) do s
#         empty!(Ax)
#         if s == "MMS"
#             pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.MMS_soln[end])
#         elseif s == "Computed"
#             pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.comp_soln[end].u[2])
#         elseif s == "Error"
#             pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.MMS_soln[end] .- O2_DirichletMMS.comp_soln[end].u[2])
#         end
#     end
#     notify(menu1.selection)
#     f
# end

# plot(O4_DirichletMMS.comp_soln[1].grid.grid,O4_DirichletMMS.comp_soln[1].u[2] .- O4_DirichletMMS.MMS_soln[1])
# plot!(O4_DirichletMMS.comp_soln[end].grid.grid,O4_DirichletMMS.comp_soln[end].u[2] .- O4_DirichletMMS.MMS_soln[end])


# plot(O4_DirichletMMS.comp_soln[end].grid.grid,O4_DirichletMMS.MMS_soln[end])
# plot!(O4_DirichletMMS.comp_soln[end].grid.grid,O4_DirichletMMS.comp_soln[end].u[2])
