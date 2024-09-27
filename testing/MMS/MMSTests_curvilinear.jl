using LinearAlgebra

using FaADE

inspect = true




###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

###=== MMS ===###


TestDirichlet   = false
TestRobin       = true
SaveTests       = false


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
        
        Rin = [3e-1]; Zin=[3e-1]
        Rout = [6e-1]; Zout=[6e-1]

        inner = FaADE.Grid.Torus(Rin,Zin,[1],[0])
        outer = FaADE.Grid.Torus(Rout,Zout,[1],[0])

        X,Y = FaADE.Grid.meshgrid(inner,outer,0.0,n,n)
        Dom = Grid2D(X,Y,periodicy=true)



        # X boundaries
        By0 = SAT_Periodic(Dom.Δy,order,Up,    Dom.Δx,:Curvilinear)
        ByL = SAT_Periodic(Dom.Δy,order,Down,  Dom.Δx,:Curvilinear)
        if BX0Type == Dirichlet
            Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left,  order, Dom.Δy, :Curvilinear)
            BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right, order, Dom.Δy, :Curvilinear)
        elseif BX0Type == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,    1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,   1,order)
        elseif BX0Type == Robin
            Bx0 = SAT_Robin(BoundaryX0,Dom.Δx,Left, order, Δy = Dom.Δy, coord=:Curvilinear)
            BxL = SAT_Robin(BoundaryXL,Dom.Δx,Right,order, Δy = Dom.Δy, coord=:Curvilinear)
        end
        BD = (Bx0,BxL,By0,ByL)


        Δt = dt_scale*min(Dom.Δx,Dom.Δy)^2
        t_f = 1e-2
        nt = round(t_f/Δt)
        Δt = t_f/nt

        # Kx(x,y) = kx
        # Ky(x,y) = ky

        P = Problem2D(order,ũ₀,kx,ky,Dom,BD,F,nothing)

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
npts = collect(21:10:101)

θ = 0.5

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

# Dirichlet
if TestDirichlet
    println("=====")
    println("Dirichlet")
    cx=0.0
    cy=0.0
    ωx=5.0
    ωy=5.0
    ωt=1.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    BxLũ(X,t)           = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0
    BxRũ(X,t;Lx=1.0)    = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=Lx

    order = 2
    println("order=",order)
    O2_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS.relerr)
    println("Order 4 relative error=",O4_DirichletMMS.relerr)

    # pD = plot(axis=:log,minorgrid=true)
    # plot!(pD,  O2_DirichletMMS.npts,   O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
    # plot!(pD,  O2_DirichletMMS.npts,   O2_DirichletMMS.npts.^2,     label=L"$\mathcal{O}(h^2)$", markershape=:circle)
    # plot!(pD,  O4_DirichletMMS.npts,   O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:circle)
    # plot!(pD,  O4_DirichletMMS.npts,   O4_DirichletMMS.npts.^4,     label=L"$\mathcal{O}(h^4)$", markershape=:circle)    
    # savefig(pD,"2DMMSDirichlet.png")


    println("=====")
end



if TestRobin
    println("=====")
    println("Dirichlet")
    cx=0.0
    cy=0.0
    ωx=1.0
    ωy=1.0
    ωt=1.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    α = 1.0
    BxLũ(X,t)           = cos(2π*ωt*t) * sin(2π*X[2]*ωy + cy) * (α*sin(2π*X[1]*ωx + cx) - 2π*ωx*K*cos(2π*X[1]*ωx + cx))
    BxRũ(X,t;Lx=1.0)    = cos(2π*ωt*t) * sin(2π*X[2]*ωy + cy) * (α*sin(2π*X[1]*ωx + cx) + 2π*ωx*K*cos(2π*X[1]*ωx + cx))

    order = 2
    println("order=",order)
    O2_RobinMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Robin,BxRũ,Robin,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_RobinMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Robin,BxRũ,Robin,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Robin convergence rates=",O2_RobinMMS.conv_rate)
    println("Order 4 Robin convergence rates=",O4_RobinMMS.conv_rate)

    println("Order 2 relative error=",O2_RobinMMS.relerr)
    println("Order 4 relative error=",O4_RobinMMS.relerr)

    # pD = plot(axis=:log,minorgrid=true)
    # plot!(pD,  O2_DirichletMMS.npts,   O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$", markershape=:circle)
    # plot!(pD,  O2_DirichletMMS.npts,   O2_DirichletMMS.npts.^2,     label=L"$\mathcal{O}(h^2)$", markershape=:circle)
    # plot!(pD,  O4_DirichletMMS.npts,   O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$", markershape=:circle)
    # plot!(pD,  O4_DirichletMMS.npts,   O4_DirichletMMS.npts.^4,     label=L"$\mathcal{O}(h^4)$", markershape=:circle)    
    # savefig(pD,"2DMMSDirichlet.png")


    println("=====")
end





if SaveTests
    O2Conv = (n=npts,
        conv_D = O2_DirichletMMS.conv_rate,
        # conv_N = O2_NeumannMMS.conv_rate,
        # conv_P = O2_PeriodicMMS.conv_rate,
        relerr_D = O2_DirichletMMS.relerr,
        # relerr_N = O2_NeumannMMS.relerr,
        # relerr_P = O2_PeriodicMMS.relerr
        )

    O4Conv = (n=npts,
        conv_D = O4_DirichletMMS.conv_rate,
        # conv_N = O4_NeumannMMS.conv_rate,
        # conv_P = O4_PeriodicMMS.conv_rate,
        relerr_D = O4_DirichletMMS.relerr,
        # relerr_N = O4_NeumannMMS.relerr,
        # relerr_P = O4_PeriodicMMS.relerr
        )

    using JLD2
    jldsave("testing/MMS/FullMMS.jld2";O2Conv,O4Conv)


    using DelimitedFiles

    nameappend=string("timeconv")

    open(string("testing/MMS/MMS_Tests_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr])
    end
    open(string("testing/MMS/MMS_Rates_O2",nameappend,".csv"),"w") do io
        writedlm(io,[O2_DirichletMMS.conv_rate O2_NeumannMMS.conv_rate])
    end

    open(string("testing/MMS/MMS_Tests_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr])
    end
    open(string("testing/MMS/MMS_Rates_O4",nameappend,".csv"),"w") do io
        writedlm(io,[O4_DirichletMMS.conv_rate O4_NeumannMMS.conv_rate])
    end
end





if inspect
    using GLMakie

    f = Figure()

    Ax = Axis3(f[1,1])
    menu1 = Menu(f[2,1], options = ["MMS","Computed","Error"], default="MMS")

    on(menu1.selection) do s
        empty!(Ax)
        if s == "MMS"
            pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.MMS_soln[end])
        elseif s == "Computed"
            pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.comp_soln[end].u[2])
        elseif s == "Error"
            pobj = surface!(Ax, O2_DirichletMMS.grids[end].gridx, O2_DirichletMMS.grids[end].gridy, O2_DirichletMMS.MMS_soln[end] .- O2_DirichletMMS.comp_soln[end].u[2])
        end
    end
    notify(menu1.selection)
    f
end