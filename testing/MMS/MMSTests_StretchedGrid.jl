using LinearAlgebra

using Plots

using FaADE






###=== GLOBAL PROPS ===###
𝒟x = [0.0,1.0]
𝒟y = [0.0,1.0]

###=== MMS ===###


TestDirichlet   = false
TestNeumann     = false
TestPeriodic    = true


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

        grid_fn_x(x) = 0.05x
        grid_fn_y(y) = y

        𝒟x,𝒟y= FaADE.Grid.meshgrid(grid_fn_x.(LinRange(Dx[1],Dx[2],n)),grid_fn_y.(LinRange(Dy[1],Dy[2],n)))
        Dom = Grid2D(𝒟x,𝒟y,ymap=false)

        # Dom = Grid2D(Dx,Dy,n,n)

        # X boundaries
        if BX0Type == Periodic
            Bx0 = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Left)
            BxL = FaADE.SATs.SAT_Periodic(Dom.Δx,1,order,Right)
            By0 = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Up)
            ByL = FaADE.SATs.SAT_Periodic(Dom.Δy,2,order,Down)
        elseif BX0Type == Dirichlet
            Bx0 = FaADE.SATs.SAT_Dirichlet(BoundaryX0,Dom.Δx,Left,  order, 0.0, :Cartesian)
            BxL = FaADE.SATs.SAT_Dirichlet(BoundaryXL,Dom.Δx,Right, order, 0.0, :Cartesian)
            By0 = FaADE.SATs.SAT_Dirichlet(BoundaryY0,Dom.Δy,Up,    order, 0.0, :Cartesian)
            ByL = FaADE.SATs.SAT_Dirichlet(BoundaryYL,Dom.Δy,Down,  order, 0.0, :Cartesian)
        elseif BX0Type == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,    1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,   1,order)
            By0 = FaADE.SATs.SAT_Neumann(BoundaryY0,Dom.Δy,Up,      2,order)
            ByL = FaADE.SATs.SAT_Neumann(BoundaryYL,Dom.Δy,Down,    2,order)
        end
        BD = FaADE.Inputs.SATBoundaries(Bx0,BxL,By0,ByL)


        Δt = dt_scale*Dom.Δx^2
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
    cx=1.0
    cy=0.0
    ωx=7.5
    ωy=5.0
    ωt=1.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    BxLũ(y,t)           = cos(2π*ωt*t) * sin(cx) * sin(2π*y*ωy + cy) #Boundary condition x=0
    BxRũ(y,t;Lx=1.0)    = cos(2π*ωt*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy) #Boundary condition x=Lx
    ByLũ(x,t)           = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(cy) #Boundary condition y=0
    ByRũ(x,t;Ly=1.0)    = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*Ly*ωy + cy) #Boundary condition y=Ly

    order = 2
    println("order=",order)
    O2_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        ByLũ,Dirichlet,ByRũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Dirichlet,BxRũ,Dirichlet,
        ByLũ,Dirichlet,ByRũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

    println("=====")
end



# Neumann
if TestNeumann
    println("=====")
    println("Neumann")

    cx=1.0
    cy=0.0
    ωx=7.5
    ωy=6.0
    ωt=1.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K=K)

    BxLũ(y,t) =         2π*ωx * K * cos(2π*ωt*t) * cos(cx)             * sin(2π*y*ωy + cy) #Boundary condition x=0
    BxRũ(y,t;Lx=1.0) =  2π*ωx * K * cos(2π*ωt*t) * cos(2π*Lx*ωx + cx)  * sin(2π*y*ωy + cy) #Boundary condition x=Lx
    ByLũ(x,t) =         2π*ωy * K * cos(2π*ωt*t) * sin(2π*x*ωx + cx)   * cos(cy) #Boundary condition y=0
    ByRũ(x,t;Ly=1.0) =  2π*ωy * K * cos(2π*ωt*t) * sin(2π*x*ωx + cx)   * cos(2π*Ly*ωy + cy) #Boundary condition y=Ly

    order = 2
    println("order=",order)
    O2_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Neumann,BxRũ,Neumann,
        ByLũ,Neumann,ByRũ,Neumann,
        FD,analytic,IC,order,
        kx=K, ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
        BxLũ,Neumann,BxRũ,Neumann,
        ByLũ,Neumann,ByRũ,Neumann,
        FD,analytic,IC,order,
        kx=K, ky=K,θ=θ)

    println("Order 2 Neumann convergence rates=",O2_NeumannMMS.conv_rate)
    println("Order 4 Neumann convergence rates=",O4_NeumannMMS.conv_rate)

    println("=====")
end


# Periodic
if TestPeriodic
    println("=====")
    println("Periodic")

    cx=0.0
    cy=0.0
    ωx=1.0
    ωy=1.0
    ωt=9.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K=K)

    order = 2
    O2_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
        nothing,Periodic,nothing,Periodic,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K, ky=K,θ=θ)

    order = 4
    O4_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
        nothing,Periodic,nothing,Periodic,
        nothing,Periodic,nothing,Periodic,
        FD,analytic,IC,order,
        kx=K, ky=K,θ=θ)

    println("Order 2 Periodic convergence rates=",O2_PeriodicMMS.conv_rate)
    println("Order 4 Periodic convergence rates=",O4_PeriodicMMS.conv_rate)

    println("=====")
end


if TestDirichlet == TestNeumann == TestPeriodic
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
    jldsave("testing/MMS/FullMMS_stretchgrid.jld2";O2Conv,O4Conv)


    using DelimitedFiles

    nameappend=string("timeconv")

    open(string("testing/MMS/MMS_Tests_O2_stretchgrid",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr O2_PeriodicMMS.relerr])
    end
    open(string("testing/MMS/MMS_Rates_O2_stretchgrid",nameappend,".csv"),"w") do io
        writedlm(io,[O2_DirichletMMS.conv_rate O2_NeumannMMS.conv_rate O2_PeriodicMMS.conv_rate])
    end

    open(string("testing/MMS/MMS_Tests_O4_stretchgrid",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr O4_PeriodicMMS.relerr])
    end
    open(string("testing/MMS/MMS_Rates_O4_stretchgrid",nameappend,".csv"),"w") do io
        writedlm(io,[O4_DirichletMMS.conv_rate O4_NeumannMMS.conv_rate O4_PeriodicMMS.conv_rate])
    end
end







