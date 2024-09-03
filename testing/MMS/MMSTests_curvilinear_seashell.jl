using LinearAlgebra
using Revise
using FaADE

inspect = true



###=== MMS ===###

TestDirichlet   = false
TestNeumann     = true
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





function comp_MMS(D,npts,
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
        
        Dom = D(n)

        # X boundaries
        if BX0Type == Dirichlet
            Bx0 = SAT_Dirichlet(BoundaryX0,Dom.Δx,Left, order, Dom.Δy, :Curvilinear)
            BxL = SAT_Dirichlet(BoundaryX0,Dom.Δx,Right,order, Dom.Δy, :Curvilinear)
            By0 = SAT_Dirichlet(BoundaryY0,Dom.Δy,Down, order, Dom.Δx, :Curvilinear)
            ByL = SAT_Dirichlet(BoundaryYL,Dom.Δy,Up,   order, Dom.Δx, :Curvilinear)
        elseif BX0Type == Neumann
            Bx0 = SAT_Neumann(BoundaryX0,Dom.Δx,Left,    order, Dom.Δy, :Curvilinear)
            BxL = SAT_Neumann(BoundaryXL,Dom.Δx,Right,   order, Dom.Δy, :Curvilinear)
            By0 = SAT_Neumann(BoundaryY0,Dom.Δy,Down,    order, Dom.Δx, :Curvilinear)
            ByL = SAT_Neumann(BoundaryYL,Dom.Δy,Up,      order, Dom.Δx, :Curvilinear)
        end
        BD = (Bx0,BxL,By0,ByL)

        Δt = dt_scale*min(Dom.Δx,Dom.Δy)^2
        t_f = 1e-1
        nt = round(t_f/Δt)
        Δt = t_f/nt

        P = Problem2D(order,ũ₀,kx,ky,Dom,BD,F,nothing)

        println("Solving n=",Dom.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f)
        # soln = solve(P,Dom,Δt,Δt*100)

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
# npts = collect(21:10:101)
npts = [21,41,81]

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
    ωx=3.5
    ωy=2.5
    ωt=3.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    Bũ(X,t)           = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0

    
    Dx(n) = Grid2D(u->u*[0.5, 0.0] - [0.0,0.0],
        v->v*[cos(3π/4),sin(3π/4)],
        v->v*[cos(π/4),sin(π/4)] + [0.5,0.0],
        u->[cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)] + 
            u*[0.5,0.0], # top of D1
        n,n)

    order = 2
    println("order=",order)
    O2_DirichletMMS_x = comp_MMS(Dx,npts,
        Bũ,Dirichlet,Bũ,Dirichlet,
        Bũ,Dirichlet,Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS_x = comp_MMS(Dx,npts,
        Bũ,Dirichlet,Bũ,Dirichlet,
        Bũ,Dirichlet,Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS_x.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS_x.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS_x.relerr)
    println("Order 4 relative error=",O4_DirichletMMS_x.relerr)



    Dy(n) = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
        v->[0.0, v/2 - 0.25],
        v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
            (1-v)*[0.0,0.25] + v*[0.0,0.25],
        u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
        n,n)

    order = 2
    println("order=",order)
    O2_DirichletMMS_y = comp_MMS(Dx,npts,
        Bũ,Dirichlet,Bũ,Dirichlet,
        Bũ,Dirichlet,Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS_y = comp_MMS(Dx,npts,
        Bũ,Dirichlet,Bũ,Dirichlet,
        Bũ,Dirichlet,Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS_y.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS_y.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS_y.relerr)
    println("Order 4 relative error=",O4_DirichletMMS_y.relerr)
    

    println("=====")
end





# Dirichlet
if TestNeumann
    println("=====")
    println("Neumann")
    cx=0.0
    cy=0.0
    ωx=1.0
    ωy=1.0
    ωt=1.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    # Bũ(X,t)           = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0
    DBũx(X,t)          = 2π*ωx * K * cos(2π*ωt*t) * cos(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0
    DBũy(X,t)          = 2π*ωy * K * cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * cos(2π*X[2]*ωy + cy) #Boundary condition y=0

    
    Dx(n) = Grid2D(u->u*[0.5, 0.0] - [0.0,0.0],
        v->v*[cos(3π/4),sin(3π/4)],
        v->v*[cos(π/4),sin(π/4)] + [0.5,0.0],
        u->[cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)] + 
            u*[0.5,0.0], # top of D1
        n,n)

    order = 2
    println("order=",order)
    O2_NeumannMMS_x = comp_MMS(Dx,npts,
        DBũx,Neumann,DBũx,Neumann,
        DBũy,Neumann,DBũy,Neumann,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    # order = 4
    # println("order=",order)
    # O4_NeumannMMS_x = comp_MMS(Dx,npts,
    #     DBũx,Neumann,DBũx,Neumann,
    #     DBũy,Neumann,DBũy,Neumann,
    #     FD,analytic,IC,order,
    #     kx=K,ky=K,θ=θ)

    println("Order 2 Neumann convergence rates=",O2_NeumannMMS_x.conv_rate)
    println("Order 4 Neumann convergence rates=",O4_NeumannMMS_x.conv_rate)

    println("Order 2 relative error=",O2_NeumannMMS_x.relerr)
    println("Order 4 relative error=",O4_NeumannMMS_x.relerr)



    # Dy(n) = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
    #     v->[0.0, v/2 - 0.25],
    #     v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
    #         (1-v)*[0.0,0.25] + v*[0.0,0.25],
    #     u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
    #     n,n)

    # order = 2
    # println("order=",order)
    # O2_NeumannMMS_y = comp_MMS(Dx,npts,
    #     DBũx,Neumann,DBũx,Neumann,
    #     DBũy,Neumann,DBũy,Neumann,
    #     FD,analytic,IC,order,
    #     kx=K,ky=K,θ=θ)

    # order = 4
    # println("order=",order)
    # O4_NeumannMMS_y = comp_MMS(Dx,npts,
    #     DBũx,Neumann,DBũx,Neumann,
    #     DBũy,Neumann,DBũy,Neumann,
    #     FD,analytic,IC,order,
    #     kx=K,ky=K,θ=θ)

    # println("Order 2 Neumann convergence rates=",O2_NeumannMMS_y.conv_rate)
    # println("Order 4 Neumann convergence rates=",O4_NeumannMMS_y.conv_rate)

    # println("Order 2 relative error=",O2_NeumannMMS_y.relerr)
    # println("Order 4 relative error=",O4_NeumannMMS_y.relerr)
    

    println("=====")
end











if SaveTests
    using DelimitedFiles

    nameappend=string("timeconv")

    open(string("testing/MMS/out/MMS_Tests_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS_x.relerr O2_DirichletMMS_y.relerr])# O2_NeumannMMS.relerr])
    end

    open(string("testing/MMS/out/MMS_Tests_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS_x.relerr O4_DirichletMMS_y.relerr])# O4_NeumannMMS.relerr])
    end
end



