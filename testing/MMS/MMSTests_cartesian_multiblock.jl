"""
MMS tests in a multiblock domain of 4 blocks stacked like
12
34
"""
using LinearAlgebra
using Revise
using FaADE





TestDirichlet   = true
TestNeumann     = false
SaveTests       = false

# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::GridMultiBlock,t::Float64)
    u_MMS = [zeros(size(grid.Grids[I])) for I in eachindex(grid.Grids)]

    for I in eachindex(grid.Grids)
        for J in eachindex(grid.Grids[I])
            u_MMS[I][J] = MMS(grid.Grids[I][J]...,t)
        end
    end

    return u_MMS
end




function comp_MMS(npts,
    Bxy,BType,
    F,ũ,ũ₀,order;
    dt_scale=1.0,t_f=0.1,kx=1.0,ky=kx,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    # relerr = []
    relerr = zeros(length(npts))

    # Loop
    # for n in npts
    for I in eachindex(npts)
        n = npts[I]

        D1 = Grid2D([0.0,0.5],[0.0,0.5],n,n)
        D2 = Grid2D([0.5,1.0],[0.0,0.5],n,n)
        D3 = Grid2D([0.0,0.5],[0.5,1.0],n,n)
        D4 = Grid2D([0.5,1.0],[0.5,1.0],n,n)

        joints = ((Joint(2,Right),Joint(3,Up)),
            (Joint(1,Left),Joint(4,Up)),
            (Joint(4,Right),Joint(1,Down)),
            (Joint(3,Left),Joint(2,Down)))

        Dom = GridMultiBlock((D1,D2,D3,D4),joints)


        if BType == Dirichlet
            
            # Block 1 BCs
            Dl1 = FaADE.SATs.SAT_Dirichlet(Bxy["Bx0"],D1.Δx,Left,   order)
            Dd1 = FaADE.SATs.SAT_Dirichlet(Bxy["By0"],D1.Δy,Down,   order)
            # Block 2 BCs
            Dd2 = FaADE.SATs.SAT_Dirichlet(Bxy["By0"],D2.Δy,Down,   order)
            Dr2 = FaADE.SATs.SAT_Dirichlet(Bxy["Bxl"],D2.Δx,Right,  order) # Bottom of missing block
            # Block 3 BCs
            Dl3 = FaADE.SATs.SAT_Dirichlet(Bxy["Bx0"],D3.Δx,Left,   order)
            Du3 = FaADE.SATs.SAT_Dirichlet(Bxy["Byl"],D3.Δy,Up,     order)
            # Block 4 BCs
            Du4 = FaADE.SATs.SAT_Dirichlet(Bxy["Byl"],D4.Δy,Up,     order)
            Dr4 = FaADE.SATs.SAT_Dirichlet(Bxy["Bxl"],D4.Δx,Right,  order) # Left of missing block

            BD= Dict(1 => (Dl1,Dd1), 2=> (Dr2,Dd2), 3 => (Dl3,Du3), 4 => (Dr4,Du4))


        elseif BType == Neumann

            # Block 1 BCs
            Nl1 = FaADE.SATs.SAT_Neumann(Bxy["Bx0"],D1.Δx,Left,   order)
            Nd1 = FaADE.SATs.SAT_Neumann(Bxy["By0"],D1.Δy,Down,   order)
            # Block 2 BCs
            Nd2 = FaADE.SATs.SAT_Neumann(Bxy["By0"],D2.Δy,Down,   order)
            Nr2 = FaADE.SATs.SAT_Neumann(Bxy["Bxl"],D2.Δx,Right,  order) # Bottom of missing block
            # Block 3 BCs
            Nl3 = FaADE.SATs.SAT_Neumann(Bxy["Bx0"],D3.Δx,Left,   order)
            Nu3 = FaADE.SATs.SAT_Neumann(Bxy["Byl"],D3.Δy,Up,     order)
            # Block 4 BCs
            Nu4 = FaADE.SATs.SAT_Neumann(Bxy["Byl"],D4.Δy,Up,     order)
            Nr4 = FaADE.SATs.SAT_Neumann(Bxy["Bxl"],D4.Δx,Right,  order) # Left of missing block

            BD= Dict(1 => (Nl1,Nd1), 2=> (Nr2,Nd2), 3 => (Nl3,Nu3), 4 => (Nr4,Nu4))
        end


        Δt = dt_scale*min(D1.Δx,D1.Δy)^2
        t_f = 1e-1
        nt = round(t_f/Δt)
        Δt = t_f/nt

        P = Problem2D(order,ũ₀,K,K,Dom,BD,source=F)

        println("Solving n=",n," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2])

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)

        tmpnume = 0.0
        tmpdeno = 0.0

        for J in eachindex(Dom.Grids)
            tmpnume += norm(soln.u[2][J] .- u_MMS[J])
            tmpdeno += norm(u_MMS[J])
        end
        relerr[I] = tmpnume/tmpdeno

    end

    # conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))
    # conv_rate = log2.(relerr[1][1:end-1] ./ relerr[1][2:end], relerr[2][1:end-1] ./ relerr[2][2:end])
    @show relerr
    @show conv_rate = log2.(relerr[1:end-1] ./ relerr[2:end])

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end














########################################################################################
########################################################################################
# CALLING
########################################################################################
########################################################################################


###=== MMS TESTS ===###

npts = [21,41,81]
# npts = 21:10:61
θ = 0.5


K = 1.0



exact(x,y,t;
    ωt=1.0,ωx=1.0,ωy=1.0,cx=0.0,cy=0.0) = cos(2π*ωt*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)

u₀(x,y;ωt=1.0,ωx=1.0,ωy=1.0,cx=0.0,cy=0.0) = exact(x,y,0.0,ωt=ωt,ωx=ωx,ωy=ωy,cx=cx,cy=cy)

F(x,y,t;
    ωt=1.0,ωx=1.0,ωy=1.0,cx=0.0,cy=0.0,K=1.0) = begin
    -2π*ωt*sin(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωx^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
    K * 4π^2 * ωy^2 * cos(2π*ωt*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy)
end




if TestDirichlet
    println("=====")
    println("Dirichlet")
    cx=1.0
    cy=0.5
    ωx=16.5
    ωy=13.0
    ωt=2.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = exact(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = u₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    Bx0(X,t) = cos(2π*ωt*t) * sin(cx)               * sin(2π*ωy*X[2] + cy)  #Boundary condition x=0
    Bxl(X,t) = cos(2π*ωt*t) * sin(2π*ωx + cx)       * sin(2π*X[2]*ωy + cy)  #Boundary condition x=Lx
    By0(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx)  * sin(cy)               #Boundary condition y=0
    Byl(X,t) = cos(2π*ωt*t) * sin(2π*X[1]*ωx + cx)  * sin(2π*ωy + cy)       #Boundary condition y=Ly
    

    Bxy = Dict(
        "Bx0" => Bx0,
        "Bxl" => Bxl,
        "By0" => By0,
        "Byl" => Byl
    )

    order = 2
    println("order=",order)
    O2_DirichletMMS = comp_MMS(npts,
        Bxy,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)
    
    order = 4
    println("order=",order)
    O4_DirichletMMS = comp_MMS(npts,
        Bxy,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)
end





if TestNeumann
    println("=====")
    println("Neumann")
    cx=0.0
    cy=0.0
    ωx=0.0
    ωy=0.0
    ωt=2.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = exact(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = u₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    Bx0(X,t) = cos(2π*ωt*t) * 2π*ωx * cos(cx)           * sin(2π*ωy*X[2] + cy)      #Boundary condition x=0
    Bxl(X,t) = cos(2π*ωt*t) * 2π*ωx * cos(2π*ωx + cx)   * sin(2π*X[2]*ωy + cy)      #Boundary condition x=Lx
    By0(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx)      * 2π*ωy * cos(cy)           #Boundary condition y=0
    Byl(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx)      * 2π*ωy * cos(2π*ωy + cy)   #Boundary condition y=Ly
    

    Bxy = Dict(
        "Bx0" => Bx0,
        "Bxl" => Bxl,
        "By0" => By0,
        "Byl" => Byl
    )

    order = 2
    println("order=",order)
    O2_NeumannMMS = comp_MMS(npts,
        Bxy,Neumann,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)
    
    order = 4
    println("order=",order)
    O4_NeumannMMS = comp_MMS(npts,
        Bxy,Neumann,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)
end













if SaveTests
    using DelimitedFiles

    nameappend=string("conv")

    open(string("./2D/data/MMS_MV4_Tests_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS.relerr O2_NeumannMMS.relerr])
    end

    open(string("./2D/data/MMS_MV4_Tests_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS.relerr O4_NeumannMMS.relerr])
    end
end






#====== 8 Block wedge ======#
#   ---------
#   | 6 | 7 |
#   -------------
#   | 4 |   | 5 |
#   -------------
#   | 1 | 2 | 3 |
#   -------------   
#




