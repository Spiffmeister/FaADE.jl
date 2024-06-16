using LinearAlgebra

using FaADE

inspect = true



###=== MMS ===###

TestDirichlet   = true
SaveTests       = false


# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::GridMultiBlock,t::Float64)
    D1 = grid.Grids[1]
    D2 = grid.Grids[2]
    u_MMS1 = zeros(size(D1))
    u_MMS2 = zeros(size(D2))
    for I in eachindex(D1)
        u_MMS1[I] = MMS(D1.gridx[I],D1.gridy[I],t)
    end
    for I in eachindex(D2)
        u_MMS2[I] = MMS(D2.gridx[I],D2.gridy[I],t)
    end
    return u_MMS1, u_MMS2
end





function comp_MMS_y(npts,
        Bxy,BType,
        F,ũ,ũ₀,order;
        dt_scale=1.0,t_f=0.1,kx=1.0,ky=kx,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    # relerr = []
    relerr = [zeros(length(npts)),zeros(length(npts))]

    # Loop
    # for n in npts
    for I in eachindex(npts)
        n = npts[I]

        D1 = Grid2D(u->u*[cos(7π/4), sin(7π/4)] - [0.0,0.25],
                v->[0.0, v/2 - 0.25],
                v->[cos(v*(9π/4 - 7π/4) + 7π/4), sin(v*(9π/4 - 7π/4) + 7π/4)] - 
                    (1-v)*[0.0,0.25] + v*[0.0,0.25],
                u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # top of D1
            n,n)
    
        D2 = Grid2D(u->u*[cos(π/4), sin(π/4)] + [0.0, 0.25], # bottom - same as D1 top -- done
                v->[cos(v*(3π/4 - 5π/4) + 5π/4), sin(v*(3π/4 - 5π/4) + 5π/4)] + v*[0.0, 0.5] - 
                        [cos(5π/4), sin(5π/4) - 0.25], # left - shifted up 0.25 -- done
                v->v*[0.0, 0.5] + [cos(π/4), sin(π/4) + 0.25], # shifted to top of D1 -- done
                u->u*[cos(3π/4), sin(3π/4)] + [cos(π/4), sin(π/4) + 0.75], # top - shifted up 0.25
            n,n)

        joints = ((Joint(2,Up),), (Joint(1,Down),),)
        Dom = GridMultiBlock((D1,D2),joints)

        if BType == Dirichlet
            Dl1 = SAT_Dirichlet(Bxy, D1.Δx, Left, order, D1.Δy, :Curvilinear) # Block 3 BCs
            Dr1 = SAT_Dirichlet(Bxy, D1.Δx, Right,order, D1.Δy, :Curvilinear) # 
            Dd1 = SAT_Dirichlet(Bxy, D1.Δy, Down, order, D1.Δx, :Curvilinear) # Block 5 BCs
            
            Dr2 = SAT_Dirichlet(Bxy, D2.Δx, Right,order, D2.Δy, :Curvilinear) # Block 2 BCs
            Dl2 = SAT_Dirichlet(Bxy, D2.Δx, Left, order, D2.Δy, :Curvilinear) # Block 5 BCs
            Du2 = SAT_Dirichlet(Bxy, D2.Δy, Up,   order, D2.Δx, :Curvilinear) # Block 3 BCs

            BD = Dict(1 => (Dl1,Dr1,Dd1), 2 => (Dr2,Dl2,Du2))
        elseif BType == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,    1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,   1,order)
        end

        Δt = dt_scale*min(D1.Δx,D1.Δy)^2
        t_f = 1e-1
        nt = round(t_f/Δt)
        Δt = t_f/nt

        P = Problem2D(order,ũ₀,kx,ky,Dom,BD,F,nothing)

        println("Solving n=",D1.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2])

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        # push!(relerr, norm(u_MMS[1] .- soln.u[2][1])/norm(u_MMS[1]) + norm(u_MMS[2] .- soln.u[2][2])/norm(u_MMS[2]))
        @show relerr[1][I] = norm(u_MMS[1] .- soln.u[2][1])/norm(u_MMS[1])
        @show relerr[2][I] = norm(u_MMS[2] .- soln.u[2][2])/norm(u_MMS[2])
    end

    # conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))
    conv_rate = log2.(relerr[1][1:end-1] ./ relerr[1][2:end], relerr[2][1:end-1] ./ relerr[2][2:end])

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end




function comp_MMS_x(npts,
        Bxy,BType,
        F,ũ,ũ₀,order;
        dt_scale=1.0,t_f=0.1,kx=1.0,ky=kx,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    relerr = []

    relerr = [zeros(length(npts)),zeros(length(npts))]

    # Loop
    for I in eachindex(npts)
        n = npts[I]

        D1 = Grid2D(u->u*[0.5, 0.0] - [0.0,0.0],
            v->v*[cos(3π/4),sin(3π/4)] + [0.0,0.0],
            v->v*[cos(π/4),sin(π/4)] + [0.5,0.0],
            u->[cos(u*(π/4 - 3π/4) + 3π/4), sin(u*(π/4 - 3π/4) + 3π/4)] + 
                u*[0.5,0.0], # top of D1
            n,n)

        D2 = Grid2D(u->[cos(u*(-5π/4 + 7π/4) + 5π/4), sin(u*(-5π/4 + 7π/4) + 5π/4)] - 
                [cos(7π/4), sin(7π/4)] + u*[0.5,0.0] - [0.5,0.0], # bottom - same as D1 top -- done
            v->v*[cos(π/4),sin(π/4)] + [cos(5π/4)-cos(7π/4)-0.5,sin(5π/4)-sin(7π/4)], # shifted to top of D1 -- done
            v->v*[cos(3π/4),sin(3π/4)] + [0.0,0.0], # left - shifted up 0.25 -- done
            u->u*[0.5, 0.0] + [cos(3π/4) - 0.5, sin(3π/4)], # top - shifted up 0.25
            n,n)

        joints = ((Joint(2,Left),),(Joint(1,Right),),)

        Dom = GridMultiBlock((D1,D2),joints)



        if BType == Dirichlet
            Dr1 = SAT_Dirichlet(Bxy, D1.Δx, Right,  order, D1.Δy, :Curvilinear) #
            Dd1 = SAT_Dirichlet(Bxy, D1.Δy, Down,   order, D1.Δx, :Curvilinear) #
            Du1 = SAT_Dirichlet(Bxy, D1.Δy, Up,     order, D1.Δx, :Curvilinear) #

            Dl2 = SAT_Dirichlet(Bxy, D2.Δx, Left, order, D2.Δy, :Curvilinear) #
            Du2 = SAT_Dirichlet(Bxy, D2.Δy, Up,   order, D2.Δx, :Curvilinear) #
            Dd2 = SAT_Dirichlet(Bxy, D2.Δy, Down, order, D2.Δx, :Curvilinear) #

            BD = Dict(1 => (Dr1,Dd1,Du1), 2 => (Dl2,Du2,Dd2))
        elseif BType == Neumann
            Bx0 = FaADE.SATs.SAT_Neumann(BoundaryX0,Dom.Δx,Left,    1,order)
            BxL = FaADE.SATs.SAT_Neumann(BoundaryXL,Dom.Δx,Right,   1,order)
        end

        Δt = dt_scale*min(D1.Δx,D1.Δy)^2
        t_f = 1e-1
        nt = round(t_f/Δt)
        Δt = t_f/nt

        P = Problem2D(order,ũ₀,kx,ky,Dom,BD,F,nothing)

        println("Solving n=",D1.nx," case with Δt=",Δt)
        soln = solve(P,Dom,Δt,t_f)

        u_MMS = generate_MMS(ũ,Dom,soln.t[2])

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        # push!(relerr, ((norm(u_MMS[1] .- soln.u[2][1]))/norm(u_MMS[1]) + 
                        # (norm(u_MMS[2] .- soln.u[2][2]))/norm(u_MMS[2])))
        relerr[1][I] = norm(u_MMS[1] .- soln.u[2][1])/norm(u_MMS[1])
        relerr[2][I] = norm(u_MMS[2] .- soln.u[2][2])/norm(u_MMS[2])
    end

    # conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (2npts[1:end-1].-1))./(1 ./ (2npts[2:end].-1) ))
    conv_rate = (relerr[1][1:end-1] ./ relerr[1][2:end], relerr[2][1:end-1] ./ relerr[2][2:end])

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end



########################################################################################
########################################################################################
# CALLING
########################################################################################
########################################################################################


###=== MMS TESTS ===###
# npts = collect(21:10:81)
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
    cx=1.0
    cy=0.5
    ωx=4.5
    ωy=3.5
    ωt=3.0

    println("ωx=",ωx,"  ωy=",ωy,",  cx=",cx,",  cy=",cy,", ωt=",ωt," θ=",θ)

    analytic(x,y,t) = ũ(x,y,t, ωt=ωt , ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    IC(x,y) = ũ₀(x,y, ωx=ωx, cx=cx, ωy=ωy, cy=cy)
    FD(X,t) = F(X[1],X[2],t, ωt=ωt, ωx=ωx, cx=cx, ωy=ωy, cy=cy, K = K)

    Bũ(X,t)           = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0

    order = 2
    println("order=",order)
    O2_DirichletMMS_y = comp_MMS_y(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS_y = comp_MMS_y(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS_y.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS_y.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS_y.relerr)
    println("Order 4 relative error=",O4_DirichletMMS_y.relerr)



    order = 2
    println("order=",order)
    O2_DirichletMMS_x = comp_MMS_x(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS_x = comp_MMS_x(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS_x.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS_x.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS_x.relerr)
    println("Order 4 relative error=",O4_DirichletMMS_x.relerr)
    

    println("=====")
end







if SaveTests
    using DelimitedFiles

    nameappend=string("spatial")

    open(string("testing/MMS/out/MMS_DoubleSeashell_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS_x.relerr O2_DirichletMMS_y.relerr])# O2_NeumannMMS.relerr])
    end

    open(string("testing/MMS/out/MMS_DoubleSeashell_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS_x.relerr O4_DirichletMMS_y.relerr])# O4_NeumannMMS.relerr])
    end
end



