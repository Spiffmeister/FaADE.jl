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
    D3 = grid.Grids[3]
    D4 = grid.Grids[4]
    D5 = grid.Grids[5]
    u_MMS1 = zeros(size(D1))
    u_MMS2 = zeros(size(D1))
    u_MMS3 = zeros(size(D1))
    u_MMS4 = zeros(size(D1))
    u_MMS5 = zeros(size(D1))
    for I in eachindex(D1)
        u_MMS1[I] = MMS(D1.gridx[I],D1.gridy[I],t)
        u_MMS2[I] = MMS(D2.gridx[I],D2.gridy[I],t)
        u_MMS3[I] = MMS(D3.gridx[I],D3.gridy[I],t)
        u_MMS4[I] = MMS(D4.gridx[I],D4.gridy[I],t)
        u_MMS5[I] = MMS(D5.gridx[I],D5.gridy[I],t)
    end
    return u_MMS1, u_MMS2, u_MMS3, u_MMS4, u_MMS5
end





function comp_MMS(npts,
        Bxy,BType,
        F,ũ,ũ₀,order;
        dt_scale=1.0,t_f=0.1,kx=1.0,ky=kx,θ=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    # relerr = []
    relerr = [zeros(length(npts)) for i in 1:5]

    # Loop
    # for n in npts
    for I in eachindex(npts)
        n = npts[I]

        D1 = Grid2D(u->[u*0.5 - 0.25,-0.25],
            v->[-0.25,v*0.5 - 0.25],
            v->[0.25,v*0.5 - 0.25],
            u->[u*0.5 - 0.25, 0.25],
            n,n)

        T = FaADE.Grid.Torus([1.0],[1.0],[1],[0])

        # Right domain
        D2 = Grid2D(u->[0.25, -u*0.5 + 0.25], # Bottom
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Left
            v->v*(T(7π/4,0.0) + [-0.25, 0.25]) + [0.25, -0.25], # Right
            u->T(u*(7π/4 - 9π/4) + 9π/4,0.0), # Top
            n,n)

        # Top domain
        D3 = Grid2D(u->[u*0.5 - 0.25, 0.25], # Bottom
            v->v*(T(3π/4,0.0) + [0.25,-0.25]) + [-0.25,0.25], # Left
            v->v*(T(π/4,0.0) - [0.25,0.25]) + [0.25,0.25], # Right
            u->T(u*(π/4 - 3π/4) + 3π/4,0.0), # Top
            n,n)

        # Left domain
        D4 = Grid2D(u->[-0.25,u*0.5 - 0.25],
            v->v*(T(5π/4,0.0) - [-0.25, -0.25]) + [-0.25, -0.25],
            v->v*(T(3π/4,0.0) - [-0.25,0.25]) + [-0.25,0.25],
            u->T(u*(3π/4 - 5π/4) + 5π/4,0.0),
            n,n)

        # Bottom domain
        D5 = Grid2D(u->[-u*0.5 + 0.25, -0.25],
            v->v*(T(7π/4,0.0) - [0.25,-0.25]) + [0.25, -0.25],
            v->v*(T(5π/4,0.0) - [-0.25,-0.25]) + [-0.25, -0.25],
            u->T(u*(5π/4 - 7π/4) + 7π/4, 0.0),
            n,n)

        joints = ((Joint(2,Right),Joint(3,Up),Joint(4,Left),Joint(5,Down)),
            (Joint(1,Down),Joint(3,Left),Joint(5,Right)),
            (Joint(1,Down),Joint(4,Left),Joint(2,Right)),
            (Joint(1,Down),Joint(5,Left),Joint(3,Right)),
            (Joint(1,Down),Joint(2,Left),Joint(4,Right)))

        Dom = GridMultiBlock((D1,D2,D3,D4,D5),joints)

        if BType == Dirichlet
            Dr = FaADE.SATs.SAT_Dirichlet(Bxy, D2.Δy, Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
            Du = FaADE.SATs.SAT_Dirichlet(Bxy, D3.Δy, Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
            Dl = FaADE.SATs.SAT_Dirichlet(Bxy, D4.Δy, Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
            Dd = FaADE.SATs.SAT_Dirichlet(Bxy, D5.Δy, Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

            BD = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))
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

        u_MMS = generate_MMS(ũ,Dom,t_f)

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        # push!(relerr, norm(u_MMS[1] .- soln.u[2][1])/norm(u_MMS[1]) + norm(u_MMS[2] .- soln.u[2][2])/norm(u_MMS[2]))
        @show relerr[1][I] = norm(u_MMS[1] .- soln.u[2][1])/norm(u_MMS[1])
        @show relerr[2][I] = norm(u_MMS[2] .- soln.u[2][2])/norm(u_MMS[2])
        @show relerr[3][I] = norm(u_MMS[3] .- soln.u[2][3])/norm(u_MMS[3])
        @show relerr[4][I] = norm(u_MMS[4] .- soln.u[2][4])/norm(u_MMS[4])
        @show relerr[5][I] = norm(u_MMS[5] .- soln.u[2][5])/norm(u_MMS[5])
    end

    # conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))
    # conv_rate = log2.(relerr[1][1:end-1] ./ relerr[1][2:end], relerr[2][1:end-1] ./ relerr[2][2:end])
    conv_rate = [log2.(relerr[I][1:end-1] ./ relerr[I][2:end]) for I in 1:5]

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

    Bũ(X,t) = cos(2π*ωt*t) * sin(2π*ωx*X[1] + cx) * sin(2π*X[2]*ωy + cy) #Boundary condition x=0

    order = 2
    println("order=",order)
    O2_DirichletMMS = comp_MMS(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    order = 4
    println("order=",order)
    O4_DirichletMMS = comp_MMS(npts,
        Bũ,Dirichlet,
        FD,analytic,IC,order,
        kx=K,ky=K,θ=θ)

    println("Order 2 Dirichlet convergence rates=",O2_DirichletMMS.conv_rate)
    println("Order 4 Dirichlet convergence rates=",O4_DirichletMMS.conv_rate)

    println("Order 2 relative error=",O2_DirichletMMS.relerr)
    println("Order 4 relative error=",O4_DirichletMMS.relerr)

    println("=====")
end

for I in 1:5
    println("Order 2 convergence rate in block $I ",O2_DirichletMMS.conv_rate[I])
end
for I in 1:5 
    println("Order 2 relative error in block $I ",O2_DirichletMMS.relerr[I])
end
for I in 1:5
    println("Order 4 convergence rate in block $I ",O4_DirichletMMS.conv_rate[I])
end
for I in 1:5 
    println("Order 4 relative error in block $I ",O4_DirichletMMS.relerr[I])
end




if SaveTests
    using DelimitedFiles

    nameappend=string("spatial")

    open(string("testing/MMS/out/MMS_DoubleSeashell_O2",nameappend,".csv"),"w") do io
        writedlm(io,[npts O2_DirichletMMS.relerr])
    end

    open(string("testing/MMS/out/MMS_DoubleSeashell_O4",nameappend,".csv"),"w") do io
        writedlm(io,[npts O4_DirichletMMS.relerr])
    end
end



