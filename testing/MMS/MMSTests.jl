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
npts = [21,31,41,51,61]

###=== MMS ===###



# Solution
ũ(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*y*ωy + cy)

# Initial condition
ũ₀(x,y;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = sin(2π*ωx*x + cx) * sin(2π*ωy*y + cy)
    
# Dirichlet boundaries
Dx0_Lũ(y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) =        cos(2π*t) * sin(cx) * sin(2π*y*ωy + cy)
DxL_Rũ(y,t;
    ωx=1.0,cx=0.0,Lx=1.0,
    ωy=1.0,cy=0.0) =        cos(2π*t) * sin(2π*Lx*ωx + cx) * sin(2π*y*ωy + cy)

Dy0_Lũ(x,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) =        cos(2π*t) * sin(2π*x*ωx + cx) * sin(cy)
DyL_Rũ(x,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0,Ly=1.0) = cos(2π*t) * sin(2π*x*ωx + cx) * sin(2π*Ly*ωy + cy)


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




# Source term based on MMS
#F = ∂ₜũ - K∇ũ
K = 1.0
F(x,y,t;
    ωx=1.0,cx=0.0,
    ωy=1.0,cy=0.0) = 
        -2π*sin(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) + 
        K * 4π^2 * (ωx^2 + ωy^2) * cos(2π*t)*sin(2π*x*ωx + cx)*sin(2π*y*ωy + cy) 
            





# Generates the exact MMS solution
function generate_MMS(MMS::Function,grid::SBP_operators.Helpers.Grid2D,t::Float64;ωx=0.0,cx=0.0,ωy=0.0,cy=0.0)
    u_MMS = zeros(grid.nx,grid.ny)
    for j = 1:grid.ny
        for i = 1:grid.nx
            u_MMS[i,j] = MMS(grid.gridx[i],grid.gridy[j],t,ωx=ωx,cx=cx,ωy=ωy,cy=cy)
        end
    end
    return u_MMS
end





function comp_MMS(Dx,Dy,npts,
        BoundaryX0,BX0Type,BoundaryXL,BXLType,
        BoundaryY0,BY0Type,BoundaryYL,BYLType,
        F,ũ,ũ₀,order;
        dt_scale=0.01,t_f=0.1,kx=1.0,ky=1.0,
        ωx=1.0,ωy=1.0,cx=0.0,cy=0.0,Lx=1.0,Ly=1.0)

    comp_soln = []
    MMS_soln = []
    grids = []
    relerr = []
    # X boundaries
    if BX0Type != Periodic
        Bx0 = Boundary(BX0Type,(y,t)->BoundaryX0(y,t,ωx=ωx,cx=cx,ωy=ωy,cy=cy),Left,1)
        BxL = Boundary(BXLType,(y,t)->BoundaryXL(y,t,ωx=ωx,cx=cx,ωy=ωy,cy=cy,Lx=Lx),Right,1)
    else
        Bx0L = PeriodicBoundary(1)
    end
    # Y boundaries
    if BY0Type != Periodic
        By0 = Boundary(BY0Type,(x,t)->BoundaryY0(x,t,ωx=ωx,cx=cx,ωy=ωy,cy=cy),Up,2)
        ByL = Boundary(BYLType,(x,t)->BoundaryYL(x,t,ωx=ωx,cx=cx,ωy=ωy,cy=cy,Ly=Ly),Down,2)
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
        soln = solve(P,Dom,Δt,t_f,:cgie,source=(x,y,t)->F(x,y,t,ωx=ωx,cx=cx,ωy=ωy,cy=cy))

        u_MMS = generate_MMS(ũ,Dom,t_f,ωx=ωx,cx=cx,ωy=ωy,cy=cy)

        push!(comp_soln,soln)
        push!(grids,Dom)
        push!(MMS_soln,u_MMS)
        push!(relerr, norm(u_MMS .- soln.u[2])/norm(MMS_soln))
    end

    conv_rate = log.(relerr[1:end-1]./relerr[2:end]) ./ log.( (1 ./ (npts[1:end-1].-1))./(1 ./ (npts[2:end].-1) ))

    return (comp_soln=comp_soln,MMS_soln=MMS_soln,grids=grids,relerr=relerr,conv_rate=conv_rate,npts=npts)
end




###=== MMS TESTS ===###

order = 2
println("=====")
println("Order=",order)
println("=====")
# Dirichlet
println("Dirichlet")
O2_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    Dx0_Lũ,Dirichlet,DxL_Rũ,Dirichlet,
    Dy0_Lũ,Dirichlet,DyL_Rũ,Dirichlet,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)
#=
# Neumann
println("Neumann")
O2_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    Nx0_Lũ,Neumann,NxL_Rũ,Neumann,
    Ny0_Lũ,Neumann,NyL_Rũ,Neumann,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)

# Dirichlet x Neumann y
println("Dirichlet x Neumann y")
O2_DirichXNeuYMMS = comp_MMS(𝒟x,𝒟y,npts,
    Dx0_Lũ,Dirichlet,DxL_Rũ,Dirichlet,
    Ny0_Lũ,Neumann,NyL_Rũ,Neumann,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)

# Periodic
println("Periodic")
O2_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    nothing,Periodic,nothing,Periodic,
    nothing,Periodic,nothing,Periodic,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)
=#

### Order 4
order = 4
println("=====")
println("Order=",order)
println("=====")
println("Dirichlet")
O4_DirichletMMS = comp_MMS(𝒟x,𝒟y,npts,
    Dx0_Lũ,Dirichlet,DxL_Rũ,Dirichlet,
    Dy0_Lũ,Dirichlet,DyL_Rũ,Dirichlet,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)
#=
# Neumann
println("Neumann")
O4_NeumannMMS = comp_MMS(𝒟x,𝒟y,npts,
    Nx0_Lũ,Neumann,NxL_Rũ,Neumann,
    Ny0_Lũ,Neumann,NyL_Rũ,Neumann,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)

# Dirichlet x Neumann y
println("Dirichlet x Neumann y")
O4_DirichXNeuYMMS = comp_MMS(𝒟x,𝒟y,npts,
    Dx0_Lũ,Dirichlet,DxL_Rũ,Dirichlet,
    Ny0_Lũ,Neumann,NyL_Rũ,Neumann,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)

# Periodic
println("Periodic")
O4_PeriodicMMS = comp_MMS(𝒟x,𝒟y,npts,
    nothing,Periodic,nothing,Periodic,
    nothing,Periodic,nothing,Periodic,
    F,ũ,ũ₀,order,
    ωx=1.0,cx=0.0,
    ωy=2.5,cy=1.0)
=#



###=== PLOTTING ===###
println("Order 2 convergence rates=",O2_DirichletMMS.conv_rate)#,O2_NeumannMMS.conv_rate,O2_DirichXNeuYMMS.conv_rate,O2_PeriodicMMS.conv_rate)
println("Order 4 convergence rates=",O4_DirichletMMS.conv_rate)#,O4_NeumannMMS.conv_rate,O4_DirichXNeuYMMS.conv_rate,O4_PeriodicMMS.conv_rate)


###=== PLOTTING ===###

p = plot(   log.(O2_DirichletMMS.npts),     O2_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^2)$")
plot!(p,    log.(O4_DirichletMMS.npts),     O4_DirichletMMS.relerr,     label=L"Dirichlet $\mathcal{O}(h^4)$")
#=
plot!(p,    log.(O2_NeumannMMS.npts),       O2_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^2)$")
plot!(p,    log.(O2_DirichXNeuYMMS.npts),   O2_DirichXNeuYMMS.relerr,   label=L"Dirichlet and Neumann $\mathcal{O}(h^2)$")
plot!(p,    log.(O2_PeriodicMMS.npts),      O2_PeriodicMMS.relerr,      label=L"Periodic $\mathcal{O}(h^2)$")

plot!(p,    log.(O4_NeumannMMS.npts),       O4_NeumannMMS.relerr,       label=L"Neumann $\mathcal{O}(h^4)$")
plot!(p,    log.(O4_DirichXNeuYMMS.npts),   O4_DirichXNeuYMMS.relerr,   label=L"Dirichlet and Neumann $\mathcal{O}(h^4)$")
plot!(p,    log.(O4_PeriodicMMS.npts),      O4_PeriodicMMS.relerr,      label=L"Periodic $\mathcal{O}(h^4)$")


savefig(p,".//testing//MMS//MMSTests.eps")
=#