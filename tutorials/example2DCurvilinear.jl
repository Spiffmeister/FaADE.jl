# # 2D on a circular domain with parallel map
#
# For this example we will solve the head equation with a mangetic field aligned with the grid in a circular domain
#
# ```math
#   \mathbf{B} = (0,0,1)
# ```
# In this case we expect the parallel operator to do nothing since ``\mathbf{P}_f=\mathbf{P}_b=I``
#
using Revise
using FaADE

# 
# We first need to create a domain to solve the PDE using 
#

nx = ny = 21;
order = 2;

# The interior grid will be a square, 

g0lw = 0.25

D1 = Grid2D(
    u -> [-g0lw, -g0lw] + u * ([g0lw, -g0lw] - [-g0lw, -g0lw]),
    v -> [-g0lw, -g0lw] + v * ([-g0lw, g0lw] - [-g0lw, -g0lw]),
    v -> [g0lw, -g0lw] + v * ([g0lw, g0lw] - [g0lw, -g0lw]),
    u -> [-g0lw, g0lw] + u * ([g0lw, g0lw] - [-g0lw, g0lw]),
    nx, ny)

# For the surrounding domains we need the shape of the circular domain. There is an inbuilt function for generating toroidal surfaces which we can use to construct the boundary curves,

R = 1.0

Tor = FaADE.Grid.Torus([R],[R],[1],[0])


# Right domain
D2 = Grid2D(#u->[0.25, -u*0.5 + 0.25], # Bottom
    u -> [g0lw, g0lw] + u * ([g0lw, -g0lw] - [g0lw, g0lw]),
    v -> v * (Tor(π / 4, 0.0) - [g0lw, g0lw]) + [g0lw, g0lw], # Left
    v -> v * (Tor(7π / 4, 0.0) + [-g0lw, g0lw]) + [g0lw, -g0lw], # Right
    u -> Tor(u * (7π / 4 - 9π / 4) + 9π / 4, 0.0), # Top
    nx, ny)

# Top domain
D3 = Grid2D(#u->[u*0.5 - g0lw, g0lw], # Bottom
    u -> [-g0lw, g0lw] + u * ([g0lw, g0lw] - [-g0lw, g0lw]),
    v -> v * (Tor(3π / 4, 0.0) + [g0lw, -g0lw]) + [-g0lw, g0lw], # Left
    v -> v * (Tor(π / 4, 0.0) - [g0lw, g0lw]) + [g0lw, g0lw], # Right
    u -> Tor(u * (π / 4 - 3π / 4) + 3π / 4, 0.0), # Top
    nx, ny)

# Left domain
D4 = Grid2D(#u->[-g0lw,u*0.5 - g0lw],
    u -> [-g0lw, -g0lw] + u * ([-g0lw, g0lw] - [-g0lw, -g0lw]),
    v -> v * (Tor(5π / 4, 0.0) - [-g0lw, -g0lw]) + [-g0lw, -g0lw],
    v -> v * (Tor(3π / 4, 0.0) - [-g0lw, g0lw]) + [-g0lw, g0lw],
    u -> Tor(u * (3π / 4 - 5π / 4) + 5π / 4, 0.0),
    nx, ny)


# Bottom domain
D5 = Grid2D(#u->[-u*0.5 + g0lw, -g0lw],
    u -> [-g0lw, -g0lw] + u * ([g0lw, -g0lw] - [-g0lw, -g0lw]),
    v -> v * (Tor(7π / 4, 0.0) - [g0lw, -g0lw]) + [g0lw, -g0lw],
    v -> v * (Tor(5π / 4, 0.0) - [-g0lw, -g0lw]) + [-g0lw, -g0lw],
    u -> Tor(u * (5π / 4 - 7π / 4) + 7π / 4, 0.0),
    nx, ny)



# Finally, we need to tell the grid how everything is connected, this can be done by setting a tuple of boundary objects,


joints = ((Joint(2, Right), Joint(3, Up), Joint(4, Left), Joint(5, Down)), #Domain 1
    (Joint(1, Down), Joint(3, Left), Joint(5, Right)), #Domain 2
    (Joint(1, Down), Joint(4, Left), Joint(2, Right)), #Domain 3
    (Joint(1, Down), Joint(5, Left), Joint(3, Right)), #Domain 4
    (Joint(1, Down), Joint(2, Left), Joint(4, Right))) #Domain 5

# We can then build the grid,

Dom = GridMultiBlock((D1, D2, D3, D4, D5), joints)



# To plot the domain,

using Plots

gridfig = plot()
for grid in Dom.Grids
    scatter!(grid.gridx[:],grid.gridy[:],markersize=2)
end
gridfig


# Before continuing we'll set up the boundary conditions,

Bxy(X, t) = 0.0

Dr = SAT_Dirichlet(Bxy, D2.Δy, FaADE.Up, order, D2.Δx, :Curvilinear) # Block 2 BCs
Du = SAT_Dirichlet(Bxy, D3.Δy, FaADE.Up, order, D3.Δx, :Curvilinear) # Block 3 BCs
Dl = SAT_Dirichlet(Bxy, D4.Δy, FaADE.Up, order, D4.Δx, :Curvilinear) # Block 4 BCs
Dd = SAT_Dirichlet(Bxy, D5.Δy, FaADE.Up, order, D5.Δx, :Curvilinear) # Block 5 BCs

# For multiblock problems we store the boundaries in a dictionary, where the key is the index for the grid,

BCs = Dict(2 => (Dr,), 3 => (Du,), 4 => (Dl,), 5 => (Dd,))



# Now for setting up the parallel map, first we specifiy the magnetic field,

function Bfield(X, x, p, t)
    r = x[1]^2 + x[2]^2

    X[1] = 2x[2] * exp(1 - r)
    X[2] = -2x[1] * exp(1 - r)

    if (x[1] == 0.0) && (x[2] == 0.0)
        X[1] = 0.0
        X[2] = 0.0
    end
end


# this needs to be written in the format that Julias [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) takes,

dH(X, x, params, t) = Bfield(X, x, params, t)

# The field lines are defined in ``(s,\theta)`` space, but our grid is in Cartesian space. Therefore we need to define the coordinate maps which take us from the grid to the field coordinates and then back

XtoB(x, y) = [sqrt(x^2 + y^2), atan(y, x)]
BtoX(r, θ) = [r * cos(θ), r * sin(θ)]

# now we can construct the parallel grids,

gridoptions = Dict("xbound" => [-10.0, 10.0], "ybound" => [-10.0, 10.0])
gdata = construct_grid(dH, Dom, [-1.0, 1.0])




# For the parallel map, we'll specify one last function called `intercept`, which will give the field lines a value if they leave the domain just in case,

intercept(u, x, y, t) = begin
    if (sqrt(x^2 + y^2) ≈ 1.0)
        tmp = zero(typeof(x))
        return tmp
    else
        return u
    end
end


# then we create a `Dict` of the options for the parallel map which tell the constructor to use `CubicHermiteSpline`s for the interpolant and specify the intercept

interpotions = Dict("interpolant" => :chs, "intercept" => intercept)

# now we construct the map

PData = ParallelMultiBlock(gdata, Dom, order, κ=1.0, interpopts=interpotions)


# Since we're in a circular domain with ``u=0`` boundary conditions, we need a source term for something to happen, we'll add a circular source term

function F(X,t)
    x,y = X
    tmp = (2*(2x^2 - 1) + 2*(2y^2 - 1)) * exp(1 - (x^2 + y^2)) * (1 - exp(-2π^2*t))
    
    tmp2 = -2π^2 * (exp(1 - (x^2 + y^2)) - 1) * exp(-2π^2*t)
    
    return -(tmp + tmp2)
end

# Finally we can set the initial condition, time step and construct the `Problem2D` object

κ_perp = 1.0

u₀(x, y) = 0.0



Δt = 1.0e-4;
t_f = 100Δt;


P = Problem2D(order,u₀,κ_perp,κ_perp,Dom,BCs,source=F,parallel=PData)

# now we can call the solver,

soln = solve(P,Dom,Δt,t_f)


# Finally we can plot the solution, since it's symmetric about the origin we'll just plot a contour across the centre line,

fig = plot()
plot!(Dom.Grids[1].gridx[:,11], soln.u[2][1][:,11])
plot!(Dom.Grids[2].gridx[11,:], soln.u[2][2][11,:])
plot!(Dom.Grids[4].gridx[11,:], soln.u[2][4][11,:])
fig

