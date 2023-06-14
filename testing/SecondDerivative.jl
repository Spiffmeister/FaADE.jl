# == TESTING FIRST DERIVATIVE ==#

using Test
using LinearAlgebra

using Pkg; Pkg.activate("."); using SBP_operators


function buildgrid1d(n)
    x = collect(range(0,stop=1,length=n+1))
    Δx = x[2]-x[1]
    return n+1, x, Δx
end

function buildgrid2d(n)
    x = collect(range(0,stop=1,length=n+1))
    y = collect(range(0,stop=1,length=n+1))
    Δx = x[2]-x[1]
    Δy = y[2]-y[1]
    return n+1, x,y, Δx,Δy
end





#=

##======##
# FOURTH ORDER
##======##

### Constant coefficient
# Linear Function
n, x, Δx = buildgrid1d(15)
c = ones(n);
u = x;
∂ₓₓuₑ = zeros(n); ∂ₓₓu = zeros(n);
SBP_operators.Derivatives.Dₓₓ!(∂ₓₓu,u,c,n,Δx,4)

@test all(∂ₓₓuₑ .- ∂ₓₓu .≤ 1e-10)



# Quadratic function
n, x, Δx = buildgrid1d(15)
c = ones(n);
u = x.^2;

∂ₓₓuₑ = zeros(n);
∂ₓₓuₑ .= 2.0;

∂ₓₓu = zeros(n);
SBP_operators.Derivatives.Dₓₓ!(∂ₓₓu,u,c,n,Δx,4)

@test all(∂ₓₓuₑ .- ∂ₓₓu .≤ 1e-10)



# Cubic function
n, x, Δx = buildgrid1d(15)
c = ones(n);
u = x.^3;

∂ₓₓuₑ = zeros(n); 
∂ₓₓuₑ = 6x

∂ₓₓu = zeros(n);
SBP_operators.Derivatives.Dₓₓ!(∂ₓₓu,u,c,n,Δx,4)

@test all(∂ₓₓuₑ .- ∂ₓₓu .≤ 1e-10)





#### Linear Function
n, x,y,Δx,Δy = buildgrid2d(15)
cx = cy = ones(n,n);

u = zeros(n,n); ∂∂ue = zeros(n,n);
∂∂u2 = zeros(n,n); ∂∂u4 = zeros(n,n);

u₀(x,y) = x
for i = 1:n, j = 1:n
    u[i,j] = u₀(x[i],y[j])
end

SBP_operators.Derivatives.Dₓₓ!(∂∂u2,u,cx,cy,n,n,Δx,Δy,2,2)
SBP_operators.Derivatives.Dₓₓ!(∂∂u4,u,cx,cy,n,n,Δx,Δy,4,4)

norm(∂∂ue .- ∂∂u2)
norm(∂∂ue .- ∂∂u4)




#### Quadratic Function
n, x,y,Δx,Δy = buildgrid2d(15)
cx = cy = ones(n,n);

u = zeros(n,n);
∂∂u2 = zeros(n,n); ∂∂u4 = zeros(n,n);

∂∂ue = zeros(n,n);
for i = 1:n, j = 1:n
    ∂∂ue[i,j] = 2.0
end

u₀(x,y) = x.^2
for i = 1:n, j = 1:n
    u[i,j] = u₀(x[i],y[j])
end

SBP_operators.Derivatives.Dₓₓ!(∂∂u2,u,cx,cy,n,n,Δx,Δy,2,2)
SBP_operators.Derivatives.Dₓₓ!(∂∂u4,u,cx,cy,n,n,Δx,Δy,4,4)

norm(∂∂ue[2:end-1,:] .- ∂∂u2[2:end-1,:])
norm(∂∂ue .- ∂∂u4)



#### Cubic Function
n, x,y,Δx,Δy = buildgrid2d(15)
cx = cy = ones(n,n);

u = zeros(n,n);
∂∂u2 = zeros(n,n); ∂∂u4 = zeros(n,n);
∂∂ue = zeros(n,n);

u₀(x,y) = x.^3
for i = 1:n, j = 1:n
    u[i,j] = u₀(x[i],y[j])
end
for i = 1:n, j = 1:n
    ∂∂ue[i,j] = 6*x[i]
end

SBP_operators.Derivatives.Dₓₓ!(∂∂u2,u,cx,cy,n,n,Δx,Δy,2,2)
SBP_operators.Derivatives.Dₓₓ!(∂∂u4,u,cx,cy,n,n,Δx,Δy,4,4)

norm(∂∂ue[2:end-1,:] .- ∂∂u2[2:end-1,:])
norm(∂∂ue .- ∂∂u4)

=#





#=
n, x,y,Δx,Δy = buildgrid1d(40)
c = ones(n);

u = zeros(n);
∂∂u2 = zeros(n); ∂∂u4 = zeros(n);
∂∂ue = zeros(n);

ωx = 1.0; kx = 0.0;
ωy = 1.0; ky = 0.0;
for i = 1:n, j = 1:n
    ∂∂ue[i] = -4π^2 * ωx^2 * sin(2π*x[i]*ωx + kx)
end

u₀(x) = sin(2π*x*ωx+kx)
for i = 1:n, j = 1:n
    u[i] = u₀(x[i])
end

SBP_operators.Derivatives.Dₓₓ!(∂∂u2,u,c,n,Δx,2)
SBP_operators.Derivatives.Dₓₓ!(∂∂u4,u,c,n,Δx,4)

norm(∂∂ue .- ∂∂u2)
norm(∂∂ue .- ∂∂u4)


=#



################

n, x,y,Δx,Δy = buildgrid2d(40)
cx = cy = ones(n,n);

u = zeros(n,n);
∂∂u2 = zeros(n,n); ∂∂u4 = zeros(n,n);
∂∂ue = zeros(n,n);

ωx = 1.0; kx = 0.0;
ωy = 1.0; ky = 0.0;
for i = 1:n, j = 1:n
    ∂∂ue[i,j] = -4π^2 * (ωx^2 + ωy^2) * sin(2π*x[i]*ωx + kx)*sin(2π*y[j]*ωy + ky)
end

u₀(x,y) = sin(2π*x*ωx+kx)*sin(2π*y*ωy+ky)
for i = 1:n, j = 1:n
    u[i,j] = u₀(x[i],y[j])
end

SBP_operators.Derivatives.Dₓₓ!(∂∂u2,u,cx,cy,n,n,Δx,Δy,2,2)
SBP_operators.Derivatives.Dₓₓ!(∂∂u4,u,cx,cy,n,n,Δx,Δy,4,4)


norm(∂∂ue .- ∂∂u2)
norm(∂∂ue .- ∂∂u4)




