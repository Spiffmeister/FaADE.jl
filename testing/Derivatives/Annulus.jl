using Revise
using FaADE
using LinearAlgebra
# using GLMakie

plot = false

function buildprob2D(D,f,cx,cy,e)
    ue  = zeros(eltype(D.gridx),size(D))
    u   = zeros(eltype(D.gridx),size(D))
    kx  = zeros(eltype(D.gridx),size(D))
    ky  = zeros(eltype(D.gridx),size(D))
    Dxky = zeros(eltype(D.gridx),size(D))
    Dykx = zeros(eltype(D.gridx),size(D))
    for i in 1:D.nx
        for j in 1:D.ny
            if FaADE.Grid.coordtype(D) == CartesianMetric
                ue[i,j] = e(D.gridx[i,j],D.gridy[i,j])
                u[i,j]  = f(D.gridx[i,j],D.gridy[j])
            else
                ue[i,j] = e(D.gridx[i,j],D.gridy[i,j])
                u[i,j]  = f(D.gridx[i,j],D.gridy[i,j])
            end
        end
    end
    uxx = zeros(eltype(D.gridx),(D.nx,D.ny))
    return uxx, u, [kx, ky, Dxky, Dykx], ue
end

function setKoefficient!(K,kx,ky,G::FaADE.Grid.LocalGridType{TT,2,MET}) where {TT,MET}
    if MET == CartesianMetric
        for i in 1:G.nx
            for j in 1:G.ny
                K[1][i,j] = kx(G[i,j]...)                  # Kx
                K[2][i,j] = ky(G[i,j]...)                  # Ky
            end
        end
    else
        for i in 1:G.nx
            for j in 1:G.ny
                K[1][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]^2 + G.qy[i,j]^2)                  # Kx
                K[2][i,j] = ky(G[i,j]...) * G.J[i,j] * (G.rx[i,j]^2 + G.ry[i,j]^2)                  # Ky
                K[3][i,j] = kx(G[i,j]...) * G.J[i,j] * (G.qx[i,j]*G.rx[i,j] + G.qy[i,j]*G.ry[i,j])  # DxKy
                K[4][i,j] = ky(G[i,j]...) * G.J[i,j] * (G.qx[i,j]*G.rx[i,j] + G.qy[i,j]*G.ry[i,j])  # DyKx
            end
        end
    end
    K
end


Rin = [3e-1]; Zin=[3e-1]
Rout = [6e-1]; Zout=[6e-1]

inner = FaADE.Grid.Torus(Rin,Zin,[1],[0])
outer = FaADE.Grid.Torus(Rout,Zout,[1],[0])

nx = 101
ny = 201

X,Y = FaADE.Grid.meshgrid(inner,outer,0.0,nx,ny)


# BUILD GRID
Dom = Grid2D(X,Y,periodicy=true)


# u₀(x,y) = sqrt(x^2 + y^2)
# lapu(x,y) = 1/sqrt(x^2 + y^2)

u₀(x,y) = sin(2π*x)*cos(2π*y)
lapu(x,y) = -8*π^2*sin(2π*x)*cos(2π*y)


e = zeros(eltype(Dom),size(Dom));
for I in eachindex(Dom)
    e[I] = lapu(Dom[I]...);
end


Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,4,false,:Variable)
Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,4,true,:Variable)

D = FaADE.Derivatives.DiffusionOperatorND((Dx,Dy))

kx(x,y) = 1.0
ky(x,y) = 1.0

uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)-> 0.0)
setKoefficient!(K,(x,y) -> 1.0, (x,y) -> 1.0,Dom)



FaADE.Derivatives.mul!(uxx,u,K,D,0.0)

# || (u/J - e)√(JrΔxΔy) ||

println( norm((uxx./Dom.J .- e).*sqrt.(sqrt.(Dom.gridx.^2 + Dom.gridy.^2).*Dom.J*Dom.Δx*Dom.Δy) ) )



if plot
    using GLMakie

    f = GLMakie.Figure()
    ax1 = GLMakie.Axis3(f[1,1])#,aspect=DataAspect())
    ax2 = GLMakie.Axis3(f[1,2])#,aspect=DataAspect())
    GLMakie.surface!(ax1,Dom.gridx,Dom.gridy,uxx./Dom.J)
    GLMakie.surface!(ax2,Dom.gridx,Dom.gridy,e)
end
