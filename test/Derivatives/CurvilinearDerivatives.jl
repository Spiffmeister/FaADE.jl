
using Test
using LinearAlgebra
using FaADE




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


@testset "2D Cartesian box" begin
    cbottom(u) = [u,0.0]
    cleft(v) = [0.0,v]
    cright(v) = [1.0,v]
    ctop(u) = [u,1.0]
    
    Dom = Grid2D(cbottom,cleft,cright,ctop,31,31)
    
    u₀(x,y) = x^2
    kx(x,y) = 1.0
    ky(x,y) = 1.0
    
    uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)->2.0)
    
    setKoefficient!(K,kx,ky,Dom)
    
    @testset "Second order" begin
        order = 2
        Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,order,false,:Variable)
        Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,order,false,:Variable)
    
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)
    
        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)
    
        @test norm(uxx[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1e-10
    end

    @testset "Fourth order" begin
        order = 4
        Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,order,false,:Variable)
        Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,order,false,:Variable)
    
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)
    
        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)
    
        @test norm(uxx[7:end-6,7:end-6] .- ue[7:end-6,7:end-6]) ≤ 1e-10
    end
end    




@testset "2D rescaled rectangle" begin
    # x∈[0,2], y∈[0,1] rescaled to computational domain q∈[0,1], r=y
    cbottom(u) = [2u,0.0]
    cleft(v) = [0.0,v]
    cright(v) = [2.0,v]
    ctop(u) = [2u,1.0]
    
    Dom = Grid2D(cbottom,cleft,cright,ctop,31,31)
    
    u₀(x,y) = x^2
    kx(x,y) = 1.0
    ky(x,y) = 1.0
    
    uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)->2.0)
    
    setKoefficient!(K,kx,ky,Dom)

    @testset "Second order" begin
        order = 2

        Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,order,false,:Variable)
        Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,order,false,:Variable)
    
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)
    
        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)
    
        uxx_real = uxx./Dom.J
    
        @test norm(uxx_real[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1e-10
    end

    @testset "Fourth order" begin
        order = 4

        Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,order,false,:Variable)
        Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,order,false,:Variable)
    
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)
    
        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)
    
        uxx_real = uxx./Dom.J
    
        @test norm(uxx_real[7:end-6,7:end-6] .- ue[7:end-6,7:end-6]) ≤ 1e-10
    end

end    




@testset "2D rhombus" begin
    # x∈[0,2], y∈[0,1] rescaled to computational domain q∈[0,1], r=y
    
    cbottom(u) = [-1.0,0.0] + u*[1,-1]
    cleft(v) = [-1.0,0.0] + v*[1,1]
    cright(v) = [0.0,-1.0] + v*[1,1]
    ctop(u) = [0.0,1.0] + u*[1,-1]
    
    Dom = Grid2D(cbottom,cleft,cright,ctop,31,31)

    kx(x,y) = 1.0
    ky(x,y) = 1.0

    Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,2,false,:Variable)
    Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,2,false,:Variable)

    @testset "Quadratic function" begin
        u₀(x,y) = x^2
        uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)->2.0)
        setKoefficient!(K,kx,ky,Dom)
        
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)

        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)

        uxx_real = uxx./Dom.J
    
        @test norm(uxx_real[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1e-10
    end

    @testset "Exponential function" begin
        u₀(x,y) = exp(-(x^2 + y^2)/0.02)
        uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)-> u₀(x,y) * 4.0*(-0.02 + x^2 + y^2)/0.02^2)
        setKoefficient!(K,kx,ky,Dom)
        
        D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)

        FaADE.Derivatives.mul!(uxx,u,K,D,0.0)

        uxx_real = uxx./Dom.J
    
        @test norm(uxx_real[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1e-10 broken = true
    end

end


@testset "chevron" begin
    println("chevron TEST NEEDS FIXING")
    function cbottom(u)
        x = u
        u ≤ 0.5 ? y = -u : y = u-1
        return [x,y]
    end

    cleft(v) = [0.0,v]

    function ctop(u)
        x = u
        u ≤ 0.5 ? y = 1-u : y = u
        return [x,y]
    end

    cright(v) = [1.0,v]

    kx(x,y) = 1.0
    ky(x,y) = 1.0

    Dom = Grid2D(cbottom,cleft,cright,ctop,31,31)
    
    Dx = FaADE.Derivatives.DiffusionOperator(Dom.nx,Dom.Δx,2,false,:Variable)
    Dy = FaADE.Derivatives.DiffusionOperator(Dom.ny,Dom.Δy,2,false,:Variable)

    u₀(x,y) = y^2

    uxx, u, K, ue = buildprob2D(Dom,u₀,kx,ky,(x,y)-> 2.0)
    setKoefficient!(K,kx,ky,Dom)

    D = FaADE.Derivatives.DiffusionOperatorND(Dx,Dy)

    FaADE.Derivatives.mul!(uxx,u,K,D,0.0)
    
    uxx_real = uxx./Dom.J

    @test norm(uxx_real[2:end-1,2:end-1] .- ue[2:end-1,2:end-1]) ≤ 1e-10 broken = true
end




