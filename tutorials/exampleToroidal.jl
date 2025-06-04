# # 3D Example
#
# Constructing a toroidal grid.
#

using FaADE

# We need to specify the torus parameters, for this example we will use a simple torus with major radius 1 and minor radius 0.5
#
# ```math
# R(\theta,\zeta) = \sum R_{i,j} \cos(m_j\theta - n_i\zeta)
# Z(\theta,\zeta) = \sum Z_{i,j} \sin(m_j\theta - n_i\zeta)
# ```

Rin = [5e-1]; Zin=[5e-1]
inn = [1]; inm = [1]

Rout = [1.0]; Zout=[1.0]
outn = [1]; outm = [1]

# The torus is then constructed using the `Torus` function

T0 = Torus(Rin,Zin,inn,inm)
T1 = Torus(Rout,Zout,outn,outm)

# Providing the two tori to the `meshgrid` function will create a 2D grid in the toroidal coordinates at a given $\zeta$ plane and with a given number of radial points `nr` and poloidal points `nθ`
# ```julia
# meshgrid(inner,outer,ζ,nr,nθ)
# ```

X,Y = FaADE.Grid.meshgrid(T0,T1,0.0,11,21)

# which can then be used to create a `Grid2D` object

Dom = Grid2D(X,Y)

