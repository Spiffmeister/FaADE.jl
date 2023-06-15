

cd("..")
using Interpolations
# push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./FaADE")
using FaADE







PB = FaADE.Helpers.ParallelPlane(-2π,[0.0 0.0; 0.5 π; 1.0 2π])

PF = FaADE.Helpers.ParallelPlane(2π,[0.0 0.0; 0.5 π; 1.0 2π])

P = FaADE.Helpers.ParallelGrid(PB,PF,0.0)
