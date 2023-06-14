

cd("..")
using Interpolations
# push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./SPADE")
using SPADE







PB = SPADE.Helpers.ParallelPlane(-2π,[0.0 0.0; 0.5 π; 1.0 2π])

PF = SPADE.Helpers.ParallelPlane(2π,[0.0 0.0; 0.5 π; 1.0 2π])

P = SPADE.Helpers.ParallelGrid(PB,PF,0.0)
