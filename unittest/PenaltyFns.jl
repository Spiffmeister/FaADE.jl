

cd("..")
using Interpolations
# push!(LOAD_PATH,"./plas_diff")
push!(LOAD_PATH,"./SBP_operators")
using SBP_operators







PB = SBP_operators.Helpers.ParallelPlane(-2π,[0.0 0.0; 0.5 π; 1.0 2π])

PF = SBP_operators.Helpers.ParallelPlane(2π,[0.0 0.0; 0.5 π; 1.0 2π])

P = SBP_operators.Helpers.ParallelGrid(PB,PF,0.0)
