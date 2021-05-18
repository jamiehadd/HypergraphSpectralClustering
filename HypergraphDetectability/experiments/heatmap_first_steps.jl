##
using Revise
using HypergraphDetectability
using SparseArrays
using LinearAlgebra
using Arpack
using DataFrames
using VegaLite
##


## 
DF = gapExperiment(500, 5, 5; grid₂ = 0:0.05:1.0, grid₃ = 0:0.05:1.0)
## 

##
DF |> @vlplot(:rect, y = "p₃:o", x="p₂:o", color="mean(gap)")
##