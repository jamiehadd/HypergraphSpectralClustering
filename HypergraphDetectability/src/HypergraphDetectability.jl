module HypergraphDetectability

using Parameters
using Combinatorics
using Distributions
using SparseArrays
using LinearAlgebra
using DataFrames
using Arpack
using Kronecker
using StatsBase

include("hypergraph.jl")
include("utils.jl")
include("affinity-functions.jl")
include("HSBM.jl")
include("samplers.jl")
include("matrices.jl")
include("experiments.jl")
include("eigenstuff.jl")
include("degrees.jl")
include("edge-incidence.jl")
include("data.jl")

export hypergraph
export detectabilityData
export nonBacktrackingMatrix
export nonBacktrackingMatrices
export nonBacktrackingMatrices_
export gapExperiment

export AffinityFunction
export partitionAffinityFunction
export sampleSBM

export aggregateEigenvector

export binaryClusters
export degreeTensor

export linearizedBPMatrix

export edgeIncidence
export pointedEdge
export reducedNonBacktrackingMatrices

export reducedNonBacktrackingMatrix
export reducedBPJacobian
end # module
