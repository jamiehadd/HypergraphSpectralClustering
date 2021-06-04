module HypergraphDetectability

import Parameters
import Combinatorics
import Distributions
import SparseArrays
import LinearAlgebra
import DataFrames
import Arpack
import Kronecker


include("hypergraph.jl")
include("utils.jl")
include("affinity-functions.jl")
include("HSBM.jl")
include("samplers.jl")
include("matrices.jl")
include("experiments.jl")
include("eigenstuff.jl")
include("degrees.jl")


export hypergraph
export detectabilityData
export nonBacktrackingMatrix
export nonBacktrackingMatrices
export gapExperiment

export AffinityFunction
export partitionAffinityFunction
export sampleSBM

export aggregateEigenvector

export binaryClusters
export degreeTensor

export linearizedBPMatrix


end # module
