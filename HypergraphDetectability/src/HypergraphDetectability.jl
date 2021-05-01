module HypergraphDetectability

import Parameters
import Combinatorics
import Distributions
import SparseArrays
import LinearAlgebra
import DataFrames
import Arpack

include("hypergraph.jl")
include("utils.jl")
include("affinity-functions.jl")
include("HSBM.jl")
include("samplers.jl")
include("matrices.jl")
include("experiments.jl")



export hypergraph
export detectabilityData
export nonBacktrackingMatrix
export gapExperiment

export AffinityFunction
export partitionAffinityFunction
export sampleSBM

end # module
