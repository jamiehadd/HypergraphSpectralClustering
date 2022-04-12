function gapExperiment(n, c₂, c₃; grid₂ = 0.0:0.1:1.0, grid₃ = 0:0.1:1.0)
    DF = DataFrames.DataFrame()
    for p₂ ∈ grid₂, p₃ ∈ grid₃
        H = detectabilityData(n, c₂, c₃, p₂, p₃)

        try
            B_ = HypergraphNB.reducedNonBacktrackingMatrix(H);
            E = Arpack.eigs(B_; nev = 3);
        
            gap = abs(E[1][2]) - abs(E[1][3])

            df = DataFrames.DataFrame(
                p₂ = p₂, 
                p₃ = p₃, 
                gap = gap
            )
            DF = vcat(DF, df)
        catch e
            nothing
        end
    end
    return DF
end