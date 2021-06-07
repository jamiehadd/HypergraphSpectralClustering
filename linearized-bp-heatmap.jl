using HypergraphDetectability
using DataFrames
using Clustering
using Statistics
using RCall

n  = 100
c₂ = 5
c₃ = 5
z = 1 .+ (1:n .> n/2);

DF = DataFrame()

"""
this experiment cheats a little bit by allowing us to estimate the C matrices by using the true partition in linearizedBPMatrix, in future experiments we should alternate estimation steps. 
"""

"""
one of the main learnings of this experiment is that the eigenvector of the linearized BP matrix that we care about might be the FIRST eigenvector. This phenomenon needs to be understood, and the 
"""

for p₂ ∈ 0:0.05:1, p₃ ∈ 0.0:0.05:1, ev ∈ 1:2, rep ∈ 1:5
    H = detectabilityData(n, c₂, c₃, p₂, p₃);
    BP_mat, ix = linearizedBPMatrix(H, z);
    
    # E = Arpack.eigs(BP_mat; nev = 50, ritzvec = true)

    # plot(E[1], seriestype = :scatter)

    # v = E[2][:,1]

    # # plot(real.(v))

    # u = aggregateEigenvector(v, ix)
    # plot(real.(u))

    try 
        E = Arpack.eigs(BP_mat; nev = 2, ritzvec = true)
        v = E[2][:,ev]
        u = aggregateEigenvector(v, ix)
        if !(mean(abs.(imag.(u))) ≈ 0)
            MI = 0.0
            df = DataFrame(p_2 = p₂, p_3 = p₃, MI = MI, ev = ev, rep = rep)
            append!(DF, df)
        else
            clusters = 1 .+ (real.(u) .> 0)
            MI = mutualinfo(clusters, z)
            df = DataFrame(p_2 = p₂, p_3 = p₃, MI = MI, ev = ev, rep = rep)
            append!(DF, df)
        end
    catch e
        MI = NaN
        df = DataFrame(p_2 = p₂, p_3 = p₃, MI = MI, ev = ev, rep = rep)
        append!(DF, df)
    end
end



R"""
library(tidyverse)
library(viridis)

df <- tibble($(DF))

df <- df %>%
    group_by(p_2, p_3, ev) %>% 
    summarise(MI = mean(MI, na.rm = T))

p <- df %>% 
    ggplot() + 
    aes(x = p_2, y = p_3, fill = MI) + 
    geom_tile() + 
    theme_minimal() + 
    scale_fill_viridis(option = "inferno", limits = c(0, 1)) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    xlab(expression(p[2])) + 
    ylab(expression(p[3])) + 
    theme(panel.grid = element_blank()) + 
    coord_fixed() + 
    facet_wrap(~ev)

ggsave("fig/linearized-bp-heatmap.png", p,  width = 9, height = 4)
"""




