using Revise
using HypergraphNB

using DataFrames
using Clustering
using Statistics
using RCall
using Arpack
using MultivariateStats
using TSne
using Random
Random.seed!(1234);

# make some fake data to play with 

n_groups = 3
n_nodes = 50
N = repeat([n_nodes], n_groups);

# NaN refers to 1 edges, which don't matter
P = [NaN, 0.9, 0.1];
C = [NaN, 5, 5];
H = plantedPartitionHypergraph(N, C, P);
z = vcat([repeat([z], N[z]) for z ∈ 1:length(N)]...); # true labels

B = reducedBPJacobian(H, z);
E = Arpack.eigs(B; nev = 5000);


λ₁ = maximum(real.(E[1]))
E[1][(imag.(E[1]) .== 0 ) .& (abs.(real.(E[1])) .> sqrt(λ₁))]

R"""

library(tidyverse)
library(ggforce)
df <- tibble(R = $(real.(E[1])), C = $(imag.(E[1])))

real_df <- df %>% 
    filter(C == 0, abs(R) > 1)

other_df <- df %>% 
    filter((C =! 0) | (abs(R) <= 1))

highlight_df <- real_df %>% 
    arrange(desc(abs(R))) %>% 
    filter(row_number() <= 4)

rad <- sqrt($(λ₁))

p <- other_df %>% 
    ggplot() + 
    geom_point(aes(x = R, y = C), size = .2, pch = 21, alpha = .3) + 
    geom_point(aes(x = R, y = C), data = real_df, fill = "firebrick", pch = 21) + 
#    geom_point(aes(x = R, y = C), data = highlight_df, pch = 21, size = 4) + 
    theme_minimal() +
    xlab(expression(Re~lambda)) + 
    ylab(expression(Im~lambda)) + 
    # geom_circle(aes(x0 = 0, y0 = 0, r = rad), color = "grey", linetype = "dashed") + 
    coord_fixed() + 
    ggtitle("(a). Spectrum of BP Jacobian")

ggsave("fig/spectrum.png", p, width = 5, height = 4)

p    
"""

λ₁ = maximum(real.(E[1]))
E[1][(imag.(E[1]) .== 0 ) .& (abs.(real.(E[1])) .> sqrt(λ₁))]

E_ = E[2][:, (imag.(E[1]) .== 0 ) .& (abs.(real.(E[1])) .> 1)];


V = hcat([HypergraphNB.transform_eigenvector(real.(E_[:,i]), H) for i ∈ 1:size(E_)[2]]...)
# V = 1.0*(V .> 0)   
M = fit(PCA, V; maxoutdim=size(E_)[2])


M = tsne(V, 2, 100, 1000, 10.0);

# M = fit(PCA, V; maxoutdim=2)

R"""
library(colorspace)
df <- tibble(z = as.factor($z), m1 = $(M[:,1]), m2 = $(M[:,2])) 

q <- df %>% ggplot() + 
        aes(x = m1, y = m2, color = z, shape = z) + 
        geom_point(size = 1) + 
        theme_bw() +
        xlab("1st embedding dimension") + 
        ylab("2nd embedding dimension") + 
        scale_color_discrete_qualitative(palette = "dynamic") + 
        scale_shape_manual(values = c(0, 1, 2)) + 
        theme(panel.border = element_rect(color = NA),
              legend.position = "bottom") + 
        ggtitle("(c). t-SNE Embedding") + 
        guides(legend = guide_legend(title = expression(italic(z))))
"""

function process(i)
    U = HypergraphNB.transform_eigenvector(real.(E_[:,i]), H);
    R"""
    df <- $U %>% 
        as.data.frame() %>%
        mutate(Node = row_number()) %>% 
        pivot_longer(!Node) %>% 
        mutate(Group = as.integer(sub(".", "", name))) %>% 
        select(-name) %>% 
        mutate(value = sign(value)) %>% 
        mutate(value = ifelse(value > 0, "+", "-")) %>% 
        mutate(evec = paste0("tilde(bold(x))[", $i, "]" ))
    """
    return @rget df
end

DF = vcat([process(i) for i in 1:size(E_)[2]]...)
R"""
df <- $DF

lines_df <- tibble(yintercept = $(n_nodes)*(1:$(n_groups)) + 0.5)

r <- df %>%
    rename(Sign = value) %>% 
    ggplot() + 
    aes(x = Group, y = Node, fill = Sign) + 
    geom_tile(alpha = .8) + 
    geom_hline(aes(yintercept = yintercept), data = lines_df, color = "black") +
    theme_bw() + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0), breaks = 1:10) + 
    scale_fill_manual(values = c("#f5b895", "#0f4c81")) +
    theme(panel.grid = element_blank(), 
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 15),
          legend.key = element_rect(colour = "black")) + 
    facet_grid(~evec, labeller = label_parsed) + 
    ggtitle("(b). Leading Eigenvectors") + 
    xlab("Cluster")
"""

R"""
library(patchwork)

p + r + q

ggsave("fig/algorithm-demo.png", dpi = 300, width = 10, height = 3.5)

p / r / q

ggsave("fig/algorithm-demo-vertical.png", dpi = 300, width = 5, height = 12)
"""
