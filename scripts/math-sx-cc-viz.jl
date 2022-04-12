using DelimitedFiles
using Random 
using Revise
using RCall
using MultivariateStats
using LinearAlgebra
using TSne
using Clustering
using StatsBase
# Random.seed!(54321)


function readThroughput(base_dir, suffix)
    path = "$base_dir/$suffix"
    M = readdlm("$path/M.txt", Int64);
    labels = readdlm("$path/labels.txt", String)[:,2];
    return M, labels
end

function laplacianCluster(M, n_groups,  km_iter, n_vecs = n_groups)
    D = diagm(vec(sum(M, dims = 1)))
    L = inv(D)*(D - M)
    E = eigvecs(L)
    E = real.(E)
    V = E[:,2:(1+n_vecs)]

    obj = 1

    tot_ss = (V .- mean(V, dims = 1)).^2 |> sum
    z = rand(1:n_groups, size(M)[1]);

    for i ∈ 1:km_iter
        clus = Clustering.kmeans(V', n_groups)
        cost_ss = clus.totalcost
        if cost_ss / tot_ss < obj
            z = Clustering.assignments(clus);
            obj = cost_ss / tot_ss
        end    
    end

    v1, v2 = V[:,1], V[:,2]
    return v1, v2, z
end

# GRAPH

base_dir = "throughput/math-sx-cc"
suffix = "graph"
M, labels = readThroughput(base_dir, suffix)


n_groups = 4
v1, v2, z = laplacianCluster(M, n_groups, 1000, 4)

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)


df <- tibble(m1 = $v1, m2 = $v2, label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = factor(group))

examples <- df %>% 
    filter(number <= 30) 

# df <- df %>% 
#     filter(m2 <= max(examples$m2),
#            m2 >= min(examples$m2),
#            m1 <= -0.005)

p <- df %>% 
    ggplot() + 
    aes(x = m1, y = m2) + 
    geom_point(aes(color = group), pch = 21) + 
    theme_minimal() + 
    geom_label_repel(aes(label = label, fill = group), data = examples, force = 10, max.overlaps = 30) + 
    guides(color = F, fill = F) + 
    ggtitle("(a). Projected Graph") + 
    scale_fill_discrete_qualitative(palette = "Pastel 1") + 
    scale_color_discrete_qualitative(palette = "Pastel 1") + 
    xlab("First principal component") + 
    ylab("Second principal component") + 
    scale_y_continuous(expand = c(-0.02, 0.02))
ggsave("fig/clustering-math-graph-cc.png", p, width = 6, height = 5, dpi = 300)
"""




# # HYPERGRAPH

suffix = "hypergraph"
M, labels = readThroughput(base_dir, suffix)


v1, v2, z = laplacianCluster(M, n_groups, 1000, 4)

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)


df <- tibble(m1 = $v1, m2 = $v2, label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = factor(group))

examples <- df %>% 
    filter(number < 30) 

# df <- df %>% 
#     filter(m2 <= max(examples$m2),
#            m2 >= min(examples$m2),
#            m1 <= -0.005)

q <- df %>% 
    ggplot() + 
    aes(x = m1, y = m2) + 
    geom_point(aes(color = group), pch = 21) + 
    theme_minimal() + 
    geom_label_repel(aes(label = label, fill = group), data = examples, force = 10, max.overlaps = 30) + 
    guides(color = F, fill = F) + 
    ggtitle("(b). Hypergraph") + 
    scale_fill_discrete_qualitative(palette = "Pastel 1") + 
    scale_color_discrete_qualitative(palette = "Pastel 1") + 
    xlab("First principal component") + 
    ylab("Second principal component") + 
    scale_y_continuous(expand = c(-0.02, 0.02))
ggsave("fig/clustering-math-hypergraph-cc.png", q, width = 6, height = 5, dpi = 300)
"""

## VANILLA 

suffix = "hypergraph-vanilla-5"
M, labels = readThroughput(base_dir, suffix)
v1, v2, z = laplacianCluster(M, n_groups, 1000, 4)

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)


df <- tibble(m1 = $v1, m2 = $v2, label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = factor(group))

examples <- df %>% 
    filter(number <= 30) 

r1 <- df %>% 
    ggplot() + 
    aes(x = m1, y = m2) + 
    geom_point(aes(color = group), pch = 21) + 
    theme_minimal() + 
    geom_label_repel(aes(label = label, fill = group), data = examples, force = 100, max.overlaps = 30, force_pull = 10) + 
    guides(color = F, fill = F) + 
    ggtitle("(c). NBHSC, 5 Eigenvectors") + 
    scale_fill_discrete_qualitative(palette = "Pastel 1") + 
    scale_color_discrete_qualitative(palette = "Pastel 1") + 
    xlab("First principal component") + 
    ylab("Second principal component") + 
    scale_y_continuous(expand = c(-0.02, 0.02))
ggsave("fig/clustering-math-hypergraph-vanilla-cc-5.png", q, width = 6, height = 5, dpi = 300)
"""


suffix = "hypergraph-vanilla-10"
M, labels = readThroughput(base_dir, suffix)
v1, v2, z = laplacianCluster(M, n_groups, 1000, 4)

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)


df <- tibble(m1 = $v1, m2 = $v2, label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = factor(group))

examples <- df %>% 
    filter(number <= 30) 

r2 <- df %>% 
    ggplot() + 
    aes(x = m1, y = m2) + 
    geom_point(aes(color = group), pch = 21) + 
    theme_minimal() + 
    geom_label_repel(aes(label = label, fill = group), data = examples, force = 100, max.overlaps = 30, force_pull = 10) + 
    guides(color = F, fill = F) + 
    ggtitle("(d). NBHSC, 10 Eigenvectors") + 
    scale_fill_discrete_qualitative(palette = "Pastel 1") + 
    scale_color_discrete_qualitative(palette = "Pastel 1") + 
    xlab("First principal component") + 
    ylab("Second principal component") + 
    scale_y_continuous(expand = c(-0.02, 0.02))
ggsave("fig/clustering-math-hypergraph-vanilla-cc-10.png", q, width = 6, height = 5, dpi = 300)
"""

R"""
library(patchwork)

p + q + r1 + r2

ggsave("fig/clustering-math-cc.png", width = 12, height = 10, dpi = 300)
"""

