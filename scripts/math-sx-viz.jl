using DelimitedFiles
using Random 
using Revise
using RCall
using MultivariateStats
using TSne

function readThroughput(base_dir, suffix)
    path = "$base_dir/$suffix"
    V = readdlm("$path/V.txt");
    z = readdlm("$path/z.txt");
    z = floor.(Int, z)[:,1]
    labels = readdlm("$path/labels.txt", String)[:,2];
    return V, z, labels
end



# graph

base_dir = "throughput/math-sx"
suffix = "graph"
V, z, labels = readThroughput(base_dir, suffix)

M = fit(PCA, V'; maxoutdim=2).proj

# M = tsne(V', 2, 0, 1000, 1000.0);

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)





df <- tibble(m1 = $(M[:,1]), m2 = $(M[:,2]), label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = factor(group))



examples <- df %>% 
    filter(number < 30) 



df <- df %>% 
    filter(m2 <= max(examples$m2),
           m2 >= min(examples$m2),
           m1 <= -0.005)

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
ggsave("fig/clustering-math-graph.png", p, width = 6, height = 5, dpi = 300)
"""

# HYPERGRAPH



suffix = "hypergraph"
V, z, labels = readThroughput(base_dir, suffix)

M = fit(PCA, V'; maxoutdim=2).proj

# M = tsne(V', 2, 0, 1000, 50.0);

R"""
library(tidyverse)
library(ggrepel)
library(colorspace)

label_map <- c(2L, 1L, 3L, 4L) # depends on run of the actual clustering routine. 

df <- tibble(m1 = $(M[:,1]), m2 = $(M[:,2]), label = $(labels), group = $z) %>% 
    mutate(number = row_number(),
           group = map_int(group, ~label_map[[.x]]),
           group = factor(group))



examples <- df %>% 
    filter(number < 30) 

df <- df %>% 
    filter(m2 <= max(examples$m2),
           m2 >= min(examples$m2),
           m1 <= -0.005)

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
ggsave("fig/clustering-math-hypergraph.png", q, width = 6, height = 5, dpi = 300)
"""

# combined

R"""
library(patchwork)

p + q 

ggsave("fig/clustering-math.png", width = 12, height = 6, dpi = 300)
"""

