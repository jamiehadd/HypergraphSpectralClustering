using Revise
using RCall
using HypergraphNB
using Arpack
using MultivariateStats
using StatsBase
using DataFrames



function createViz(data_name, figure_path; nev = 30, nev_to_plot = 10, show_PCA = false)

    path      = "throughput/data-throughput/"*data_name*".csv"
    proj_path = "throughput/data-throughput/"*data_name*"-projected.csv"


    label_path = "data/"*data_name*"/label-names-"*data_name*".txt"
    
    H, z = HypergraphNB.read_hypergraph_data(data_name);
    
    B = reducedBPJacobian(H, z);    
    E = Arpack.eigs(B; nev = nev, ritzvec = true)            
    V = hcat([HypergraphNB.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:nev]...)
    
    V_ = 1.0*(V .> 0)
    M = fit(PCA, V_; maxoutdim=4)
    
    K_ = sort(collect(keys(H.E)))

    n = length(H.N)

    c, G = degreeTensor(H, z);
    q = 1/n * StatsBase.counts(z)
    G_ = zero(G)
    for i ∈ 1:length(K_)
        G_[i,:,:] = (G[i,:,:] / ((K_[i] - 1) * c[i]) .- 1) .* q
    end
    
    DF = DataFrame()

    for i ∈ 1:length(unique(z))
        id = zeros(0)
        K = zeros(0)
        value = zeros(0)
        for k ∈ 1:size(G_, 1)
            push!(id, i)
            push!(value, G_[k,i,i])
            push!(K, K_[k])
        end
        df = DataFrame(id = id, value = value, K = K)
        DF = vcat(DF, df)
    end
        
    legend_rows = length(unique(z)) > 2 ? 2 : 1 
    
    R"""
    library(tidyverse)
    library(viridis)
    library(colorspace)
    library(patchwork)
    library(ggforce)
    
    df  <- read_csv($(path)) %>% mutate(projected = "Hypergraph")
    df_ <- read_csv($(proj_path)) %>% mutate(projected = "Projected graph")
    
    df <- df %>% 
        rbind(df_) %>% 
        mutate(obj = cost_SS / tot_SS) 
    
    best_df <- df %>% 
        filter(nev == $(nev_to_plot)) %>% 
        group_by(projected) %>% 
        filter(obj == min(obj, na.rm = T)) %>% 
        filter(row_number() == 1)
    
    p <- df %>% 
        filter(nev == $(nev_to_plot)) %>% 
        group_by(ngroups, projected) %>% 
        filter(obj == min(obj, na.rm = T)) %>% 
        filter(row_number() == 1) %>% 
        ggplot() + 
        geom_vline(aes(xintercept = factor($(length(unique(z))))), color = "grey") + 
        aes(x = factor(ngroups), y = ari, group = projected) +
        geom_line(aes(linetype = projected), color = "darkslategray") +     
        geom_point(color = "darkslategray", pch = 21, fill = "white") + 
        theme_bw() + 
        theme(strip.background = element_blank(),
              legend.position = "bottom",
              panel.border = element_rect(color = NA)) + 
        scale_color_brewer(palette = "Set2") +
        guides(color = guide_legend(title = element_blank()),
               linetype = guide_legend(title = element_blank())) + 
        xlab("Number of estimated groups") + 
        ylab("Adjusted Rand Index against ground truth") +
        scale_y_continuous(limits = c(NA, 1)) + 
        ggtitle(paste("Cluster recovery with", $(nev_to_plot), "eigenvectors")) 
    

    node_labels <- read_csv($(label_path), col_names = c("name")) %>% 
        mutate(id = row_number())

    if($show_PCA){
        df <- tibble(z = $z, m1 = $(M.proj[:,1]), m2 = $(M.proj[:,2])) 

        q <- df %>% 
            left_join(node_labels, by = c("z" = "id"), copy = TRUE) %>%  
            ggplot() + 
            aes(x = m1, y = m2, color = name) + 
            # ggforce::geom_mark_hull(aes(fill = name), alpha = 0.2) + 
            geom_point(pch = 21) + 
            theme_bw() +
            xlab("1st principal component") + 
            ylab("2nd principal component") + 
            guides(color = guide_legend(title = element_blank(), nrow = $legend_rows),
                fill = FALSE) + 
            scale_color_discrete_qualitative(palette = "dynamic") + 
            theme(panel.border = element_rect(color = NA),
                legend.position = "bottom") + 
            ggtitle("Eigenspace visualization (via PCA)")
    }else{
        q <- df %>% 
        filter(nev == $(nev_to_plot)) %>% 
        group_by(ngroups, projected) %>% 
        filter(obj == min(obj, na.rm = T)) %>% 
        filter(row_number() == 1) %>% 
        ggplot() + 
        geom_vline(aes(xintercept = factor($(length(unique(z))))), color = "grey") + 
        aes(x = factor(ngroups), y = obj, group = projected) +
        geom_line(aes(linetype = projected), color = "firebrick") +     
        geom_point(color = "firebrick", pch = 21, fill = "white") + 
        theme_bw() + 
        theme(strip.background = element_blank(),
              legend.position = "bottom",
              panel.border = element_rect(color = NA), 
              legend.text = element_text(size = 11)) + 
        scale_color_brewer(palette = "Set2") +
        guides(color = guide_legend(title = element_blank()),
               linetype = guide_legend(title = element_blank())) + 
        xlab("Number of groups ℓ") + 
        ylab("Proportion of within-group variance") +
        scale_y_continuous(limits = c(NA, 1)) + 
        ggtitle(paste("Scree plot with", $(nev_to_plot), "eigenvectors")) 
    }

    r <- $(DF) %>% 
        tibble() %>% 
        left_join(node_labels, by = c("id" = "id"), copy = TRUE) %>% 
        ggplot() + 
        aes(x = K, y = value, color = name, group = name) + 
        geom_line() + 
        geom_point(pch = 21, fill = "white") + 
        theme_bw() + 
        scale_color_discrete_qualitative(palette = "dynamic") + 
        guides(color = guide_legend(title = element_blank(), nrow = $legend_rows)) + 
        theme(panel.border = element_rect(color = NA),
              legend.position = "bottom") + 
        xlab(expression(paste("Edge size ", italic(k)))) + 
        ylab(expression(paste("Affinity ", italic(c)))) + 
        ggtitle("Within-group affinities")
    

    if($show_PCA){
        P <- p + q + r
    }else{
        P <- q + p + r
    }
    
    ggsave($(figure_path), P, width = 12, height = 4, dpi = 300)
            
    """
end


for i ∈ 1:2
    show_PCA = i == 2

    data_name = "SN-congress-bills"
    createViz(data_name, "fig/$(data_name)-$i.png", nev_to_plot = 10, show_PCA = show_PCA)

    data_name = "contact-high-school-classes"
    createViz(data_name, "fig/$(data_name)-$i.png", nev_to_plot = 30, show_PCA = show_PCA)

    data_name = "contact-primary-school-classes"
    createViz(data_name, "fig/$(data_name)-$i.png", nev_to_plot = 30, show_PCA = show_PCA)
end