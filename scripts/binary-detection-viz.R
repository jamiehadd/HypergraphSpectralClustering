library(tidyverse)
library(patchwork)
# ----------------------------------------------------------
# primary ellipse (detectability threshold)
# ----------------------------------------------------------

compute_c <- function(p, k, l, C){
    q <- 1/l
    (k-1)/q * (p + (1-p)*q*(1-q^(k-2))/(1-q^(k-1)))*C
}
compute_term <- function(p, k, l, C){
    c_i <- compute_c(p, k, l, C)
    l / (2*(l - 1)^2) * ((c_i - (k-1)*C)^2)/((k-1)*C)
}

compute_term_2 <- function(p, k, l, C){
    c_i <- compute_c(p, k, l, C)
    1 / (l - 1) * (c_i - (k-1)*C)
}

compute_primary_centroid <- function(k, l){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))
    -(r - q)/ (1 - r)
}

compute_primary_radius <- function(k, l, C){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))
    sqrt(2*(l-1)^2/l * q^2/((k-1)*C)*(1/(1-r)^2))
}

compute_primary_ellipse <- function(k_1, k_2, l, C_1, C_2){
    tibble(theta = seq(0, 2*pi, 0.001)) %>% 
    mutate(x = compute_primary_radius(k_1, l, C_1)*cos(theta) + compute_primary_centroid(k_1, l),
           y = compute_primary_radius(k_2, l, C_2)*sin(theta) + compute_primary_centroid(k_2, l)) %>% 
    arrange(desc(theta))
}

# ----------------------------------------------------------
# secondary ellipse (noise from eigenvalue collision)
# ----------------------------------------------------------

compute_secondary_centroid <- function(k, l, C){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))

    a <- l/((2*(l-1)^2)*((k-1)*C))
    b <- -1/(l-1)
    y <- (k-1)*C/q*(1-r)
    z <- -(k-1)*C*(r/q - 1)
    z/y - b/(2*a*y)
}

compute_secondary_radius <- function(k, l, C, m){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))
    a <- l/((2*(l-1)^2)*((k-1)*C))
    y <- (k-1)*C/q*(1-r)
    sqrt(m/(a*y^2))
}

compute_m <- function(k_1, k_2, l, C_1, C_2){
    a_1 <- l/((2*(l-1)^2)*((k_1-1)*C_1))
    a_2 <- l/((2*(l-1)^2)*((k_2-1)*C_2))
    b   <- -1/(l-1)
    1/4*b^2*(1/a_1 + 1/a_2)
}

compute_secondary_ellipse <- function(k_1, k_2, l, C_1, C_2){
    tibble(theta = seq(0, 2*pi, 0.01)) %>% 
    mutate(m = compute_m(k_1, k_2, l, C_1, C_2),
           x = compute_secondary_radius(k_1, l, C_1, m)*cos(theta) + compute_secondary_centroid(k_1, l, C_1),
           y = compute_secondary_radius(k_2, l, C_2, m)*sin(theta) + compute_secondary_centroid(k_2, l, C_2)) %>% 
    arrange(desc(theta))
}

# ----------------------------------------------------------
# function for elliptical heatmaps
# ----------------------------------------------------------

make_heatmap <- function(path, l, k_1, k_2, c_1, c_2, ax1, ax2){

    ellipse <- compute_primary_ellipse(k_1, k_2, l, c_1, c_2)

    secondary_ellipse <- compute_secondary_ellipse(k_1, k_2, l, c_1, c_2)


    df <- read_csv(path)

    df <- df %>%
        group_by(P_2, P_3, P_4) %>% 
        summarise(ARI = mean(ARI, na.rm = T))

    ax1 <- rlang::enquo(ax1)
    ax2 <- rlang::enquo(ax2)

    q <- df %>% 
        ggplot() + 
        geom_tile(aes(x = !!ax1, y = !!ax2, fill = ARI)) + 
        theme_bw() + 
        viridis::scale_fill_viridis(option = "inferno", limits = c(0, 1), na.value="black") + 
        scale_x_continuous(expand = c(0,0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0,0), limits = c(0, 1)) + 
        xlab(expression(italic(p)[2])) + 
        ylab(expression(italic(p)[3])) + 
        theme(strip.background = element_blank(),
            panel.spacing.x = unit(8, "mm"),
            strip.text = element_text(size = 10)) +
        coord_fixed() 

    p <- q + 
        geom_path(aes(x = x, y = y), data = secondary_ellipse,  linetype = "dashed") + 
        geom_path(aes(x = x, y = y), data = ellipse, color = "white") 

    return(list(p, q))
}

# ----------------------------------------------------------
# functions for affine heatmaps
# ----------------------------------------------------------

compute_c <- function(p, k, l, C){
    q <- 1/l
    return((k-1)/q * (p + (1-p)*q*(1-q^(k-2))/(1-q^(k-1)))*C)
}

compute_lambda_2 <- function(p, k, l, C){
    q <- 1/l
    c_i <- compute_c(p, k, l, C)
    return((c_i - (k-1)*C) / (2*(1-q)))
}

compute_lambda_1 <- function(k, C){
    return(C*(k-1))
}

compute_intercept <- function(k, l, C){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))
    b <- (k-1)*C/q*r
    b_ <- l/(2*(l-1))*(b - (k-1)*C)
    return(b_)
}

compute_slope <- function(k, l, C){
    q <- 1/l
    r <- q*(1-q^(k-2))/(1-q^(k-1))
    a <- (k-1)/q*(1 - r)*C

    a_ <- l/(2*(l - 1))*a
    return(a_)
}

make_affine_heatmap <- function(path, l, k_1, k_2, c_1, c_2){

    df <- expand.grid(p_1 = seq(0, 1, 0.01), p_2 = seq(0, 1, 0.01)) %>% 
        tibble() %>% 
        mutate(l_1 = compute_lambda_2(p_1, k = k_1, l = l, C = c_1), 
            l_2 = compute_lambda_2(p_2, k = k_2, l = l, C = c_2),
            L = l_1 + l_2,
            thresh = sqrt(compute_lambda_1(k_1, c_1) + compute_lambda_1(k_2, c_2))) 


    a <- -compute_slope(k_1, l, c_1) / compute_slope(k_2, l, c_2)

    b1 <- compute_intercept(k_1, l, c_1)
    b2 <- compute_intercept(k_2, l, c_2)
    thresh <- sqrt(compute_lambda_1(k_1, c_1) + compute_lambda_1(k_2, c_2))

    b <- (thresh - b1 - b2)/compute_slope(k_2, l, c_2)
    b_ <- (-thresh - b1 - b2)/compute_slope(k_2, l, c_2)

    df <- read_csv(path)

    df <- df %>%
        group_by(P_2, P_3, EV) %>% 
        summarise(ARI = mean(ARI, na.rm = T))

    p <- df %>% 
        ggplot() + 
        aes(x = P_2, y = P_3, fill = ARI) + 
        geom_tile() + 
        theme_bw() + 
        viridis::scale_fill_viridis(option = "inferno", limits = c(0, 1), na.value = "black") + 
        scale_x_continuous(expand = c(0,0)) + 
        scale_y_continuous(expand = c(0,0)) + 
        xlab(expression(italic(p)[2])) + 
        ylab(expression(italic(p)[3])) + 
        theme(panel.grid = element_blank(), 
        panel.spacing = unit(2, "lines"),
        strip.background = element_blank(),
        panel.spacing.x = unit(8, "mm"),
        strip.text = element_text(size = 10)) + 
        coord_fixed() + 
        guides(fill = guide_colorbar(title = "ARI"))
        
    q <- p + 
        geom_abline(aes(slope = a, intercept = b), color = "white") + 
        geom_abline(aes(slope = a, intercept = b_), color = "white") 

    return(list(q, p))

}

# ----------------------------------------------------------
# MAIN COMPUTATION: plot construction
# ----------------------------------------------------------

# first heatmap

l    <- 2
k_1  <- 2
k_2  <- 3
c_1  <- 5
c_2  <- 5
path <- "throughput/bulk-throughput/exp-1.csv"

V_1  <- make_heatmap(path, l, k_1, k_2, c_1, c_2, P_2, P_3)

# second heatmap (based on third experiment)

l    <- 2
k_1  <- 3
k_2  <- 4

c_1  <- 3
c_2  <- 3

path <- "throughput/bulk-throughput/exp-3.csv"

V_3  <- make_heatmap(path, l, k_1, k_2, c_1, c_2, P_3, P_4)

# third heatmap (based on fourth experiment)

l    <- 2
k_1  <- 2
k_2  <- 3

c_1  <- 5
c_2  <- 50

path <- "throughput/bulk-throughput/exp-4.csv"

V_4  <- make_heatmap(path, l, k_1, k_2, c_1, c_2, P_2, P_3)

# affine heatmap 

l    <- 2
k_1  <- 2
k_2  <- 3

c_1  <- 5
c_2  <- 5

path <- "throughput/vanilla-heatmap.csv"

V_0 <- make_affine_heatmap(path, l, k_1, k_2, c_1, c_2)



V_0[[1]] + ggtitle("(a).") + 
    V_1[[1]] + ggtitle("(b).") +
    V_3[[1]] + ggtitle("(c).") +
    V_4[[1]] + ggtitle("(d).") +
    plot_layout(guides = 'collect') 

ggsave("fig/4-heatmaps.png", dpi = 300, height = 7, width = 8)