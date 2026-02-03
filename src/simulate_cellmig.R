.libPaths(new = "/mnt/nfs/simo/rpack_4_4/")
library(cellmig)
library(ggplot2)
library(ggforce)
library(patchwork)
library(parallel)

N_simulations <- 100
N_bioreps <- c(3, 6, 9, 12)
N_techreps <- c(3, 6, 9, 12)
N_cellreps <- c(25, 50, 100)
deltas <- log(c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2))
offset <- 6
chains <- 2
mc.cores <- 50

run_configure <- function(b, N_bio, N_tech, N_cells, deltas, offset, chains) {
    control <- list(
        N_biorep = N_bio,
        N_techrep = N_tech,
        N_cell = N_cells,
        delta = deltas,   # Treatment effects
        sigma_bio = 0.1,  # Biological variability
        sigma_tech = 0.1,# Technical variability
        offset = offset,  # Systematic offset (explain usage if relevant)
        prior_alpha_p_M = -0.3,    # Mean for plate effects
        prior_alpha_p_SD = 0.4,   # SD for plate effects
        prior_kappa_mu_M = 1.6,   # Mean for kappa location
        prior_kappa_mu_SD = 0.3,  # SD for kappa location
        prior_kappa_sigma_M = 0,  # Mean for kappa scale
        prior_kappa_sigma_SD = 0.1# SD for kappa scale
    )
    
    seed_id <- sample(x = 1000+1:10^6, size = 1)
    set.seed(seed_id)
    
    y_p <- gen_partial(control = control)  
    
    sim_data <- y_p$y
    
    # format simulated data to be used as input in cellmig
    sim_data$well <- as.character(sim_data$well)
    sim_data$compound <- as.character(sim_data$compound)
    sim_data$plate <- as.character(sim_data$plate)
    sim_data$offset <- 0
    sim_data$offset[sim_data$group==control$offset] <- 1
    
    o <- cellmig(x = sim_data,
                 control = list(mcmc_warmup = 300,
                                mcmc_steps = 1000,
                                mcmc_chains = chains,
                                mcmc_cores = 1,
                                mcmc_algorithm = "NUTS",
                                adapt_delta = 0.8,
                                max_treedepth = 10))
    
    og <- o$posteriors$delta_t
    og$real_delta <- deltas[-offset]
    og$fp <- (og$real_delta >= og$X2.5. & og$real_delta <= og$X97.5.) & 
        !(og$X2.5. <= 0 & og$X97.5. >= 0)
    og$b <- b
    og$N_bio <- N_bio
    og$N_tech <- N_tech
    og$N_cells <- N_cells
    og$seed_id <- seed_id
    
    return(og)
}

os <- c()
for(N_bio in N_bioreps) {
    for(N_tech in N_techreps) {
        for(N_cells in N_cellreps) {
            cat("N_bio:", N_bio, " N_tech:", N_tech, " N_cells:", N_cells, "\n")
            o <- mclapply(X = 1:N_simulations, 
                          FUN = run_configure, 
                          N_bio = N_bio, N_tech = N_tech, N_cells = N_cells, 
                          deltas = deltas, offset = offset, chains = chains, 
                          mc.cores = mc.cores)
            o <- do.call(rbind, o)
            os <- rbind(os, o)
        }
    }
}

save(os, file = "/mnt/nfs/simo/cellmig_application/results/benchmark_B100.RData", compress = T)


dim(os)