library(cellmig)
library(ggplot2)
library(ggforce)
library(patchwork)
library(parallel)
library(peakRAM)

N_simulations <- 3
N_bioreps <- c(3, 6, 9, 12)
N_techreps <- c(3, 6, 9, 12)
N_cellreps <- 80
N_groups <- c(25, 50, 100)
chains <- 1

run_configure <- function(b, N_bio, N_tech, N_cells, N_group, chains) {
  control <- list(
    N_biorep = N_bio,
    N_techrep = N_tech,
    N_cell = N_cells,
    N_group = N_group,
    prior_alpha_p_M = -0.5,
    prior_alpha_p_SD = 0.2,
    prior_kappa_mu_M = 1.5,
    prior_kappa_mu_SD = 0.2,
    prior_kappa_sigma_M = 0,
    prior_kappa_sigma_SD = 0.2,
    prior_sigma_bio_M = 0,
    prior_sigma_bio_SD = 0.1,
    prior_sigma_tech_M = 0,
    prior_sigma_tech_SD = 0.1,
    prior_sigma_delta_M = 0,
    prior_sigma_delta_SD = 0.2)
  
  seed_id <- sample(x = 1000+1:10^6, size = 1)
  set.seed(seed_id)
  
  y_p <- gen_full(control = control)  
  
  sim_data <- y_p$y
  
  # format simulated data to be used as input in cellmig
  sim_data$well <- as.character(sim_data$well)
  sim_data$compound <- as.character(sim_data$compound)
  sim_data$plate <- as.character(sim_data$plate)
  sim_data$offset <- 0
  sim_data$offset[sim_data$group==1] <- 1
  
  time_cellmig <- system.time({
    ram_cellmig <- peakRAM(
      
      o <- cellmig(x = sim_data,
                   control = list(mcmc_warmup = 1000,
                                  mcmc_steps = 2000,
                                  mcmc_chains = chains,
                                  mcmc_cores = chains,
                                  mcmc_algorithm = "NUTS",
                                  adapt_delta = 0.9,
                                  max_treedepth = 10))
      
    )})
  
  return(list(b = b,
              N_bio = N_bio, 
              N_tech = N_tech, 
              N_cells = N_cells, 
              N_group = N_group, 
              chains = chains,
              time_cellmig = time_cellmig,
              ram_cellmig = ram_cellmig,
              control = control))
}

i <- 1
o <- vector(mode = "list", length = 1)
for(N_bio in N_bioreps) {
  for(N_tech in N_techreps) {
    for(N_cells in N_cellreps) {
      for(N_group in N_groups) {
        cat(i, "\n")
        o[[i]] <- mclapply(X = 1:N_simulations, 
                           FUN = run_configure, 
                           N_bio = N_bio, 
                           N_tech = N_tech, 
                           N_cells = N_cells, 
                           N_group = N_group,
                           chains = chains,
                           mc.cores = N_simulations)
        i <- i + 1
      }
    }
  }
}

save(o, file = "compute_benchmark.RData", compress = T)



o <- get(load("compute_benchmark.RData"))
res <- do.call(rbind, lapply(X = o, FUN = function(x) {
  return(do.call(rbind, lapply(X = x, FUN = function(x) {
    return(data.frame(b = x$b, 
                      N_bio = x$N_bio, 
                      N_tech = x$N_tech, 
                      N_cells = x$N_cells,
                      N_group = x$N_group,
                      N_chains = x$chains,
                      CPU = x$time_cellmig["elapsed"],
                      Peak_RAM = x$ram_cellmig$Peak_RAM_Used_MiB))
  })))
}))
res$N_bio_label <- paste0("Nbio = ", res$N_bio)
res$N_bio_label <- factor(x = res$N_bio_label, levels = unique(res$N_bio_label))
res$N_tech_label <- paste0("Ntech = ", res$N_tech)
res$N_tech_label <- factor(x = res$N_tech_label, levels = unique(res$N_tech_label))

fig_cpu <- ggplot(data = res)+
  facet_grid(N_tech_label~N_bio_label)+
  geom_line(aes(x = N_group, y = CPU/60^2, group = b), col = "darkgray")+
  geom_point(aes(x = N_group, y = CPU/60^2), size = 2, alpha = 0.5)+
  theme_bw(base_size = 10)+
  ylab(label = "Elapsed time (hours)")+
  xlab(label = "Treatments")

fig_ram <- ggplot(data = res)+
  facet_grid(N_tech_label~N_bio_label)+
  geom_line(aes(x = N_group, y = Peak_RAM/1024, group = b), col = "darkgray")+
  geom_point(aes(x = N_group, y = Peak_RAM/1024), size = 2, alpha = 0.5)+
  theme_bw(base_size = 10)+
  ylab(label = "Peak RAM (gigabytes)")+
  xlab(label = "Treatments")

fig_compute <- (fig_cpu|fig_ram)+patchwork::plot_annotation(tag_levels = "A")

ggsave(filename = "Supplementary_compute.pdf",
       plot = fig_compute,
       device = "pdf",
       width = 8,
       height = 6)
