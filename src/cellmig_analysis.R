require(cellmig)

dir.create("results")

d <- get(load("data/dataset_1.RData"))
o <- cellmig(x = d[,c("well", "plate", "compound", "dose", "v", "offset")],
             control = list(mcmc_warmup = 1000,
                            mcmc_steps = 2000,
                            mcmc_chains = 4,
                            mcmc_cores = 4,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))
save(o, file = "results/cellmig_dataset_1.RData", compress = TRUE)
rm(o, d)


d <- get(load("data/dataset_2.RData"))
o <- cellmig(x = d[,c("well", "plate", "compound", "dose", "v", "offset")],
             control = list(mcmc_warmup = 1000,
                            mcmc_steps = 2000,
                            mcmc_chains = 4,
                            mcmc_cores = 4,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))
save(o, file = "results/cellmig_dataset_2.RData", compress = TRUE)
