library(ggplot2)
rhos <- seq(-0.99,0.99,by=0.01)
sigma_pre <- seq(0.05,10,by=0.05)
nsim <- 500
combos <- expand.grid(rhos,sigma_pre)
colnames(combos) <- c("rho","sigma_pre")
set.seed(0)
vals <- pbmcapply::pbmclapply(1:nrow(combos), function(x) {
  rho <- combos[x,1]
  sigma_pre <- combos[x,2]
  Sigma <- matrix(c(1,rho,rho,1), nrow=2)*sigma_pre
  out <- sapply(1:nsim, function(i) {
    dat <- MASS::mvrnorm(20, c(0,1), Sigma)
    dz <- mean(dat[,2]-dat[,1])/sd(dat[,2]-dat[,1])
    delta <- mean(dat[,2]-dat[,1])/sd(dat[,1])
    c(dz,delta)
  })
  data.frame(rho = rho,
             sigma = sigma_pre,
             dz = mean(out[1,]),
             var_dz = var(out[1,]),
             delta = mean(out[2,]),
             var_delta = var(out[2,]))
}, mc.cores = parallel::detectCores())

out <- do.call(rbind.data.frame, vals)
head(out)

out.long <- tidyr::gather(out, key = statistic, value = estimate, dz, delta)
out.long$variance <- ifelse(out.long$statistic == "dz", out.long$var_dz, out.long$var_delta)
out.long[,c("var_dz", "var_delta")] <- NULL
out.long <- tidyr::gather(out.long, type, value, estimate, variance)

ggplot(out.long,
       aes(x = rho,
           y = sigma,
           fill = value)) +
  geom_raster() +
  facet_grid(type ~ statistic) +
  scale_fill_viridis_c(trans = "log10") +
  labs(y = expression(sigma["pre"]),
       x = "Correlation (ρ)",
       fill = "Value")
  