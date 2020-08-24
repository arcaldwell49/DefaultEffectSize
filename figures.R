library(ggplot2)
library(ggnewscale)

# Set correlations and SDs for which to calculate SMDs
rhos <- seq(-0.99,0.99,by=0.01)
sigma_pre <- seq(0.05,10,by=0.05)
combos <- expand.grid(rhos,sigma_pre)
colnames(combos) <- c("rho","sigma_pre")

nsim <- 1000 # number of simulations per condition
set.seed(0) # set seed for reproducibility

# loop through all rho-SD_pre combinations
vals <- pbmcapply::pbmclapply(1:nrow(combos), function(x) {
  rho <- combos[x,1] # correlation between repeated measures
  sigma_pre <- combos[x,2] # SD_pre
  Sigma <- matrix(c(1,rho,rho,1), nrow=2)*sigma_pre^2 # create covariance matrix
  
  # for each combination, run nsim simulations
  out <- sapply(1:nsim, function(i) {
    # simulate from a bivariate distribution with n = 20/study
    dat <- MASS::mvrnorm(20, c(0,1), Sigma)
    rr <- cor(dat[,1],dat[,2]) # observed correlation between repeated measures
    
    diff <- mean(dat[,2]-dat[,1])
    sd_pre <- sd(dat[,1])
    
    # calculate SMDs
    dz <- diff/sd(dat[,2]-dat[,1])
    delta <- diff/sd_pre
    dav <- diff/((sd_pre + sd(dat[,2])) / 2)
    drm <- dz * sqrt(2*(1-rr))
    
    # return SMDs
    cbind(dz = dz,
          delta = delta,
          dav = dav,
          drm = drm)
  })
  
  # add all simulated SMDs to a data frame
  data.frame(rho = rho,
             sigma = sigma_pre,
             dz = out[1,],
             delta = out[2,],
             dav = out[3,],
             drm = out[4,])
}, mc.cores = parallel::detectCores())

out <- do.call(rbind.data.frame, vals) # list to data frame
saveRDS(out, "~/Documents/GitHub/DefaultEffectSize/sim_values.rds") # save all

# load all SMDs and aggregate to obtain estimates/variances
out <- readRDS("~/Documents/GitHub/DefaultEffectSize/sim_values.rds")
means <- aggregate(cbind(dz,delta,dav,drm) ~ rho + sigma, out, mean)
vars <- aggregate(cbind(dz,delta,dav,drm) ~ rho + sigma, out, var)
stacked <- rbind(cbind(means, statistic = "Estimate"),
                 cbind(vars, statistic = "Variance"))
out.long <- tidyr::gather(stacked, key = type, value = value, 3:6)
saveRDS(out.long, "~/Documents/GitHub/DefaultEffectSize/plot_values.rds")

# load estimates/variances only
out.long <- readRDS("~/Documents/GitHub/DefaultEffectSize/plot_values.rds")
# fix labels for plotting
out.long$type <- factor(out.long$type,
                        levels = unique(out.long$type),
                        labels = c(expression(d["z"]),
                                   expression(Delta),
                                   expression(d["av"]),
                                   expression(d["rm"])))
out.long$statistic <- factor(out.long$statistic,
                             levels = unique(out.long$statistic),
                             labels = c(expression(Estimate),
                                        expression(Standard~Error)))

# plot -- raster
ggplot() +
  geom_raster(data = subset(out.long, statistic == "Estimate"),
              aes(x = rho,
                  y = sigma,
                  fill = value)) +
  scale_fill_viridis_c(trans = "log10",
                       name = "Estimate",
                       option = "B",
                       breaks = c(0.1,1,10,100),
                       labels = c(0.1,1,10,100),
                       guide = guide_colorbar(order = 1,
                                              barwidth = 1)) +
  new_scale_fill() +
  geom_raster(data = subset(out.long, statistic == "Standard ~ Error"),
              aes(x = rho,
                  y = sigma,
                  fill = sqrt(value))) +
  scale_fill_gradientn(colors = pals::parula(),
                       name = "Standard Error",
                       trans = "log10",
                       #breaks = c(0.1,10,1000),
                       #labels = c(0.1,10,1000),
                       guide = guide_colorbar(order = 2,
                                              barwidth = 1)) +
  facet_grid(statistic ~ type,
             labeller = label_parsed) +
  labs(y = expression(sigma["pre"]),
       x = expression("Correlation" ~ (rho))) +
  scale_x_continuous(breaks = c(-1,0,1),
                     limits = c(-1,1)) +
  scale_y_continuous(limits = c(0,10)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

ggsave("smd_simulation.pdf",
       width = 10,
       height = 4.5)


#geom_line version version
ggplot() +
  geom_line(data = subset(out.long, statistic == "Estimate" & type %in% unique(out.long$type)),
              aes(x = rho,
                  y = value,
                  color = sigma,
                  group = sigma)) +
  scale_color_viridis_c(name = expression(sigma["pre"]),
                        trans = "log10",
                       option = "C",
                       #breaks = c(0.1,1,10,100),
                       #labels = c(0.1,1,10,100),
                       guide = guide_colorbar(order = 1,
                                              barwidth = 0.5)) +
  geom_line(data = subset(out.long, statistic == "Standard ~ Error" & type %in% unique(out.long$type)),
             aes(x = rho,
                  color = sigma,
                 group = sigma,
                  y = sqrt(value))) +
  facet_grid(statistic ~ type,
             labeller = label_parsed,
             scales = "free_y",
             switch = "y") +
  labs(y = "",
       x = expression("Correlation"~(italic(r)))) +
  scale_x_continuous(breaks = c(-1,0,1),
                     limits = c(-1,1)) +
  scale_y_continuous(trans="log10") +
  #scale_y_continuous(limits = c(0,10)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff")) +
  theme(strip.placement = "outside")


ggsave("smd_simulation_point.pdf",
       width = 8,
       height = 3)




library(pracma)
correct <- function(m) 1 - 3/(4 * m - 1)

# Set correlations and SDs for which to calculate SMDs
rhos <- seq(-0.99,0.99,by=0.01)
sigma <- pracma::logseq(0.05,10,20) #seq(0.05,10,by=0.05)
combos <- expand.grid(rhos,sigma)
colnames(combos) <- c("rho","sigma")
n <- 20
dz <- 1/sqrt(2*combos$sigma^2*(1-combos$rho))
vardz <- (n-1)/(n*(n-3)) * (1 + dz^2*n) - dz^2/correct(n-1)^2
delta <- 1/combos$sigma
vardelta <- (n-1)/(n*(n-3)) * (2*(1-combos$rho) + delta^2*n) - delta^2/correct(n-1)^2

closed_form <- rbind(cbind(combos,
                           data.frame(statistic = "Estimate",
                                      type = 'd["z"]',
                                      value = dz)),
                     cbind(combos,
                           data.frame(statistic = "Variance",
                                      type = 'd["z"]',
                                      value = vardz)),
                     cbind(combos,
                           data.frame(statistic = "Estimate",
                                      type = 'Delta["pre"]',
                                      value = delta)),
                     cbind(combos,
                           data.frame(statistic = "Variance",
                                      type = 'Delta["pre"]',
                                      value = vardelta)))

closed_form$statistic <- factor(closed_form$statistic,
                                levels = unique(closed_form$statistic),
                                labels = c(expression(Estimate),
                                           expression(Standard~Error)))

ggplot() +
  geom_line(data = subset(closed_form, statistic == "Estimate"),
            aes(x = rho,
                y = value,
                color = sigma,
                group = sigma)) +
  scale_color_viridis_c(name = expression(sigma["pre"]),
                        trans = "log10",
                        option = "C",
                        #breaks = c(0.1,1,10,100),
                        #labels = c(0.1,1,10,100),
                        guide = guide_colorbar(order = 1,
                                               barwidth = 0.5)) +
  geom_line(data = subset(closed_form, statistic == "Standard ~ Error"),
            aes(x = rho,
                color = sigma,
                group = sigma,
                y = sqrt(value))) +
  facet_grid(statistic ~ type,
             labeller = label_parsed,
             scales = "free_y",
             switch = "y") +
  labs(y = "",
       x = expression("Correlation"~(italic(r)))) +
  scale_x_continuous(breaks = c(-1,0,1),
                     limits = c(-1,1)) +
  scale_y_continuous(trans="log10") +
  #scale_y_continuous(limits = c(0,10)) +
  theme_classic() +
  #theme(panel.background = element_rect(fill = "#f9f9f9",color = "#ffffff")) +
  theme(aspect.ratio = 0.5,
        strip.placement = "outside",
        strip.background = element_rect(fill="white",color="white"),
        strip.text = element_text(face = "bold",
                                  size = 12),
        axis.title.x = element_text(face = "bold",
                                    size = 12),
        legend.title = element_text(face = "bold",
                                    size = 12),
        panel.border = element_rect(fill = NA,color = "lightgrey"))


ggsave("smd_closed_form.pdf",
       width = 6,
       height = 3)
