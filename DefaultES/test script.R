df_raw = genCorData(n=100000,
                    sigma = .20,
                    mu = c(3.9,4.2),
                    rho = .9,
                    corstr = "cs"
                    cnames = c("pre","post"))
cor.test(~pre+post,data=df_raw)

mu <- c(3, 8, 15)
sigma <- c(1, 2, 3)

corMat <- matrix(c(1, .2, .8, .2, 1, .6, .8, .6, 1), nrow = 3)

dtcor1 <- genCorData(1000, mu = mu, sigma = sigma, rho = .7, corstr = "cs")

mu2 = c(1, 3)
sigma2 = 2
corMat2 <- matrix(c(1, .8,.8,1), nrow = 2)
covMat2 = sigma2*corMat2
dtcor2 <- mvrnorm(n=10,mu2,Sigma=covMat2) %>%
  as.data.frame()


round(cor(x=dtcor2$V1,y=dtcor2$V2),2)


var_raw = .25^2
cor_raw = matrix(c(1,.8,.8,1),nrow=2)
Sigma_raw = cor_raw*var_raw
df_raw = mvrnorm(
  n = 10,
  Sigma = Sigma_raw,
  mu = c(3.9, 4.2)
  ) %>% 
  as.data.frame() %>%
  rename(pre = V1,
         post = V2)

#Add a change score column
df_raw = df_raw %>% mutate(pre = round(pre,2),
                           post = round(post,2)) %>%
  mutate(diff = post-pre)

cor_prepost = cor(df_raw$pre,df_raw$post,
                  method = "pearson")

SMD <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  sd_pre = sd(d$pre, na.rm = TRUE)
  sd_post = sd(d$post, na.rm = TRUE)
  sd_av = (sd_pre+sd_post)/2
  m_diff = mean(d$diff, na.rm = TRUE)
  sd_diff = sd(d$diff, na.rm = TRUE)
  cor_prepost = cor(d$pre,d$post,
                    method = "pearson")
  dz = m_diff/sd_diff
  dav = m_diff/sd_av
  glass = m_diff/sd_pre
  drm = m_diff/sqrt((sd_pre^2+sd_post^2)-(2*cor_prepost*sd_pre*sd_post))*sqrt(2*(1-cor_prepost))
  CLES = pnorm(dz)
  result = c(dz,drm,dav,glass,CLES)
  
} 

raw_boot = boot(df_raw, SMD, R = 2000)
dz_raw = boot.ci(raw_boot, type="bca", index=1)
