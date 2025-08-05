# -------- SET-UP ------------------------------------------------------------
library(MASS)          # for mvrnorm
set.seed(999)      # reproducible sass

# ---- user-editable population parameters -----------------------------------
mu  <- list(
  MFIS = 32.7,  sd_MFIS = 15.0,      # Grange 2025  [oai_citation:0‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12196498/)
  SDMT = 38.0,  sd_SDMT = 10.0,      # Grange 2025  [oai_citation:1‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12196498/)
  PSS  = 19.5,  sd_PSS  = 6.0,       # Novak 2023  [oai_citation:2‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9864697/)
  Age  = 43,    sd_Age  = 10,
  Edu  = 15,    sd_Edu  = 3,
  PHQ  = 5,     sd_PHQ  = 5          # Tardo 2022  [oai_citation:3‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9864697/?utm_source=chatgpt.com)
)

# Lesion-% means ± SD for the four empirical regions (Pongratz 2023)
les_mu <- c(PV = 63.4, JC = 14.4, SC = 18.1, IT = 4.2)         # %
les_sd <- c(PV = 32, JC = 7, SC = 9, IT = 2)        # crude %

# Pad eight placeholder tracts with the grand mean/SD
les_mu <- c(les_mu, rep(mean(les_mu), 8))                      # Loc5–Loc12
les_sd <- c(les_sd, rep(13, 8))   # Loc5–Loc12

# MFIS–SDMT r, plus tiny r’s elsewhere
r_mat <- matrix(0, 6, 6); colnames(r_mat) <- rownames(r_mat) <-
  c("MFIS","SDMT","PSS","Age","Edu","PHQ")
r_mat["MFIS","SDMT"] <- r_mat["SDMT","MFIS"] <- -0.31          # Maier 2023  [oai_citation:4‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10672347/?utm_source=chatgpt.com)

# Helper: sample truncated normal
rtn <- function(n, m, s) pmin(100, pmax(0, rnorm(n, m, s)))

# ---- simulation function ----------------------------------------------------
sim_power <- function(N, beta_int, reps = 2000, which_les = 1) {
  hits <- 0
  for (i in 1:reps) {
    # covariates & predictors
    Age <- rnorm(N, mu$Age,  mu$sd_Age)
    Sex <- rbinom(N, 1, 0.5)                                     # 50/50 split
    Edu <- rnorm(N, mu$Edu,  mu$sd_Edu)
    PHQ <- rnorm(N, mu$PHQ,  mu$sd_PHQ)
    PSS <- rnorm(N, mu$PSS,  mu$sd_PSS)
    Les <- rtn(N, les_mu[which_les], les_sd[which_les])/100      # 0–1 scale
    
    # outcomes (multivariate generation with preset r)
    eps  <- mvrnorm(N, mu = c(0,0), Sigma = matrix(c(1, r_mat[1,2],
                                                     r_mat[1,2], 1), 2, 2))
    MFIS <- scale(mu$MFIS + mu$sd_MFIS * eps[,1], scale = FALSE)[,1]
    SDMT <- scale(mu$SDMT + mu$sd_SDMT * eps[,2], scale = FALSE)[,1]
    
    # embed true effects (all continuous z-scored for ease)
    zPSS <- scale(PSS)[,1]; zLes <- scale(Les)[,1]
    int_term <- zPSS * zLes
    MFIS <- MFIS + beta_int * int_term
    SDMT <- SDMT + beta_int * int_term
    
    dat <- data.frame(MFIS, SDMT, PSS, Les, Age, Sex, Edu, PHQ)
    fit <- lm(cbind(MFIS, SDMT) ~ PSS + Les + PSS:Les + Age + Sex + Edu + PHQ, data = dat)
    pval <- summary(manova(fit))$stats["PSS:Les","Pr(>F)"]
    hits <- hits + (pval < .05)
  }
  hits/reps
}

# ---- power grid -------------------------------------------------------------
Ns     <- c(50,75,100,120)
betas  <- c(-0.15, -0.20, -0.30)            # small-medium, medium, large
les_ix <- 1:12                       # PV, JC, SC, IT, Loc5-Loc12

grid <- expand.grid(N = Ns, beta = betas, Lesion = les_ix)
grid$Power <- mapply(sim_power, N = grid$N, beta_int = grid$beta,
                     which_les = grid$Lesion)
print(grid)

# Save results to CSV
write.csv(grid, file = "manova_interaction_sim.csv", row.names = FALSE)