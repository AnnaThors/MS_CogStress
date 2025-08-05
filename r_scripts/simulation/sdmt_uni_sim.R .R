# ---------- Exploratory SDMT main-effect simulation (fixed column) -------------------------
library(MASS)
set.seed(999)

# --- population parameters --------------------------------------------------
mu <- list(
  SDMT = 38.0, sd_SDMT = 10,
  PSS  = 19.5, sd_PSS  = 6,
  Age  = 43,   sd_Age  = 10
)

# lesion distribution
les_mu <- c(PV = 63.4, JC = 14.4, SC = 18.1, IT = 4.2)
les_sd <- c(PV = 32,   JC = 7,    SC = 9,    IT = 2)
les_mu <- c(les_mu, rep(mean(les_mu), 8))     # Loc5–Loc12 placeholders
les_sd <- c(les_sd, rep(13, 8))

rtn <- function(n, m, s) pmin(100, pmax(0, rnorm(n, m, s)))

# --- simulation core --------------------------------------------------------
sim_effect <- function(N, beta_z, reps = 2000, tract = 1) {
  beta_hat <- numeric(reps); p_hit <- 0
  for (i in seq_len(reps)) {
    Age <- rnorm(N, mu$Age, mu$sd_Age)
    Sex <- rbinom(N, 1, 0.5)
    PSS <- rnorm(N, mu$PSS, mu$sd_PSS)
    
    Les_raw <- rtn(N, les_mu[tract], les_sd[tract]) / 100
    zLes    <- scale(Les_raw)[, 1]
    
    SDMT <- rnorm(N, mu$SDMT, mu$sd_SDMT) + beta_z * zLes
    
    fit <- lm(SDMT ~ PSS + zLes + Age + Sex)
    beta_hat[i] <- coef(fit)["zLes"]
    p_hit <- p_hit + (summary(fit)$coeff["zLes", "Pr(>|t|)"] < .05)
  }
  ci <- quantile(beta_hat, c(.025, .5, .975))
  true_f2 <- beta_z^2 / mu$sd_SDMT^2          # << fixed: theoretical f²
  data.frame(
    N,
    beta_injected = beta_z,
    median_beta   = ci[2],
    CI_low        = ci[1],
    CI_high       = ci[3],
    true_f2       = true_f2,
    power         = p_hit / reps,
    Tract         = tract
  )
}

# --- grid runner ------------------------------------------------------------
Ns     <- c(50, 75, 100, 120)
betas  <- c(-0.15, -0.20, -0.30)     # same tiers
tracts <- 1:12

results <- do.call(
  rbind,
  lapply(tracts, function(t)
    do.call(rbind,
            lapply(Ns, function(n)
              do.call(rbind,
                      lapply(betas, sim_effect,
                             N = n, tract = t)))))
)

print(results, row.names = FALSE)

write.csv(results, file = "sdmt_uni_output.csv", row.names = FALSE)