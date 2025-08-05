# ---------- SET-UP mfis_uni_sim -----------------------------------------------------------

library(MASS)      # for future tweaks
set.seed(999)

# ---- population parameters --------------------------------------------------
mu  <- list(
  MFIS = 32.7,  sd_MFIS = 10,     # tighter SD
  PSS  = 19.5, sd_PSS  = 6,
  Age  = 43,   sd_Age  = 10
)

les_mu <- c(PV = 63.4, JC = 14.4, SC = 18.1, IT = 4.2)
les_sd <- c(PV = 32,   JC = 7,    SC = 9,    IT = 2)
les_mu <- c(les_mu, rep(mean(les_mu), 8))          # Loc5–Loc12
les_sd <- c(les_sd, rep(13, 8))

# helper
rtn <- function(n, m, s) pmin(100, pmax(0, rnorm(n, m, s)))
f2_to_beta <- function(f2) sqrt(f2)                 # assumes predictors z-scored

# ---- simulation function ----------------------------------------------------
sim_power <- function(N, f2, reps = 2000, which_les = 1) {
  beta_int <- -f2_to_beta(f2)                       # negative slope for demo
  hits <- 0
  for (i in 1:reps) {
    Age <- rnorm(N, mu$Age,  mu$sd_Age)
    Sex <- rbinom(N, 1, 0.5)
    PSS <- rnorm(N, mu$PSS, mu$sd_PSS)
    Les <- rtn(N, les_mu[which_les], les_sd[which_les]) / 100
    
    MFIS <- rnorm(N, mu$MFIS, mu$sd_MFIS)
    
    MFIS <- MFIS + beta_int * scale(PSS)[,1] * scale(Les)[,1]
    
    dat <- data.frame(MFIS, PSS, Les, Age, Sex)
    pval <- summary(
      lm(MFIS ~ PSS + Les + PSS:Les + Age + Sex, data = dat)
    )$coeff["PSS:Les","Pr(>|t|)"]
    hits <- hits + (pval < .05)
  }
  hits / reps
}

# ---- power grid -------------------------------------------------------------
Ns    <- c(50, 75, 100, 120)
f2s   <- c(0.02, 0.08, 0.25)             # small, medium, large
lesix <- 1:12

grid <- expand.grid(N = Ns, f2 = f2s, Lesion = lesix)
grid$Power <- mapply(sim_power, N = grid$N, f2 = grid$f2,
                     which_les = grid$Lesion)
print(grid)

write.csv(grid, file = "mfis_uni_interaction_output.csv", row.names = FALSE)