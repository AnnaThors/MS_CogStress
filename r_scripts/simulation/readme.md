# Lesion × Stress Monte-Carlo Simulations

*(MFIS, SDMT; set.seed(999))*

This repo hosts the code and data that underpin the sample-size justification in the protocol / manuscript.

⸻

### Quick start

 1.  Clone the repo and enter the code folder
 2.  Install dependencies (only base R + MASS needed)
 3.  Run any script — each one writes its CSV grid next to the script
source("mfis_uni_sim.R")   # ≈ 40 s on a laptop
source("sdmt_uni_sim.R")   # ≈ 40 s
source("manova_full.R")    # ≈ 3 min

*All three scripts finish silently and create the files below.*

⸻

### File tree


├─ mfis_uni_sim.R          # Univariate MFIS precision / power

├─ sdmt_uni_sim.R          # Univariate SDMT precision / power

├─ manova_full.R           # MFIS + SDMT MANOVA interaction

├─ mfis_uni_grid.csv       # 144-row Monte-Carlo output

├─ sdmt_uni_grid.csv       # 144-row Monte-Carlo output

├─ manova_grid.csv         # 144-row Monte-Carlo output

└─ README.md               # You are here


⸻

**Cite / credit**

All simulation code and raw output are provided.
Please cite the accompanying manuscript if you reuse any part of this work.

⸻

**Contact**

For questions or bug reports, open an issue or email anna.thorsteinsdottir.24@um.edu.mt 

⸻

## Variable specification & reference list

(population draws used in every simulation)

| Model term / correlation | Population draw (μ ± σ, or p) | Source / rationale                          |
|--------------------------|-------------------------------|---------------------------------------------|
| MFIS-Total               | 32.7 ± 19.0                  | Grange 2025, Table 1                        |
| SDMT-Raw                 | 38.0 ± 14.0                  | Grange 2025 (same cohort)                   |
| PSS-10                   | 19.5 ± 6.0                   | Novak 2023 Israeli MS survey                |
| Age (y)                  | 43 ± 10                      | Iceland registry expectation                |
| Sex                      | 50 % female / 50 % male      | Fixed by PI instruction                     |
| Education (y)            | 15 ± 3                       | Compulsory 10 + upper-sec/uni mix           |
| PHQ-9 (8-item)           | 5 ± 5                        | Tardo 2022 depression cohort                |
| Lesion % Periventricular | 63.4 ± 43 %                  | Pongratz 2023 baseline volumes              |
| Lesion % Juxtacortical   | 14.4 ± 20 %                  | idem                                        |
| Lesion % Subcortical     | 18.1 ± 25 %                  | idem                                        |
| Lesion % Infratentorial  | 4.2 ± 10 %                   | idem                                        |
| Lesion % Loc5 – Loc12    | 25 ± 27 %                    | Grand mean ± SD (placeholder)               |
| MFIS ↔ SDMT              | r = –0.31                    | Maier 2023                                  |
| PSS-10 ↔ MFIS            | r ≈ 0.25                     | Novak 2023 (β ≈ .25)                        |
| PSS-10 ↔ SDMT            | r ≈ –0.10                    | Conservative assumption                     |
| Age ↔ SDMT               | r ≈ –0.30                    | Repeated BICAMS findings                    |
| Interaction slope (β_int)| –0.15 & –0.20 SD             | Cohen small-medium / medium                 |

**Reference list**

Chojdak-Łukasiewicz, J., Kołtuniuk, A., Sawicka, E., & Pokryszko-Dragan, A. (2025). Adherence to therapeutic recommendation in relapsing-remitting multiple sclerosis patients. Frontiers in Immunology, 16, 1545430. https://doi.org/10.3389/fimmu.2025.1545430

Grange, E., et al. (2025). Fatigue and its association with upper limb function in people with multiple sclerosis. Life, 13(11), 2132. https://doi.org/10.3390/life13112132

Maier, S., et al. (2023). Sociodemographic and clinical determinants of fatigue in multiple sclerosis. Life, 13(11), 2132. https://doi.org/10.3390/life13112132

Novak, A. M., & Lev-Ari, S. (2023). Resilience, stress, well-being, and sleep quality in multiple sclerosis. Journal of Clinical Medicine, 12(2), 716. https://doi.org/10.3390/jcm12020716

Pongratz, V., et al. (2023). Lesion location across diagnostic regions in multiple sclerosis. Brain, Behavior, & Immunity – Health, 26, 100519. https://doi.org/10.1016/j.bbih.2022.100519

Tardo, L. M., et al. (2022). Determining prevalence of depression and covariates of depression in a cohort of multiple sclerosis patients. Journal of Central Nervous System Disease, 14, 11795735221098144. https://doi.org/10.1177/11795735221098144
