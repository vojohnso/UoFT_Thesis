# Stan models and scripts

Run everything from the `Simulation/` folder (paths are relative to it).

```
stan/
  models/     the Stan programs (.stan)
  R/          helper functions: prep_data.R (data lists for Stan), plot_helpers.R, compile_model.R
  scripts/    runnable scripts: run_* (single dataset), sim_* (replicate studies), summarize_*
  compiled/   cached compiled models (.rds, created automatically, safe to delete)
output/stan/  figures, tables and saved results from these scripts
```

## Models (`stan/models/`)
| File | What it is |
|---|---|
| `gcomp_known_exposure.stan` | Poisson outcome model with the true exposure A observed (also used for the naive fit, with Z as the exposure). Spline in A, g-computation in `generated quantities`. |
| `gcomp_latent_exposure.stan` | A latent, Z observed with known error SD; error is Z = A + noise, linear exposure effect. No `D`, no spatial term. |
| `gcomp_latent_exposure_D.stan` | As above with the Option A error model Z = beta_Z1 * A + beta_Z2 * D + noise; error parameters known (data). |
| `gcomp_latent_exposure_D_priors.stan` | Option A error model with beta_Z1, beta_Z2, sigma_Z estimated, plus a CAR spatial effect (plain `multi_normal_prec` form, slow at 400 areas). |
| `gcomp_latent_exposure_D_transparent.stan` | Same error model in a "transparent" parameterization (latent exposure on the Z scale; beta_Z1 only in generated quantities). Experimental. |

## Scripts (`stan/scripts/`)
| Script | What it does |
|---|---|
| `run_known_exposure.R` | Known-exposure model vs BART and the truth, three scenarios, one dataset. |
| `run_latent_exposure.R` | Oracle vs naive vs latent model, one dataset (linear scenario, error known). |
| `run_latent_exposure_D.R` | Latent model with the `D` error model: known parameters, then estimated parameters with CAR (100 areas). |
| `sim_latent_replicates.R [n_reps]` | Replicates of oracle / naive / latent (known error). Summarise with `summarize_latent_replicates.R`. |
| `sim_latent_misspec.R [n_reps]` | Misspecified exposure model (nonlinear exposure mean). |
| `sim_latent_D_replicates.R [n_reps] [iter]` | Replicates with the `D` error model and CAR (100 areas); summary printed at the end. |

Example: `Rscript stan/scripts/sim_latent_replicates.R 50`.
Models are compiled through `compile_model("<name>")` (see `R/compile_model.R`), which reads `models/<name>.stan` and caches the result in `compiled/`.
