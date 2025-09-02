# batch_render.R  — run many parameter sets defined inline
# --------------------------------------------------------
# Requires:
#   - functions.R (your solver: make_params(), solve_two_stage(), etc.)
#   - plots.R     (the plotting helpers: draw_fig1(), draw_fig2(), draw_fig3(),
#                  and simulate_reflected_jd())
# --------------------------------------------------------

source("functions_solver.R")
source("functions_images.R")

# Global preview/save toggles for the whole batch
SHOW <- TRUE   # preview in RStudio
SAVE <- TRUE    # write PNG+PDF to disk

# Helper to define a parameter set succinctly
pset <- function(name, b, delta, r, eps, sigma, mu, u, l,
                 seed = 123, init = list(xL = -2.0, gapU = 1.0)) {
  list(name=name, b=b, delta=delta, r=r, eps=eps, sigma=sigma, mu=mu, u=u, l=l,
       seed=seed, init=init)
}

# ==== Define all parameter sets here (inline) ==================================
param_sets <- list(
  pset("baseline",
       b=-1, delta=1.0, r=2, eps=0.75, sigma=1.5, mu=1, u=1.0, l=2.0,
       seed=123, init=list(xL=-2.0, gapU=1.0))
  
  # pset("paper_params",
  #      b=-0.5, delta=2, r=2, eps=0.5, sigma=1, mu=1, u=1, l=1,
  #      seed=321, init=list(xL=-2.0, gapU=1.0))
  
  # Add more sets as needed:
  # pset("high_sigma",  b=-0.5, delta=1.0, r=1.0, eps=0.30, sigma=1.4, mu=1.5, u=1, l=1),
  # pset("low_drift",   b=-1.2, delta=0.5, r=1.0, eps=0.30, sigma=0.8, mu=1.5, u=1, l=1)
)

# ==== Batch render =============================================================
for (ps in param_sets) {
  params <- make_params(b=ps$b, delta=ps$delta, r=ps$r, eps=ps$eps,
                        sigma=ps$sigma, mu=ps$mu, u=ps$u, l=ps$l)
  
  message(sprintf("\n== Solving set '%s' ==", ps$name))
  sol <- try(solve_two_stage(params, init = ps$init), silent = TRUE)
  
  if (inherits(sol, "try-error") || !isTRUE(sol$status$converged)) {
    warning(sprintf("Skipping '%s': solver failed or did not converge.", ps$name))
    next
  }
  
  outdir <- file.path("figures", ps$name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  set.seed(ps$seed)
  sim <- simulate_reflected_jd(params = params, thresholds = sol$thresholds)
  
  # Draw all figures for this parameter set
  draw_fig1(sol, params, show = SHOW, save = SAVE, dir = outdir)
  draw_fig2(sol, params, sim = sim, show = SHOW, save = SAVE, dir = outdir)
  draw_fig3(sim, show = SHOW, save = SAVE, dir = outdir)
  
  message(sprintf("Finished '%s' -> %s", ps$name, normalizePath(outdir)))
}

message("\nBatch complete.")
