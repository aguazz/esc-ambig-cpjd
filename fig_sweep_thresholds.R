# fig_sweep_thresholds.R
# Visualize how (xL, x^kappa, x^lambda, xU) move as b varies.

source("functions_solver.R")
source("functions_images.R")

sweep_thresholds_vs_b <- function(
    b_seq = seq(-1, 1, by = 0.1),
    # hold everything else fixed (you can change these)
    delta = 1.0, r = 2, eps = 0.75, sigma = 1.5, mu = 1, u = 1.0, l = 2.0,
    init = list(xL = -2.0, gapU = 1.0),
    out_dir = file.path("figures", "sweeps"),
    out_name = NULL,
    show = interactive(), save = TRUE
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  res <- data.frame(
    b = as.numeric(b_seq),
    xL = NA_real_, x_kappa = NA_real_, x_lambda = NA_real_, xU = NA_real_,
    converged = FALSE, message = NA_character_, stringsAsFactors = FALSE
  )
  
  last_init <- init
  last_sol  <- NULL
  
  for (i in seq_along(b_seq)) {
    b <- b_seq[i]
    message(sprintf("[b sweep] %2d/%d  b = % .3f ...", i, length(b_seq), b)); flush.console()
    
    
    params <- make_params(b = b, delta = delta, r = r, eps = eps,
                          sigma = sigma, mu = mu, u = u, l = l)
    
    # Warm start from the previous successful solve (helps continuation in b)
    if (!is.null(last_sol)) {
      last_init <- list(
        xL   = last_sol$thresholds$xL,
        gapU = max(1e-4, last_sol$thresholds$xU - last_sol$thresholds$x_lambda)
      )
    }
    
    sol <- try(solve_two_stage(params, init = last_init), silent = TRUE)
    
    if (!inherits(sol, "try-error") && isTRUE(sol$status$converged)) {
      th <- sol$thresholds
      res[i, c("xL","x_kappa","x_lambda","xU")] <-
        unlist(th[c("xL","x_kappa","x_lambda","xU")])
      res[i, "converged"] <- TRUE
      res[i, "message"]   <- sol$status$message
      last_sol <- sol
    } else {
      res[i, "message"] <- if (inherits(sol, "try-error")) as.character(sol) else sol$status$message
    }
  }
  
  # Save the raw data too (handy for tables/appendix)
  if (is.null(out_name))
    out_name <- sprintf("thresholds_vs_b_%g_to_%g_by_%g",
                        min(b_seq), max(b_seq), unique(diff(b_seq))[1])
  csv_path <- file.path(out_dir, paste0(out_name, ".csv"))
  write.csv(res, csv_path, row.names = FALSE)
  
  # Plot + save using your existing renderer (shows in RStudio and writes PNG+PDF)
  ok   <- res$converged
  yall <- as.numeric(unlist(res[ok, c("xL","x_kappa","x_lambda","xU")]))
  ylm  <- range(yall, na.rm = TRUE)
  
  cols <- c("#B22222", "#1E90FF", "#8A2BE2", "#228B22")
  ltys <- c(1, 1, 1, 1)
  
  plotfun <- function() {
    op <- par(xaxs = "i"); on.exit(par(op), add = TRUE)
    plot(res$b, res$xL, type = "n", xlab = expression(b), ylab = "threshold value", ylim = ylm)
    grid()
    lines(res$b, res$xL,       lwd = 2, col = cols[1], lty = ltys[1])
    lines(res$b, res$x_kappa,  lwd = 2, col = cols[2], lty = ltys[2])
    lines(res$b, res$x_lambda, lwd = 2, col = cols[3], lty = ltys[3])
    lines(res$b, res$xU,       lwd = 2, col = cols[4], lty = ltys[4])
    legend("topleft",
           legend = expression(underline(x), x^kappa, x^lambda, bar(x)),
           col = cols, lty = ltys, lwd = 2, bty = "n", horiz = TRUE, x.intersp = 0.5, seg.len = 2)
    title("Ambiguity thresholds & barriers vs drift b")
  }
  
  render_and_save(
    fname_base = out_name,
    plotfun    = plotfun,
    show       = show,
    save       = save,
    dir        = out_dir
  )
  
  message(sprintf("Saved plot(s) + data to: %s", normalizePath(out_dir)))
  list(
    results = res,
    files = list(
      csv = csv_path,
      png = file.path(out_dir, paste0(out_name, ".png")),
      pdf = file.path(out_dir, paste0(out_name, ".pdf"))
    )
  )
}

# ==== Example run ==============================================================
out <- sweep_thresholds_vs_b(b_seq = seq(-3.5, 3.5, by = 0.01))
View(out$results)   # or read.csv(out$files$csv)
