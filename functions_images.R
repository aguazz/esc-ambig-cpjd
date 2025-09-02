# functions_images.R
# Assumes functions_sovler.R is sourced by the runner.

# shared helpers ---------------------------------------------------------------
.with_plot_margins <- function(expr) { 
  op <- par(mar=c(2,2,0.5,0.5)+0.2)
  on.exit(par(op), add=TRUE)
  force(expr) 
  }
render_and_save <- function(fname_base, plotfun, show=interactive(), save=TRUE,
                            dir="figures", width_px=2000, height_px=1400, res=300) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  if (isTRUE(show)) .with_plot_margins(plotfun())
  if (isTRUE(save)) {
    png(file.path(dir, paste0(fname_base,".png")), width=width_px, height=height_px, res=res); .with_plot_margins(plotfun()); dev.off()
    pdf(file.path(dir, paste0(fname_base,".pdf")), width=width_px/300*3.2, height=height_px/300*3.2); .with_plot_margins(plotfun()); dev.off()
  }
}

# fig1: H with thresholds ------------------------------------------------------
draw_fig1 <- function(sol, params, show=interactive(), save=TRUE, 
                      out_base="fig1_H_thresholds", dir="figures") {
  xs   <- unlist(sol$thresholds)[c("xL","x_kappa","x_lambda","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22"); ltys <- c(1,2,2,1)
  
  top_blank <- 0.075; bottom_blank <- 0.05
  frac_band <- 1 - top_blank - bottom_blank
  Hband <- params$l - (-params$u); Htot <- Hband/frac_band
  ylim_y <- c(-params$u - Htot*bottom_blank, params$l + Htot*top_blank)
  
  pad_x_frac <- 0.10; W <- xU - xL
  x_from <- xL - pad_x_frac*W; x_to <- xU + pad_x_frac*W
  
  plotfun <- function() {
    op <- par(xaxs="i"); on.exit(par(op), add=TRUE)
    curve(sol$functions$H(x), from=x_from, to=x_to, n=600,
          xlab="", ylab=expression(H(x) == V*minute*(x)),
          col="black", lwd=1.4, ylim=ylim_y)
    segments(x0=xs, y0=rep(-params$u, length(xs)), x1=xs, y1=rep(params$l, length(xs)),
             lty=ltys, col=cols, lwd=2)
    abline(h=-params$u, lty=3, col="#777777")
    abline(h=0,        lty=3, col="#BBBBBB")
    abline(h=params$l, lty=4, col="#777777")
    legend("topleft",
           legend=expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -u, l),
           lty   =c(1,       1,       2,        2,        1,       3,   4),
           col   =c("black", cols,    "#777777", "#777777"),
           lwd   =c(1.4,     2,       2,        2,        2,       2,   2),
           bty="n", horiz=TRUE, x.intersp=0.5, seg.len=2)
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}

# simulation + fig2/fig3 -------------------------------------------------------
simulate_reflected_jd <- function(T=8, dt=0.001, x0=NULL, params, thresholds) {
  n <- as.integer(round(T/dt)); tt <- seq(0, T, length.out=n+1)
  xL <- thresholds$xL; xU <- thresholds$xU; if (is.null(x0)) x0 <- thresholds$x_kappa
  EY <- -1/params$mu; drift <- params$b - params$r * EY
  X <- numeric(n+1); X[1] <- x0; U_proc <- L_proc <- numeric(n+1)
  jumped <- logical(n+1); dJ <- numeric(n+1)
  for (i in 1:n) {
    dW <- sqrt(dt)*rnorm(1); Nj <- rpois(1, params$r*dt)
    dJ[i+1] <- if (Nj>0) -sum(rexp(Nj, rate=params$mu)) else 0
    x_star <- X[i] + drift*dt + params$sigma*dW + dJ[i+1]
    U_inc <- L_inc <- 0
    if (x_star < xL) { U_inc <- xL - x_star; X[i+1] <- xL }
    else if (x_star > xU) { L_inc <- x_star - xU; X[i+1] <- xU }
    else X[i+1] <- x_star
    U_proc[i+1] <- U_proc[i] + U_inc; L_proc[i+1] <- L_proc[i] + L_inc
    jumped[i+1] <- Nj > 0
  }
  list(time=tt, X=X, U=U_proc, L=L_proc, jumped=jumped, dJ=dJ)
}

draw_fig2 <- function(sol, params, sim=NULL, seed=123, show=interactive(), save=TRUE,
                      out_base="fig2_reflected_path", dir="figures",
                      top_blank=0.075, bottom_blank=0.05) {
  xs <- sol$thresholds
  if (is.null(sim)) { 
    set.seed(seed); 
    sim <- simulate_reflected_jd(params=params, thresholds=xs) 
    }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$x_kappa; xl <- xs$x_lambda
  Hcorr <- xU - xL; frac <- 1 - top_blank - bottom_blank; Htot <- Hcorr/frac
  ylim <- c(xL - Htot*bottom_blank, xU + Htot*top_blank)
  plotfun <- function() {
    plot(t, X, type="n", xlab="", ylab="", ylim=ylim)
    usr <- par("usr"); rect(usr[1], xL, usr[2], xU, col=adjustcolor("gray85", 0.6), border=NA)
    lines(t, X, lwd=1.2)
    abline(h=c(xL, xU), col=c("#B22222", "#228B22"), lty=c(1,1), lwd=2)
    abline(h=c(xk, xl), col=c("#1E90FF", "#8A2BE2"), lty=c(2,2), lwd=2)
    pts <- which(sim$jumped); if (length(pts)) points(t[pts], X[pts], pch=16, cex=0.35)
    legend("topleft",
           legend=expression(bar(X)[t], underline(x), x^kappa, x^lambda, bar(x)),
           lty   =c(1,            1,          2,        2,        1),
           col   =c("black",      "#B22222",  "#1E90FF","#8A2BE2","#228B22"),
           lwd   =c(1.2,          2,          2,        2,        2),
           bty="n", horiz=TRUE, x.intersp=0.5, seg.len=2)
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}

draw_fig3 <- function(sim, show=interactive(), save=TRUE, out_base="fig3_controls_combined", dir="figures") {
  plotfun <- function() {
    plot(sim$time, sim$U, type="s", xlab="time", ylab="cumulative push", lwd=1.6)
    lines(sim$time, sim$L, type="s", lwd=1.6, lty=2)
    legend("topleft", legend=c(expression(U[t]~"(pushes up)"), expression(L[t]~"(pushes down)")),
           lty=c(1,2), lwd=2, bty="n")
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}
