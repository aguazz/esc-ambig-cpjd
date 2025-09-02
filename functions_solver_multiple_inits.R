# ==============================================================================
# Two-stage solver for ergodic singular control with ambiguous jump diffusion
# Stage 1 (inner): given xL, compute gamma(xL), solve x_kappa and x_lambda
#                  using a scaled IH(x_lambda)=0 and a provisional H2 that
#                  matches level & slope at x_kappa (stable).
# Stage 2 (outer): solve for (xL, xU) so that smooth fit holds at x_kappa, x_lambda.
# ==============================================================================

# ----------------------- PARAMETERS -------------------------------------------
make_params <- function(b, delta, r, eps, sigma, mu, u, l) {
  stopifnot(mu > 0, sigma > 0, r > 0, eps >= 0, eps < 1, delta >= 0, u >= 0, l >= 0)
  EY <- -1/mu
  list(b=b, delta=delta, r=r, eps=eps, sigma=sigma, mu=mu, EY=EY, u=u, l=l)
}

# -------- Heuristic seeds for (xL, gapU) based on model scales ----------
smart_init <- function(p) {
  # "pull" in region 1 (worst-case there): b - δσ + r/μ
  a_pull <- p$b - p$delta*p$sigma + p$r/p$mu
  
  # natural state-length scales: jump length 1/μ, and diffusion length σ/|a_pull|
  s_jump <- 1 / p$mu
  s_diff <- if (abs(a_pull) > 1e-8) p$sigma / abs(a_pull) else 2 / p$mu
  S <- max(s_jump, s_diff)
  
  # drift sign biases which side the lower reflector prefers
  bias <- sign(-(p$b - p$delta*p$sigma))  # negative b and δ>0 => bias > 0 ⇒ start xL < 0
  
  xL0   <- -bias * 0.8 * S
  gapU0 <- 0.9 * S + 0.4 * s_jump
  
  list(xL = xL0, gapU = max(0.2 * s_jump, gapU0))
}

# ----------------------- REGION CONTROLS & COEFFS ------------------------------
region_ambiguity <- function(i, p) {
  if (i == 1)      list(kappa=-p$delta, lambda=p$r*(1+p$eps))
  else if (i == 2) list(kappa=+p$delta, lambda=p$r*(1+p$eps))
  else if (i == 3) list(kappa=+p$delta, lambda=p$r*(1-p$eps))
  else stop("Region must be 1,2,3.")
}

region_a_coeffs <- function(i, p) {
  rc   <- region_ambiguity(i, p)
  a_st <- p$b + p$sigma*rc$kappa + p$r/p$mu
  list(
    a1 = p$mu*a_st - rc$lambda,
    a2 = a_st + 0.5*p$mu*p$sigma^2,
    a3 = 0.5*p$sigma^2,
    a_star = a_st,
    lambda_wc = rc$lambda, kappa_wc = rc$kappa
  )
}

lambda_pm <- function(a_coeffs) {
  disc <- a_coeffs$a2^2 - 4*a_coeffs$a3*a_coeffs$a1
  if (disc <= 0) stop("Non-positive discriminant; adjust parameters.")
  root <- sqrt(disc)
  list(lam_plus = (a_coeffs$a2 + root)/(2*a_coeffs$a3),
       lam_minus= (a_coeffs$a2 - root)/(2*a_coeffs$a3),
       disc=disc)
}

# ----------------------- POLYNOMIAL PART p_i(x) -------------------------------
# Forcing from c(x)=x^2 -> cp(x)=2x; closed-form coefficients for p_i (deg=2)
poly_coeffs <- function(a_coeffs, gamma, p) {
  b2 <- p$mu; b1 <- 1; b0 <- -p$mu*gamma
  if (abs(a_coeffs$a1) > 1e-10) {
    c2 <- - b2 / a_coeffs$a1
    c1 <-  2*(b2*a_coeffs$a2 - a_coeffs$a1) / a_coeffs$a1^2
    c0 <- ( b2*gamma*a_coeffs$a1^2 + 2*a_coeffs$a1*a_coeffs$a2 + 2*b2*a_coeffs$a1*a_coeffs$a3 - 2*b2*a_coeffs$a2^2 ) / a_coeffs$a1^3
    list(c2=c2, c1=c1, c0=c0, deg3=FALSE)
  } else {
    stop("a_coeffs$a1≈0 (rare); extend code for deg-3 fallback if needed.")
  }
}
p_eval  <- function(pc, x) pc$c2*x^2 + pc$c1*x + pc$c0
pp_eval <- function(pc, x) 2*pc$c2*x + pc$c1

# ----------------------- GAMMA ------------------------
gamma_from_xL <- function(xL, p) {
  # Simplified closed form from the lower reflector conditions
  -p$u * (p$b - p$delta*p$sigma) + p$u * p$eps * p$r / p$mu + xL^2
}

# ----------------------- REFLECTOR CONSTANTS c± IN REGIONS 1 & 3 --------------
cpm_from_reflector <- function(xi, lam_pm_i, poly_i, u_or_l, at_lower) {
  # Enforce H(xi)=±u_or_l and H'(xi)=0 to solve for (c-, c+)
  beta1 <- if (at_lower) (-u_or_l - p_eval(poly_i, xi)) else ( u_or_l - p_eval(poly_i, xi))
  beta2 <- pp_eval(poly_i, xi)
  den   <- lam_pm_i$lam_plus - lam_pm_i$lam_minus
  c_minus <- exp(lam_pm_i$lam_minus * xi) * (lam_pm_i$lam_plus * beta1 - beta2) / den
  c_plus  <- exp(lam_pm_i$lam_plus  * xi) * (beta2 - lam_pm_i$lam_minus * beta1) / den
  list(c_minus=c_minus, c_plus=c_plus)
}

make_H_region <- function(lam_pm_i, poly_i, cpm_i) {
  # H
  H <- function(x) {
    cpm_i$c_minus*exp(-lam_pm_i$lam_minus*x) + cpm_i$c_plus*exp(-lam_pm_i$lam_plus*x) + p_eval(poly_i, x)
  }
  # H' (H prime)
  Hp <- function(x) {
    -lam_pm_i$lam_minus*cpm_i$c_minus*exp(-lam_pm_i$lam_minus*x) -
      lam_pm_i$lam_plus *cpm_i$c_plus *exp(-lam_pm_i$lam_plus *x) + pp_eval(poly_i, x)
  }
  list(H=H, Hp=Hp)
}

# ----------------------- H2 CONSTANTS (two variants) --------------------------
# (A) Provisional H2 for Stage 1: match *level & slope* at xk (stable for IH solve)
cpm2_from_xk_C1 <- function(xk, lam_pm2, poly2, H1_xk, H1p_xk) {
  # Solve for (c-, c+) from:
  #   c_- e^{-λ_- xk} + c_+ e^{-λ_+ xk} + p2(xk) = H1(xk) (=0)
  #  -λ_- c_- e^{-λ_- xk} - λ_+ c_+ e^{-λ_+ xk} + p2'(xk) = H1'(xk)
  A <- rbind(
    c(exp(-lam_pm2$lam_minus*xk), exp(-lam_pm2$lam_plus*xk)),
    c(-lam_pm2$lam_minus*exp(-lam_pm2$lam_minus*xk), -lam_pm2$lam_plus*exp(-lam_pm2$lam_plus*xk))
  )
  b <- c(H1_xk - p_eval(poly2, xk), H1p_xk - pp_eval(poly2, xk))
  co <- solve(A, b)
  list(c_minus=co[1], c_plus=co[2])
}

# (B) Final H2 for Stage 2: match *levels* at xk and xlmb (smooth-fit enforced outside)
cpm2_from_levels <- function(xk, xlmb, lam_pm2, poly2, H1_xk, H3_xlmb) {
  A <- rbind(
    c(exp(-lam_pm2$lam_minus*xk),   exp(-lam_pm2$lam_plus*xk)),
    c(exp(-lam_pm2$lam_minus*xlmb), exp(-lam_pm2$lam_plus*xlmb))
  )
  b <- c(H1_xk - p_eval(poly2, xk),
         H3_xlmb - p_eval(poly2, xlmb))
  co <- solve(A, b)
  list(c_minus=co[1], c_plus=co[2])
}

# ----------------------- Closed-form IH at x_lambda (scaled) ------------------
Q_poly <- function(pc, y, mu) {
  pc$c2*( y^2/mu - 2*y/mu^2 + 2/mu^3 ) + pc$c1*( y/mu - 1/mu^2 ) + pc$c0*(1/mu)
}

IH_scaled_at_xlambda <- function(xL, xk, xlmb, LP1, P1, C1, LP2, P2, C2, p) {
  mu <- p$mu
  # guard small denominators
  d <- function(z) { if (abs(z) < 1e-8) (if (z >= 0) 1 else -1)*1e-8 else z }
  
  # (−u) part up to xL
  term0 <- -p$u * exp(mu*(xL - xlmb))
  
  # Region 1 contributions on (xL, xk)
  mm1 <- mu - LP1$lam_minus; mp1 <- mu - LP1$lam_plus
  T1m <- mu * C1$c_minus * exp(-LP1$lam_minus*xlmb) * (exp(mm1*(xk - xlmb)) - exp(mm1*(xL - xlmb))) / d(mm1)
  T1p <- mu * C1$c_plus  * exp(-LP1$lam_plus *xlmb) * (exp(mp1*(xk - xlmb)) - exp(mp1*(xL - xlmb))) / d(mp1)
  
  # Region 2 contributions on (xk, xlmb)
  mm2 <- mu - LP2$lam_minus; mp2 <- mu - LP2$lam_plus
  T2m <- mu * C2$c_minus * exp(-LP2$lam_minus*xlmb) * (1 - exp(mm2*(xk - xlmb))) / d(mm2)
  T2p <- mu * C2$c_plus  * exp(-LP2$lam_plus *xlmb) * (1 - exp(mp2*(xk - xlmb))) / d(mp2)
  
  # Polynomial parts
  Q1_k <- Q_poly(P1, xk, mu);   Q1_L <- Q_poly(P1, xL, mu)
  Q2_l <- Q_poly(P2, xlmb, mu); Q2_k <- Q_poly(P2, xk, mu)
  poly_term <- mu * exp(-mu*xlmb) * ( exp(mu*xk)*Q1_k - exp(mu*xL)*Q1_L + exp(mu*xlmb)*Q2_l - exp(mu*xk)*Q2_k )
  
  # Return the *scaled* residual: multiply by exp(-mu*(xlmb - xk))
  (term0 + T1m + T1p + T2m + T2p + poly_term) * exp(-mu*(xlmb - xk))
}

# ----------------------- ROOT BRACKETING UTILITIES -----------------------------
find_bracket <- function(f, a, step=0.25, max_expand=80, dir=+1) {
  fa <- f(a)
  b  <- a + dir*step; fb <- f(b)
  n  <- 0
  while (fa*fb > 0 && n < max_expand) {
    step <- step*1.6
    a <- b; fa <- fb
    b <- b + dir*step; fb <- f(b)
    n <- n + 1
  }
  if (fa*fb > 0) return(NULL)
  sort(c(a,b))
}

safe_uniroot <- function(f, a, b) {
  uniroot(f, lower=a, upper=b, tol=1e-10, extendInt="yes")$root
}

# ======================= STAGE 1: given xL -> (xk, xlmb) ======================
stage1_inner <- function(xL, p) {
  # 0) gamma first
  gamma <- gamma_from_xL(xL, p)
  
  # Region 1
  A1 <- region_a_coeffs(1, p); LP1 <- lambda_pm(A1); P1 <- poly_coeffs(A1, gamma, p)
  C1 <- cpm_from_reflector(xL, LP1, P1, u_or_l=p$u, at_lower=TRUE)
  H1 <- make_H_region(LP1, P1, C1)
  
  # 1) Solve H1(xk) = 0  (xk > xL)
  f1 <- function(x) H1$H(x)
  br1 <- find_bracket(f1, xL + 1e-6, step=0.05, max_expand=1000, dir=+1)
  if (is.null(br1)) stop("Failed to bracket x_kappa.")
  xk <- safe_uniroot(f1, br1[1], br1[2])
  
  # Region 2 (provisional H2 by C^1 at xk)
  A2 <- region_a_coeffs(2, p); LP2 <- lambda_pm(A2); P2 <- poly_coeffs(A2, gamma, p)
  C2 <- cpm2_from_xk_C1(xk, LP2, P2, H1_xk=H1$H(xk), H1p_xk=H1$Hp(xk))
  H2 <- make_H_region(LP2, P2, C2)
  
  # 2) Solve scaled IH(xlmb) = 0 for xlmb > xk
  f2 <- function(xl) IH_scaled_at_xlambda(
    xL, xk, xl, LP1, P1, C1, LP2, P2, C2, p
  )
  br2 <- find_bracket(f2, xk + 1e-3, step=0.15, max_expand=140, dir=+1)
  if (is.null(br2)) stop("Failed to bracket x_lambda.")
  xlmb <- safe_uniroot(f2, br2[1], br2[2])
  
  list(gamma=gamma,
       A1=A1, LP1=LP1, P1=P1, C1=C1, H1=H1,
       A2=A2, LP2=LP2, P2=P2, C2=C2, H2=H2,
       xk=xk, xlmb=xlmb)
}

# ======================= STAGE 2: solve (xL, xU) ==============================
# Residuals: r3 = H2(xlmb) - H3(lmb); r4 = H2'(xlmb) - H3'(xlmb)
# Scale both by S = 1 + max(|lambda_i^±|) to reduce dynamic range.
outer_residuals <- function(theta, p) {
  # theta: (xL, sU) with xU = xlmb + exp(sU) to enforce xU > xlmb
  xL <- theta[1]
  st1 <- try(stage1_inner(xL, p), silent=TRUE)
  if (inherits(st1, "try-error")) return(c(1e4, 1e4))
  
  xk   <- st1$xk; xlmb <- st1$xlmb
  xU   <- xlmb + exp(theta[2])
  
  # Region 3 from upper reflector at xU
  A3 <- region_a_coeffs(3, p); LP3 <- lambda_pm(A3); P3 <- poly_coeffs(A3, st1$gamma, p)
  C3 <- cpm_from_reflector(xU, LP3, P3, u_or_l=p$l, at_lower=FALSE)
  H3 <- make_H_region(LP3, P3, C3)
  
  # Final H2: match *levels* at xk and xlmb (depends on H3 level at xlmb)
  # C2 <- cpm2_from_levels(xk, xlmb, st1$LP2, st1$P2, st1$H1$H(xk), H3$H(xlmb))
  # H2 <- make_H_region(st1$LP2, st1$P2, C2)
  
  # Residuals (smooth fit) with slope scaling
  S <- 1 + max(abs(c(st1$LP1$lam_minus, st1$LP1$lam_plus,
                     st1$LP2$lam_minus, st1$LP2$lam_plus,
                     LP3$lam_minus, LP3$lam_plus)))
  r3 <- (st1$H2$H(xlmb) - H3$H(xlmb))   #/ S
  r4 <- (st1$H2$Hp(xlmb) - H3$Hp(xlmb)) / S
  c(r3, r4)
}

# Solve the 2x2 outer system (prefer a root finder; fallback to SSR if needed)
solve_two_stage <- function(params,
                            init = list(xL = -1.2, gapU = 1.0),
                            use_nleqslv = TRUE) {
  theta0 <- c(init$xL, log(init$gapU))
  
  if (use_nleqslv && requireNamespace("nleqslv", quietly=TRUE)) {
    f <- function(th) outer_residuals(th, params)
    out <- nleqslv::nleqslv(theta0, f, method="Broyden",
                            control=list(xtol=1e-8, ftol=1e-10, maxit=1000))
    theta <- out$x
    conv  <- out$termcd == 1
    msg   <- out$message
  } else {
    # fallback: minimize SSR
    obj <- function(th) sum(outer_residuals(th, params)^2)
    fit <- optim(theta0, obj, method="BFGS",
                 control=list(reltol=1e-10, maxit=600))
    theta <- fit$par
    conv  <- (fit$convergence == 0)
    msg   <- fit$message
  }
  
  # Reconstruct solution pieces at the optimum
  xL <- theta[1]
  st1 <- stage1_inner(xL, params)
  xk   <- st1$xk; xlmb <- st1$xlmb
  xU   <- xlmb + exp(theta[2])
  
  # Region 3 & final H2
  A3 <- region_a_coeffs(3, params); LP3 <- lambda_pm(A3); P3 <- poly_coeffs(A3, st1$gamma, params)
  C3 <- cpm_from_reflector(xU, LP3, P3, u_or_l=params$l, at_lower=FALSE)
  H3 <- make_H_region(LP3, P3, C3)
  # C2 <- cpm2_from_levels(xk, xlmb, st1$LP2, st1$P2, st1$H1$H(xk), H3$H(xlmb))
  # H2 <- make_H_region(st1$LP2, st1$P2, C2)
  
  # Assemble H (with saturation outside reflectors)
  H_all <- function(x) {
    ifelse(x <= xL, -params$u,
           ifelse(x < xk, st1$H1$H(x),
                  ifelse(x < xlmb, st1$H2$H(x),
                         ifelse(x < xU, H3$H(x), params$l))))
  }
  Hp_all <- function(x) {
    ifelse(x <= xL, 0,
           ifelse(x < xk, st1$H1$Hp(x),
                  ifelse(x < xlmb, st1$H2$Hp(x),
                         ifelse(x < xU, H3$Hp(x), 0))))
  }
  V_all <- function(x) {
    v <- sapply(x, function(xx) {
      if (xx >= 0) integrate(H_all, lower=0, upper=xx, rel.tol=1e-7)$value
      else        -integrate(H_all, lower=xx, upper=0, rel.tol=1e-7)$value
    })
    v
  }
  
  list(
    status = list(converged = conv, message = msg),
    thresholds = list(xL=xL, x_kappa=xk, x_lambda=xlmb, xU=xU),
    gamma = st1$gamma,
    functions = list(H=H_all, Hp=Hp_all, V=V_all),
    diagnostics = list(outer_residuals = outer_residuals(c(xL, log(xU - xlmb)), params))
  )
}

# Build a list of {xL, gapU} initial guesses from vectors (optional helper)
make_inits <- function(xL_vals, gapU_vals) {
  grid <- expand.grid(xL = xL_vals, gapU = gapU_vals, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, ]))
}

# Auto-restart wrapper around your existing solve_two_stage()
solve_two_stage_autorestart <- function(params,
                                        init_list,
                                        use_nleqslv = TRUE,
                                        verbose = TRUE,
                                        accept_nonconverged = FALSE) {
  stopifnot(is.list(init_list), length(init_list) >= 1)
  
  best <- NULL; best_ssr <- Inf; last_err <- NULL
  
  for (k in seq_along(init_list)) {
    init <- init_list[[k]]
    if (verbose) cat(sprintf("[attempt %d/%d] trying init: xL=%.6g, gapU=%.6g\n",
                             k, length(init_list), init$xL, init$gapU))
    
    res <- try(solve_two_stage(params, init = init, use_nleqslv = use_nleqslv),
               silent = TRUE)
    
    if (inherits(res, "try-error")) {
      if (verbose) cat("  -> error; switching to next initial guess.\n")
      last_err <- attr(res, "condition")
      next
    }
    
    ssr <- sum(res$diagnostics$outer_residuals^2)
    
    if (isTRUE(res$status$converged)) {
      if (verbose) cat("  -> converged.\n")
      return(res)
    } else {
      if (verbose) cat("  -> not converged (", res$status$message, ").\n", sep = "")
      if (ssr < best_ssr) { best <- res; best_ssr <- ssr }
      if (accept_nonconverged) {
        if (verbose) cat("  -> returning best-so-far by request.\n")
        return(res)
      }
    }
  }
  
  if (!is.null(best)) {
    if (verbose) cat("No converged solution; returning best-so-far.\n")
    return(best)
  }
  
  msg <- if (!is.null(last_err)) paste("All initial guesses failed; last error:", conditionMessage(last_err))
  else "All initial guesses failed."
  stop(msg)
}

# ===================== Diagnostics function =====================
# Prints and returns:
# - H(x^kappa), (IH)(x^lambda) [scaled & unscaled]
# - C^0/C^1 continuity at xL, x^kappa, x^lambda, xU
# - Boundary conditions at xL (H=-u, H'=0) and xU (H=l, H'=0)
diagnose_solution <- function(sol, params) {
  stopifnot(is.list(sol), !is.null(sol$thresholds))
  xL   <- sol$thresholds$xL
  xk   <- sol$thresholds$x_kappa
  xlmb <- sol$thresholds$x_lambda
  xU   <- sol$thresholds$xU
  
  # --- gamma from xL (your preferred order)
  gamma <- gamma_from_xL(xL, params)
  
  # --- Region 1 (from lower reflector)
  A1  <- region_a_coeffs(1, params);  LP1 <- lambda_pm(A1)
  P1  <- poly_coeffs(A1, gamma, params)
  C1  <- cpm_from_reflector(xL, LP1, P1, u_or_l = params$u, at_lower = TRUE)
  H1  <- make_H_region(LP1, P1, C1)
  
  # --- Region 2 (provisional, C^1 at x^kappa) for IH root evaluation
  A2  <- region_a_coeffs(2, params);  LP2 <- lambda_pm(A2)
  P2  <- poly_coeffs(A2, gamma, params)
  C2 <- cpm2_from_xk_C1(xk, LP2, P2, H1_xk = H1$H(xk), H1p_xk = H1$Hp(xk))
  H2  <- make_H_region(LP2, P2, C2)
  
  # --- (IH)(x^lambda): scaled and unscaled
  IH_scaled <- IH_scaled_at_xlambda(
    xL, xk, xlmb,
    LP1, P1, C1,
    LP2, P2, C2,
    params
  )
  IH_unscaled <- IH_scaled * exp(params$mu * (xlmb - xk))
  
  # --- Region 3 (from upper reflector) and final Region 2 (level-matched)
  A3  <- region_a_coeffs(3, params);  LP3 <- lambda_pm(A3)
  P3  <- poly_coeffs(A3, gamma, params)
  C3  <- cpm_from_reflector(xU, LP3, P3, u_or_l = params$l, at_lower = FALSE)
  H3  <- make_H_region(LP3, P3, C3)
  
  # C2_final <- cpm2_from_levels(xk, xlmb, LP2, P2, H1$H(xk), H3$H(xlmb))
  # H2  <- make_H_region(LP2, P2, C2_final)
  
  # --- Evaluate continuity and boundaries
  diagnostics <- list(
    thresholds = list(xL = xL, x_kappa = xk, x_lambda = xlmb, xU = xU),
    gamma = gamma,
    # Root checks
    H_xk = H1$H(xk),
    IH_scaled_at_xlambda   = IH_scaled,
    IH_unscaled_at_xlambda = IH_unscaled,
    # C^0 / C^1 at xL
    H_xL_plus_u   = H1$H(xL) + params$u,   # should be ~0
    Hp_xL         = H1$Hp(xL),             # should be ~0
    # Continuity at x^kappa
    val_jump_at_xk   = H1$H(xk)  - H2$H(xk),
    slope_jump_at_xk = H1$Hp(xk) - H2$Hp(xk),
    # Continuity at x^lambda
    val_jump_at_xlambda   = H2$H(xlmb)  - H3$H(xlmb),
    slope_jump_at_xlambda = H2$Hp(xlmb) - H3$Hp(xlmb),
    # C^0 / C^1 at xU
    H_xU_minus_l = H3$H(xU) - params$l,   # should be ~0
    Hp_xU        = H3$Hp(xU)              # should be ~0
  )
  
  # --- Pretty print
  fmt <- function(x) sprintf("% .6e", x)
  cat("====== Diagnostics ======\n")
  cat(sprintf("xL = % .6f, x^kappa = % .6f, x^lambda = % .6f, xU = % .6f\n\n",
              xL, xk, xlmb, xU))
  cat("Roots:\n")
  cat("  H(x^kappa)                =", fmt(diagnostics$H_xk), "\n")
  cat("  (IH)(x^lambda) [scaled]   =", fmt(diagnostics$IH_scaled_at_xlambda), "\n")
  cat("  (IH)(x^lambda) [unscaled] =", fmt(diagnostics$IH_unscaled_at_xlambda), "\n\n")
  cat("C^1 & boundary checks (differences from target):\n")
  cat("  At xL:     H(xL)+u        =", fmt(diagnostics$H_xL_plus_u), "\n")
  cat("             H'(xL)         =", fmt(diagnostics$Hp_xL), "\n")
  cat("  At x^kappa: H1-H2         =", fmt(diagnostics$val_jump_at_xk), "\n")
  cat("              H1'-H2'       =", fmt(diagnostics$slope_jump_at_xk), "\n")
  cat("  At x^lambda: H2-H3        =", fmt(diagnostics$val_jump_at_xlambda), "\n")
  cat("               H2'-H3'      =", fmt(diagnostics$slope_jump_at_xlambda), "\n")
  cat("  At xU:     H(xU)-l        =", fmt(diagnostics$H_xU_minus_l), "\n")
  cat("             H'(xU)         =", fmt(diagnostics$Hp_xU), "\n")
  invisible(diagnostics)
}


# ======================= EXAMPLE: your "hard" parameters ======================
# b     <- -1.773496
b     <- 1.99515
delta <- 0.20
r     <- 1.00
eps   <- 0.30
sigma <- 0.60
mu    <- 1.50
u     <- 1.00
l     <- 1.00

params <- make_params(b, delta, r, eps, sigma, mu, u, l)

set.seed(1)
# You control the initial guesses here (order matters; it will try them in sequence)
inits <- list(
  list(xL = -2.5, gapU = 1.0),
  list(xL = -2.0, gapU = 1.05),
  list(xL = -1.6, gapU = 1.1),
  list(xL = -1.2, gapU = 1.15),
  list(xL = -1.0, gapU = 1.2),
  list(xL = -0.5, gapU = 1.25),
  list(xL = -0.2, gapU = 1.3),
  list(xL = -0.15, gapU = 0.9),
  list(xL = -0.05, gapU = 1.4)
)

sol <- solve_two_stage_autorestart(params, init_list = inits, verbose = TRUE)
# sol2 <- solve_two_stage(params, init = smart_init(params))

print(sol$thresholds)
cat("Gamma from xL:", sol$gamma, "\n")
cat("Converged:", sol$status$converged, "-", sol$status$message, "\n")

diag <- diagnose_solution(sol, params)

# Plot (with colored verticals as you requested earlier)
xs   <- unlist(sol$thresholds)[c("xL","x_kappa","x_lambda","xU")]
cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22")
ltys <- c(1,2,2,1)

curve(sol$functions$H(x),
      from = xs[1] - 0.5, to = xs[4] + 0.5, n = 600,
      xlab = "x", ylab = "H(x) = V'(x)")
abline(v = xs, lty = ltys, col = cols, lwd = 2)
abline(h = c(-u,0,l), lty = 3, col = "gray50")
legend("topleft",
       legend = expression(underline(x), x^kappa, x^lambda, bar(x)),
       lty = ltys, col = cols, lwd = 2, bty = "n")
