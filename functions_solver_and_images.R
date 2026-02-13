# ================================== STABLE CORE ==================================
# ============================== Sections 5.1 -- 5.3 ==============================
# Implements the anchored inner solver (u_i^±) and optimal solver with robust numerics.
# Differences vs. previous template:
#   • PER-ROOT anchoring (x_{a,i,−}, x_{a,i,+}) from §5.1: anchor left if λ≥0, right if λ<0.
#   • Region 2: IH(x_λ)=0 uses the anchored + stable form (Eq. (49)) via φ_1mexp.
#   • All H, H' constructions and optimality residuals updated to use per-root anchors.
# References: eqs. (25)–(31), (32), (35), (37)–(41) and their STABLE forms (46)–(51),
#             plus optimality conditions (42)–(45) and stable analogues (52)–(55).
# Paper: "Ergodic singular control for ambiguous compound-Poisson jump diffusion processes"
#        (Sections 4–5; esp. 5.1–5.3).  [Anchoring, stability, and γ updates]
# ---------------------------------------------------------------------------------

# --------------------------------- Utilities -------------------------------------
nz <- function(x, tol=1e-16) abs(x) > tol

# Stable "relative exponential": exprel(z) = (exp(z) - 1)/z with series near z=0
exprel <- function(z) {
  out <- z
  small <- abs(z) < 1e-6
  out[!small] <- expm1(z[!small]) / z[!small]
  # 1 + z/2 + z^2/6  when |z| is tiny (guarding cancellation)
  out[ small] <- 1 + 0.5*z[ small] + (z[ small]^2)/6
  out
}

# Stable φ(a, d) = (1 - exp(-a d)) / a  with series near a=0  (used in Eq. (49))
phi_1mexp <- function(a, d) {
  ad <- a * d
  d * exprel(-ad)
}

# Stable safe division
safe_div <- function(num, den, tol=1e-12, msg="Division by ~0") {
  if (any(abs(den) < tol)) stop(msg)
  num / den
}

# Stable Horner polynomial evaluation for quadratic: c2 x^2 + c1 x + d0
poly2_eval  <- function(c2, c1, d0, x) (c2 * x + c1) * x + d0
poly2_prime <- function(c2, c1, x) 2*c2 * x + c1

# ---------------------------- Parameters & regions -------------------------------
make_params <- function(b, delta, r, eps, sigma, mu, u, l) {
  stopifnot(mu > 0, sigma > 0, r > 0, eps >= 0, eps <= 1, delta >= 0, u >= 0, l >= 0)
  EY <- -1/mu
  list(b=b, delta=delta, r=r, eps=eps, sigma=sigma, mu=mu, EY=EY, u=u, l=l)
}

# Worst-case ambiguity per region (Sec. 4.2, eq. (21))
region_ambiguity <- function(i, p) {
  if (i == 1)      list(kappa=-p$delta, lambda=p$r*(1+p$eps))
  else if (i == 2) list(kappa=+p$delta, lambda=p$r*(1+p$eps))
  else if (i == 3) list(kappa=+p$delta, lambda=p$r*(1-p$eps))
  else stop("Region must be 1,2,3.")
}

# a_{i,·} (Eq. (26))
region_a_coeffs <- function(i, p) {
  rc <- region_ambiguity(i, p)
  a_st <- p$b + p$sigma*rc$kappa + p$r/p$mu      # a* = b + σκ + r/μ
  list(
    a1 = p$mu*a_st - rc$lambda,                   # a_{i,1}
    a2 = a_st + 0.5*p$mu*p$sigma^2,               # a_{i,2}
    a3 = 0.5*p$sigma^2,                           # a_{i,3} = σ^2/2
    a_star = a_st,
    lambda_wc = rc$lambda, kappa_wc = rc$kappa
  )
}

# λ_i^± (28) with robust checks (S2)
lambda_pm <- function(a_coeffs, tol_lambda=1e-12) {
  a1 <- a_coeffs$a1; a2 <- a_coeffs$a2; a3 <- a_coeffs$a3
  disc <- a2^2 - 4*a3*a1
  if (disc <= 0) stop("Discriminant ≤ 0: parameters yield no real roots (check σ>0, λ*>0).")
  root <- sqrt(disc)
  lam_plus  <- (a2 + root)/(2*a3)
  lam_minus <- (a2 - root)/(2*a3)
  if (abs(lam_plus - lam_minus) <= tol_lambda) {
    # (S2) repeated/near-coincident roots: not implemented here
    stop(paste0("S2: near-coincident roots in λ± (|Δ| small). ",
                "Use repeated-root branch or relax tol_lambda; ",
                "current |λ+−λ−|=", signif(abs(lam_plus-lam_minus), 6)))
  }
  list(lam_plus=lam_plus, lam_minus=lam_minus, disc=disc)
}

# ------------------------ q_i(x) coefficients (Sec. 4.1.1) ----------------------
# Assumption 2: a_{i,1} != 0 -> quadratic polynomial; otherwise (S1) abort with message
q_poly_coeffs <- function(a_coeffs, p, tol_a1=1e-12) {
  a1 <- a_coeffs$a1; a2 <- a_coeffs$a2; a3 <- a_coeffs$a3; mu <- p$mu
  if (abs(a1) <= tol_a1) {
    # (S1) cubic fallback (Remark 1) not implemented in anchored γ-free decomposition
    stop("S1: |a_{i,1}|≈0 detected. Use the cubic (Remark 1) branch for p_i(x).")
  }
  # (29) with γ-free shift qi = pi - (μ/a1)γ   (qi independent of γ)
  c2 <- - mu / a1
  c1 <-  2*(mu*a2 - a1)/a1^2
  d0 <- ( 2*a1*a2 + 2*mu*a1*a3 - 2*mu*a2^2 )/a1^3
  list(c2=c2, c1=c1, d0=d0, gamma_fac = mu/a1)
}
q_eval      <- function(qc, x) poly2_eval(qc$c2, qc$c1, qc$d0, x)
qprime_eval <- function(qc, x) poly2_prime(qc$c2, qc$c1, x)

# ----------------------------- Anchoring (Sec. 5.1) -----------------------------
# PER-ROOT anchors: for each root λ_i^±, choose left endpoint if λ≥0, right endpoint if λ<0.
anchors_for_region <- function(lam_minus, lam_plus, x_left, x_right) {
  list(
    minus = if (lam_minus >= 0) x_left else x_right,
    plus  = if (lam_plus  >= 0) x_left else x_right
  )
}

# Helpers: anchored exponentials g(x) = exp(-λ (x - x_a)), g'(x) = -λ * g(x)
g_eval  <- function(lambda, x, xa)  exp(-lambda * (x - xa))
gp_eval <- function(lambda, x, xa) -lambda * exp(-lambda * (x - xa))

# ----------------------------- γ via (46) (stable) ------------------------------
# γ = 0.5 σ^2 H1'(x) + A   with A from §4.2.1; derivative uses PER-ROOT anchors (Eq. (46))
gamma_from_xL_stable <- function(xL, xa1m, xa1p, u1m, u1p, lam1, q1c, p) {
  sigma <- p$sigma; mu <- p$mu; r <- p$r; u <- p$u; delta <- p$delta; eps <- p$eps; b <- p$b
  A <- -u*(b - delta*p$sigma + r/mu) + u*(1+eps)*r/mu + xL^2
  H1p <- gp_eval(lam1$lam_minus, xL, xa1m) * u1m +
    gp_eval(lam1$lam_plus , xL, xa1p) * u1p +
    qprime_eval(q1c, xL)
  0.5*sigma^2 * H1p + A
}

# ---------------------------- Matrix solve (2x2, LU) ----------------------------
solve_2x2_LU <- function(M, b, tag="") {
  # Condition estimates (W2)
  k2 <- tryCatch(kappa(M), error=function(e) Inf)
  if (is.finite(k2) && k2 >= 1e8 && k2 <= 1e12)
    warning(sprintf("W2: cond2(M)≈%.3e for %s (ill-conditioned)", k2, tag))
  if (is.infinite(k2) || k2 > 1e12)
    warning(sprintf("W2 (strong): cond2(M)≈%s for %s", format(k2), tag))
  as.numeric(solve(M, b))
}

# -------------- Region 2 integrals for IH(xλ) (stable, eq. (49)) ---------------
# Qi(y) helper: polynomial contribution used inside (49); depends on q_i and γ via (μ/a_{i,1})γ term
Q_poly <- function(y, qi, ai1, gamma, mu) {
  c2 <- qi$c2; c1 <- qi$c1; d0 <- qi$d0
  term2 <- c2*( y^2/mu - 2*y/(mu^2) + 2/(mu^3) )
  term1 <- c1*( y/mu  - 1/(mu^2) )
  term0 <- ( d0 + (mu/ai1)*gamma )*(1/mu)
  term2 + term1 + term0
}

# ============================== Inner solver (stable) ===========================
# Builds suboptimal H using PER-ROOT anchored u_i^± per Sec. 5.2 (46)–(51)
build_suboptimal_H <- function(p, xL, xk, xl, xU,
                               tol_gap=1e-10, tol_a1=1e-12, tol_lambda=1e-12,
                               kappa3_scale_thresh = 50,
                               tol_regime_switch = 1e-3) {
  # NEW: allow crossing xl ~ xU safely (avoid degenerate Region 3)
  tol_sw <- tol_regime_switch * (1 + abs(xU) + abs(xl))
  regime2 <- (xl > xU)
  if (!regime2 && (xU - xl) <= tol_sw) regime2 <- TRUE
  
  # effective split point for the *piecewise H construction*
  xlam_eff <- if (regime2) xU else xl
  
  
  # Gaps and warnings (S3/W1)
  d1 <- xk - xL
  d2 <- xl - xk
  d3eff <- xU - xlam_eff
  
  if (d1 <= tol_gap) warning("S3/W1: small gap on I1 (xk - xL) → near-collinearity.")
  if (d2 <= tol_gap) warning("S3/W1: small gap on I2 (xl - xk) → near-collinearity.")
  if (!regime2 && d3eff <= tol_gap) warning("S3/W1: small gap on I3 (xU - xl) → near-collinearity.")
  
  
  # Per-region coefficients (a_{i,·}), λ±, and q_i
  a1c <- region_a_coeffs(1, p)
  a2c <- region_a_coeffs(2, p)
  a3c <- region_a_coeffs(3, p)
  
  lam1 <- lambda_pm(a1c, tol_lambda=tol_lambda)
  lam2 <- lambda_pm(a2c, tol_lambda=tol_lambda)
  lam3 <- lambda_pm(a3c, tol_lambda=tol_lambda)
  
  q1c <- q_poly_coeffs(a1c, p, tol_a1)
  q2c <- q_poly_coeffs(a2c, p, tol_a1)
  q3c <- q_poly_coeffs(a3c, p, tol_a1)
  
  # PER-ROOT anchors
  anc1 <- anchors_for_region(lam1$lam_minus, lam1$lam_plus, xL, xk)
  anc2 <- anchors_for_region(lam2$lam_minus, lam2$lam_plus, xk, xlam_eff)
  
  # Region 3 only exists in Regime 1 (xl <= xU)
  if (!regime2) {
    anc3 <- anchors_for_region(lam3$lam_minus, lam3$lam_plus, xl, xU)
  } else {
    anc3 <- list(minus = NA_real_, plus = NA_real_)
  }
  
  # Convenience closures for anchored exponentials
  g1m  <- function(x) g_eval (lam1$lam_minus, x, anc1$minus);  gp1m <- function(x) gp_eval(lam1$lam_minus, x, anc1$minus)
  g1p  <- function(x) g_eval (lam1$lam_plus , x, anc1$plus );  gp1p <- function(x) gp_eval(lam1$lam_plus , x, anc1$plus )
  g2m  <- function(x) g_eval (lam2$lam_minus, x, anc2$minus);  gp2m <- function(x) gp_eval(lam2$lam_minus, x, anc2$minus)
  g2p  <- function(x) g_eval (lam2$lam_plus , x, anc2$plus );  gp2p <- function(x) gp_eval(lam2$lam_plus , x, anc2$plus )
  g3m  <- function(x) g_eval (lam3$lam_minus, x, anc3$minus);  gp3m <- function(x) gp_eval(lam3$lam_minus, x, anc3$minus)
  g3p  <- function(x) g_eval (lam3$lam_plus , x, anc3$plus );  gp3p <- function(x) gp_eval(lam3$lam_plus , x, anc3$plus )
  
  # ------------------------ Region 1: solve u1^± (stable Eq. (47)) ----------------
  u <- p$u; mu <- p$mu; sigma2 <- p$sigma^2; a11 <- a1c$a1
  
  # m rows follow §5.2.2 with derivative term evaluated at xL in both equations
  m11m <- g1m(xL) + (mu/a11)*(sigma2/2)*gp1m(xL)
  m11p <- g1p(xL) + (mu/a11)*(sigma2/2)*gp1p(xL)
  m12m <- g1m(xk) + (mu/a11)*(sigma2/2)*gp1m(xL)
  m12p <- g1p(xk) + (mu/a11)*(sigma2/2)*gp1p(xL)
  
  A   <- -u*(p$b - p$delta*p$sigma + p$r/mu) + u*(1+p$eps)*p$r/mu + xL^2
  b11 <- -u - q_eval(q1c, xL) - (mu/a11)*((sigma2/2)*qprime_eval(q1c, xL) + A)
  b12 <-    0 - q_eval(q1c, xk) - (mu/a11)*((sigma2/2)*qprime_eval(q1c, xL) + A)
  
  M1_raw <- rbind(c(m11m, m11p), c(m12m, m12p))
  b1_raw <- c(b11, b12)
  
  # Condition number before scaling (to decide whether to scale)
  k1_raw <- tryCatch(kappa(M1_raw), error = function(e) Inf)
  
  if (is.finite(k1_raw) && k1_raw > kappa3_scale_thresh) {
    # Simple and effective: row scaling by row norms (or max-abs)
    s1r <- max(1, sqrt(sum(M1_raw[1,]^2)))  # guard with max(1, ·) to avoid over-scaling tiny rows
    s2r <- max(1, sqrt(sum(M1_raw[2,]^2)))
    M1  <- rbind(M1_raw[1,] / s1r, M1_raw[2,] / s2r)
    b1  <- c(b1_raw[1] / s1r, b1_raw[2] / s2r)
    s1  <- solve_2x2_LU(M1, b1, tag = "Region 1 (u1±) [row-scaled]")
  } else {
    M1  <- M1_raw
    b1  <- b1_raw
    s1  <- solve_2x2_LU(M1, b1, tag = "Region 1 (u1±)")
  }
  
  u1m <- s1[1]; u1p <- s1[2]
  
  # γ from (46) using anchored H1′ at xL
  gamma <- gamma_from_xL_stable(xL, anc1$minus, anc1$plus, u1m, u1p, lam1, q1c, p)
  
  # ------------------------ Region 2: solve u2^± ----------------
  if (!regime2) {
    # ===== Regime 1 (xl <= xU): your original IH(xλ)=0 system =====
    
    # Row 1: H-continuity at xκ
    m21m <- g2m(xk)
    m21p <- g2p(xk)
    b21  <- -(mu/a2c$a1)*gamma - q_eval(q2c, xk)
    
    # Row 2: IH(xλ)=0 (anchored, stable form, Eq. (49))
    m22m <- mu * exp(-lam2$lam_minus * (xl - anc2$minus)) * phi_1mexp(mu - lam2$lam_minus, (xl - xk))
    m22p <- mu * exp(-lam2$lam_plus  * (xl - anc2$plus )) * phi_1mexp(mu - lam2$lam_plus , (xl - xk))
    
    # RHS β(2)_2: region 1 homogeneous + polynomials (anchored per root)
    if (abs(mu - lam1$lam_minus) < 10^-6) {
      term_c1_minus <-
        mu * u1m * exp(-mu*(xl - xk)) * exp(-lam1$lam_minus*(xk - anc1$minus)) *
        phi_1mexp(mu - lam1$lam_minus, (xk - xL))
    } else {
      term_c1_minus <-
        mu * u1m * (
          exp(-mu*(xl - xk)) * exp(-lam1$lam_minus*(xk - anc1$minus)) -
            exp(-mu*(xl - xL)) * exp(-lam1$lam_minus*(xL - anc1$minus))
        ) / (mu - lam1$lam_minus)
    }
    
    if (abs(mu - lam1$lam_plus) < 10^-6) {
      term_c1_plus <-
        mu * u1p * exp(-mu*(xl - xk)) * exp(-lam1$lam_plus*(xk - anc1$plus)) *
        phi_1mexp(mu - lam1$lam_plus, (xk - xL))
    } else {
      term_c1_plus <-
        mu * u1p * (
          exp(-mu*(xl - xk)) * exp(-lam1$lam_plus*(xk - anc1$plus)) -
            exp(-mu*(xl - xL)) * exp(-lam1$lam_plus*(xL - anc1$plus))
        ) / (mu - lam1$lam_plus)
    }
    
    term_c1 <- term_c1_minus + term_c1_plus
    
    Q1_xk <- Q_poly(xk, q1c, a1c$a1, gamma, mu)
    Q1_xL <- Q_poly(xL, q1c, a1c$a1, gamma, mu)
    Q2_xl <- Q_poly(xl, q2c, a2c$a1, gamma, mu)
    Q2_xk <- Q_poly(xk, q2c, a2c$a1, gamma, mu)
    
    term_poly <- mu*exp(-mu*xl) * (exp(mu*xk)*Q1_xk - exp(mu*xL)*Q1_xL + exp(mu*xl)*Q2_xl - exp(mu*xk)*Q2_xk)
    b22 <- + u*exp(mu*(xL - xl)) - term_c1 - term_poly
    
    M2_raw <- rbind(c(m21m, m21p), c(m22m, m22p))
    b2_raw <- c(b21, b22)
    
    k2_raw <- tryCatch(kappa(M2_raw), error = function(e) Inf)
    if (is.finite(k2_raw) && k2_raw > kappa3_scale_thresh) {
      s1r <- max(1, sqrt(sum(M2_raw[1,]^2)))
      s2r <- max(1, sqrt(sum(M2_raw[2,]^2)))
      M2  <- rbind(M2_raw[1,] / s1r, M2_raw[2,] / s2r)
      b2  <- c(b2_raw[1] / s1r, b2_raw[2] / s2r)
      s2  <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [row-scaled]")
    } else {
      M2 <- M2_raw
      b2 <- b2_raw
      s2 <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±)")
    }
    
    u2m <- s2[1]; u2p <- s2[2]
    
  } else {
    # ===== Regime 2 (xl > xU): NO region 3 inside the band. Solve u2± by H2(xk)=0 and H2(xU)=l =====
    
    # Row 1: H2(xk)=0
    m21m <- g2m(xk)
    m21p <- g2p(xk)
    b21  <- -(mu/a2c$a1)*gamma - q_eval(q2c, xk)
    
    # Row 2: H2(xU)=l  (since H is constant = l for x>=xU)
    m22m <- g2m(xU)
    m22p <- g2p(xU)
    b22  <- p$l - (mu/a2c$a1)*gamma - q_eval(q2c, xU)
    
    M2_raw <- rbind(c(m21m, m21p), c(m22m, m22p))
    b2_raw <- c(b21, b22)
    
    k2_raw <- tryCatch(kappa(M2_raw), error = function(e) Inf)
    if (is.finite(k2_raw) && k2_raw > kappa3_scale_thresh) {
      s1r <- max(1, sqrt(sum(M2_raw[1,]^2)))
      s2r <- max(1, sqrt(sum(M2_raw[2,]^2)))
      M2  <- rbind(M2_raw[1,] / s1r, M2_raw[2,] / s2r)
      b2  <- c(b2_raw[1] / s1r, b2_raw[2] / s2r)
      s2  <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [Regime 2, row-scaled]")
    } else {
      M2 <- M2_raw
      b2 <- b2_raw
      s2 <- solve_2x2_LU(M2, b2, tag = "Region 2 (u2±) [Regime 2]")
    }
    
    u2m <- s2[1]; u2p <- s2[2]
  }
  
  
  # ------------------------ Region 3: solve u3^± -------------------
  if (!regime2) {
    # (your original Region 3 block stays EXACTLY the same)
    
    # Row 1: H3(xU) = l
    m31m <- g3m(xU)
    m31p <- g3p(xU)
    b31  <- p$l - (mu/a3c$a1)*gamma - q_eval(q3c, xU)
    
    # Row 2: continuity at xλ (use anchors xa3 for region 3, xa2 for region 2)
    m32m <- g3m(xl)
    m32p <- g3p(xl)
    b32  <- u2m * g2m(xl) + u2p * g2p(xl) +
      mu*gamma*(1/a2c$a1 - 1/a3c$a1) + q_eval(q2c, xl) - q_eval(q3c, xl)
    
    M3_raw <- rbind(c(m31m, m31p), c(m32m, m32p))
    b3_raw <- c(b31, b32)
    
    # Condition number before scaling (to decide whether to scale)
    k3_raw <- tryCatch(kappa(M3_raw), error = function(e) Inf)
    
    if (is.finite(k3_raw) && k3_raw > kappa3_scale_thresh) {
      # Simple and effective: row scaling by row norms (or max-abs)
      s1 <- max(1, sqrt(sum(M3_raw[1,]^2)))  # guard with max(1, ·) to avoid over-scaling tiny rows
      s2 <- max(1, sqrt(sum(M3_raw[2,]^2)))
      M3  <- rbind(M3_raw[1,] / s1, M3_raw[2,] / s2)
      b3  <- c(b3_raw[1] / s1, b3_raw[2] / s2)
      s3  <- solve_2x2_LU(M3, b3, tag = "Region 3 (u3±) [row-scaled]")
      scaled_row_factors <- c(s1, s2)
    } else {
      M3  <- M3_raw
      b3  <- b3_raw
      s3  <- solve_2x2_LU(M3, b3, tag = "Region 3 (u3±)")
      scaled_row_factors <- c(1, 1)
    }
    
    u3m <- s3[1]; u3p <- s3[2]
  
  } else {
    # Regime 2: no Region 3 inside [xL, xU]. We keep placeholders for compatibility.
    u3m <- 0; u3p <- 0
    M3  <- matrix(NA_real_, 2, 2)
  }
  
  # ----------------------------- Build H and H' -----------------------------------
  H1_fun  <- function(x) u1m*g1m(x) + u1p*g1p(x) + q_eval(q1c, x) + q1c$gamma_fac*gamma
  H2_fun  <- function(x) u2m*g2m(x) + u2p*g2p(x) + q_eval(q2c, x) + q2c$gamma_fac*gamma
  
  if (!regime2) {
    H3_fun  <- function(x) u3m*g3m(x) + u3p*g3p(x) + q_eval(q3c, x) + q3c$gamma_fac*gamma
    Hp3_fun <- function(x) gp3m(x)*u3m + gp3p(x)*u3p + qprime_eval(q3c, x)
  } else {
    H3_fun  <- function(x) rep(p$l, length(x))
    Hp3_fun <- function(x) rep(0, length(x))
  }
  
  Hp1_fun <- function(x) gp1m(x)*u1m + gp1p(x)*u1p + qprime_eval(q1c, x)
  Hp2_fun <- function(x) gp2m(x)*u2m + gp2p(x)*u2p + qprime_eval(q2c, x)
  
  
  H <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    left  <- which(x <= xL)
    mid1  <- which(x >  xL & x <  xk)
    mid2  <- which(x >= xk & x <  xlam_eff)
    mid3  <- which(x >= xlam_eff & x <  xU)
    right <- which(x >= xU)
    if (length(left))  out[left]  <- -p$u
    if (length(mid1))  out[mid1]  <- H1_fun(x[mid1])
    if (length(mid2))  out[mid2]  <- H2_fun(x[mid2])
    if (length(mid3))  out[mid3]  <- H3_fun(x[mid3])
    if (length(right)) out[right] <-  p$l
    out
  }
  
  Hp <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    left  <- which(x <= xL)
    mid1  <- which(x >  xL & x <  xk)
    mid2  <- which(x >= xk & x <  xlam_eff)
    mid3  <- which(x >= xlam_eff & x <  xU)
    right <- which(x >= xU)
    if (length(left))  out[left]  <- 0
    if (length(mid1))  out[mid1]  <- Hp1_fun(x[mid1])
    if (length(mid2))  out[mid2]  <- Hp2_fun(x[mid2])
    if (length(mid3))  out[mid3]  <- Hp3_fun(x[mid3])
    if (length(right)) out[right] <- 0
    out
  }
  
  list(
    H=H, Hp=Hp, gamma=gamma,
    regime2 = regime2,
    pieces=list(
      region1=list(u_minus=u1m, u_plus=u1p, a=a1c, lam=lam1, q=q1c,
                   anchor_minus=anc1$minus, anchor_plus=anc1$plus, M=M1),
      region2=list(u_minus=u2m, u_plus=u2p, a=a2c, lam=lam2, q=q2c,
                   anchor_minus=anc2$minus, anchor_plus=anc2$plus, M=M2),
      region3=list(
        u_minus=u3m, u_plus=u3p, a=a3c, lam=lam3, q=q3c, 
        anchor_minus=anc3$minus, anchor_plus=anc3$plus, M=M3
      )
    ),
    params=p, x=list(xL=xL, xk=xk, xl=xl, xU=xU)
  )
}

# ------------------------ IH(xλ) residual (stable, eq. (49)) --------------------
# ------------------------ IH(xλ) residual (stable, generalized) ----------------
# NEW: works both when xλ<=xU (Regime 1) and when xλ>xU (Regime 2).
# In Regime 2, it correctly adds the constant-tail contribution from [xU, xλ].
IH_at_xlambda_stable <- function(sol, xla = NULL) {
  if (is.null(xla)) xla <- sol$x$xl
  
  p  <- sol$params
  mu <- p$mu
  xs <- sol$x
  xL <- xs$xL; xk <- xs$xk; xU <- xs$xU
  
  # last point that is inside the band (integration upper limit for the interior pieces)
  B <- min(xla, xU)
  
  R1 <- sol$pieces$region1
  R2 <- sol$pieces$region2
  
  q1c <- R1$q; q2c <- R2$q
  a1c <- R1$a; a2c <- R2$a
  lam1 <- R1$lam; lam2 <- R2$lam
  xa1m <- R1$anchor_minus; xa1p <- R1$anchor_plus
  xa2m <- R2$anchor_minus; xa2p <- R2$anchor_plus
  u1m <- R1$u_minus; u1p <- R1$u_plus
  u2m <- R2$u_minus; u2p <- R2$u_plus
  gamma <- sol$gamma
  
  # ----- Region 1 homogeneous contribution over [xL, xk], evaluated at xla -----
  if (abs(mu - lam1$lam_minus) < 1e-6) {
    term_c1_minus <-
      mu * u1m * exp(-mu*(xla - xk)) * exp(-lam1$lam_minus*(xk - xa1m)) *
      phi_1mexp(mu - lam1$lam_minus, (xk - xL))
  } else {
    term_c1_minus <-
      mu * u1m * (
        exp(-mu*(xla - xk)) * exp(-lam1$lam_minus*(xk - xa1m)) -
          exp(-mu*(xla - xL)) * exp(-lam1$lam_minus*(xL - xa1m))
      ) / (mu - lam1$lam_minus)
  }
  
  if (abs(mu - lam1$lam_plus) < 1e-6) {
    term_c1_plus <-
      mu * u1p * exp(-mu*(xla - xk)) * exp(-lam1$lam_plus*(xk - xa1p)) *
      phi_1mexp(mu - lam1$lam_plus, (xk - xL))
  } else {
    term_c1_plus <-
      mu * u1p * (
        exp(-mu*(xla - xk)) * exp(-lam1$lam_plus*(xk - xa1p)) -
          exp(-mu*(xla - xL)) * exp(-lam1$lam_plus*(xL - xa1p))
      ) / (mu - lam1$lam_plus)
  }
  
  term_c1 <- term_c1_minus + term_c1_plus
  
  # ----- Region 2 homogeneous contribution over [xk, B], evaluated at xla -----
  # NOTE: extra factor exp(-mu*(xla-B)) appears if xla>B (i.e. Regime 2)
  shift_B <- exp(-mu * (xla - B))
  
  term_c2 <-
    shift_B * (
      mu * u2m * exp(-lam2$lam_minus*(B - xa2m)) * phi_1mexp(mu - lam2$lam_minus, (B - xk)) +
        mu * u2p * exp(-lam2$lam_plus *(B - xa2p)) * phi_1mexp(mu - lam2$lam_plus , (B - xk))
    )
  
  # ----- Polynomial contributions over [xL, xk] and [xk, B] -----
  Q1_xk <- Q_poly(xk, q1c, a1c$a1, gamma, mu)
  Q1_xL <- Q_poly(xL, q1c, a1c$a1, gamma, mu)
  Q2_B  <- Q_poly(B , q2c, a2c$a1, gamma, mu)
  Q2_xk <- Q_poly(xk, q2c, a2c$a1, gamma, mu)
  
  term_poly <-
    mu * exp(-mu*xla) *
    (exp(mu*xk)*Q1_xk - exp(mu*xL)*Q1_xL + exp(mu*B)*Q2_B - exp(mu*xk)*Q2_xk)
  
  # ----- Constant tail when xla > xU: contribution from [xU, xla] where H=l -----
  tail <- if (xla > xU) p$l * (1 - exp(mu*(xU - xla))) else 0
  
  # ----- Full residual (should be 0 at the true xλ) -----
  -p$u * exp(mu*(xL - xla)) + term_c1 + term_c2 + term_poly + tail
}

# ================================= Diagnostics ==================================
diagnose <- function(sol, tol=1e-8) {
  p  <- sol$params; xs <- sol$x
  xL <- xs$xL; xk <- xs$xk; xl <- xs$xl; xU <- xs$xU
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2; R3 <- sol$pieces$region3
  gamma <- sol$gamma
  
  # Closures with PER-ROOT anchors
  g  <- function(lambda, x, xa)  exp(-lambda * (x - xa))
  gp <- function(lambda, x, xa) -lambda * exp(-lambda * (x - xa))
  
  H1_fun <- function(x) R1$u_minus*g(R1$lam$lam_minus,x,R1$anchor_minus) + R1$u_plus*g(R1$lam$lam_plus,x,R1$anchor_plus) +
    q_eval(R1$q, x) + R1$q$gamma_fac*gamma
  H2_fun <- function(x) R2$u_minus*g(R2$lam$lam_minus,x,R2$anchor_minus) + R2$u_plus*g(R2$lam$lam_plus,x,R2$anchor_plus) +
    q_eval(R2$q, x) + R2$q$gamma_fac*gamma
  H3_fun <- function(x) R3$u_minus*g(R3$lam$lam_minus,x,R3$anchor_minus) + R3$u_plus*g(R3$lam$lam_plus,x,R3$anchor_plus) +
    q_eval(R3$q, x) + R3$q$gamma_fac*gamma
  
  H1p_at <- function(x) gp(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus + gp(R1$lam$lam_plus,x,R1$anchor_plus)*R1$u_plus + qprime_eval(R1$q, x)
  H2p_at <- function(x) gp(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus + gp(R2$lam$lam_plus,x,R2$anchor_plus)*R2$u_plus + qprime_eval(R2$q, x)
  H3p_at <- function(x) gp(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus + gp(R3$lam$lam_plus,x,R3$anchor_plus)*R3$u_plus + qprime_eval(R3$q, x)
  
  cont_xL <- H1_fun(xL) + p$u
  cont_xk <- H1_fun(xk) - H2_fun(xk)
  cont_xl <- H2_fun(xl) - H3_fun(xl)
  cont_xU <- H3_fun(xU) - p$l
  
  root_k <- H1_fun(xk)               # H(xκ)=0
  IH_xl  <- IH_at_xlambda_stable(sol) # IH(xλ)=0
  
  dcont_xL <- H1p_at(xL) - 0
  dcont_xk <- H2p_at(xk) - H1p_at(xk)
  dcont_xl <- H3p_at(xl) - H2p_at(xl)
  dcont_xU <- 0 - H3p_at(xU)
  
  checks <- data.frame(
    check = c("continuity@xL", "continuity@xk", "continuity@xl", "continuity@xU",
              "H(xk)=0", "IH(xl)=0",
              "H'(jump)@xL", "H'(jump)@xk", "H'(jump)@xl", "H'(jump)@xU"),
    value = c(cont_xL, cont_xk, cont_xl, cont_xU, root_k, IH_xl,
              dcont_xL, dcont_xk, dcont_xl, dcont_xU),
    type  = c(rep("continuity", 4), rep("ambiguity_condition", 2), rep("derivative_jump", 4)),
    pass  = abs(c(cont_xL, cont_xk, cont_xl, cont_xU, root_k, IH_xl,
                  dcont_xL, dcont_xk, dcont_xl, dcont_xU)) <= tol
  )
  rownames(checks) <- NULL
  
  kappa2 <- function(M) {
    tryCatch(kappa(M), error=function(e) Inf)
  }
  kappas <- data.frame(
    block   = c("M1 (Region1)", "M2 (Region2)", "M3 (Region3)"),
    kappa2  = c(kappa2(R1$M), kappa2(R2$M), kappa2(R3$M))
  )
  list(checks=checks, condition_numbers=kappas)
}

# =========================== Optimal barriers (stable) ===========================
# NEW: allow xl and xU to be on either side, as long as both are > xk.
z_from_x <- function(xL, xk, xl, xU) {
  stopifnot(xL < xk, xk < xl, xk < xU)
  c(
    xL,
    log(xk - xL),   # > 0
    log(xl - xk),   # > 0
    log(xU - xk)    # > 0
  )
}

x_from_z <- function(z) {
  stopifnot(length(z) == 4)
  xL <- z[1]
  d1 <- exp(z[2])     # xk - xL
  dL <- exp(z[3])     # xl - xk
  dU <- exp(z[4])     # xU - xk
  xk <- xL + d1
  xl <- xk + dL
  xU <- xk + dU
  c(xL = xL, xk = xk, xl = xl, xU = xU)
}

# Residuals of optimality equations with Regime-2 support (xl may exceed xU)
opt_conditions_residuals <- function(z, p, tol_build=1e-10, tol_regime_switch = 1e-3) {
  xs  <- x_from_z(z)
  sol <- build_suboptimal_H(p, xs["xL"], xs["xk"], xs["xl"], xs["xU"],
                            tol_regime_switch = tol_regime_switch)
  
  xL <- xs["xL"]; xk <- xs["xk"]; xl <- xs["xl"]; xU <- xs["xU"]
  R1 <- sol$pieces$region1; R2 <- sol$pieces$region2
  
  H1p <- function(x) gp_eval(R1$lam$lam_minus,x,R1$anchor_minus)*R1$u_minus +
    gp_eval(R1$lam$lam_plus ,x,R1$anchor_plus )*R1$u_plus  +
    qprime_eval(R1$q, x)
  
  H2p <- function(x) gp_eval(R2$lam$lam_minus,x,R2$anchor_minus)*R2$u_minus +
    gp_eval(R2$lam$lam_plus ,x,R2$anchor_plus )*R2$u_plus  +
    qprime_eval(R2$q, x)
  
  if (!isTRUE(sol$regime2)) {
    # Regime 1: your original 4 smooth-fit equations (needs Region 3)
    R3 <- sol$pieces$region3
    H3p <- function(x) gp_eval(R3$lam$lam_minus,x,R3$anchor_minus)*R3$u_minus +
      gp_eval(R3$lam$lam_plus ,x,R3$anchor_plus )*R3$u_plus  +
      qprime_eval(R3$q, x)
    
    r1 <- H1p(xL)                 # (52)
    r2 <- H3p(xU)                 # (53)
    r3 <- H1p(xk) - H2p(xk)       # (54)
    r4 <- H2p(xl) - H3p(xl)       # (55)
    c(r1, r2, r3, r4)
  } else {
    # Regime 2: no Region 3 inside band. Use:
    #   (i) smooth-fit at xL, xU, xk
    #   (ii) definition of xλ via IH(xλ)=0 (now xλ may be > xU)
    r1 <- H1p(xL)
    r2 <- H2p(xU)
    r3 <- H1p(xk) - H2p(xk)
    r4 <- IH_at_xlambda_stable(sol, xla = xl)
    c(r1, r2, r3, r4)
  }
}


# ------------------------ Outer solver: Broyden → Newton ------------------------
# Switch to a pure Newton method (with FD Jacobian) when close to a solution.
# Nearness is detected by max|f| <= switch_ftol (and after at least switch_iter_min iterations).
# Uses the stable inner builder from Secs. 5.1–5.3.

solve_optimal_barriers <- function(
    p, xL0, xk0, xl0, xU0,
    method = "Broyden",
    control = list(ftol=1e-11, xtol=1e-11, maxit=1e6, stepmax=1, trace=0),
    tol_build = 1e-10,
    verbose = FALSE,
    tol_regime_switch = 1e-3,
    # --- hybrid controls ---
    switch_to_newton = FALSE,
    switch_ftol = 5e-2,          # threshold to trigger Newton (on max|f|)
    switch_iter_min = 2L,        # don't switch on the very first step
    newton_control = list(ftol=1e-12, xtol=1e-12, maxit=1e5, stepmax=1, trace=0),
    fd_step = 1e-6               # relative step for finite-difference Jacobian
) {
  if (!requireNamespace("nleqslv", quietly=TRUE))
    stop("Please install.packages('nleqslv')")
  if (!(xL0 < xk0 && xk0 < xl0 && xk0 < xU0))
    stop("Initial guess must satisfy xL0 < xk0, and xk0 < xl0, xk0 < xU0 (xl0 may be > xU0).")
  
  # ---- reparam z <-> x (keeps barriers ordered) ----
  z0 <- z_from_x(xL0, xk0, xl0, xU0)
  
  # Residuals (Eqs. (52)–(55) in stable form)
  fn <- function(z) {
    r <- opt_conditions_residuals(z, p, tol_build=tol_build, tol_regime_switch = tol_regime_switch)
    if (isTRUE(verbose)) {
      xs <- x_from_z(z)
      cat(sprintf(
        "xL=%.6g  xk=%.6g  xl=%.6g  xU=%.6g  |  r=[% .3e, % .3e, % .3e, % .3e]\n",
        xs["xL"], xs["xk"], xs["xl"], xs["xU"], r[1], r[2], r[3], r[4]
      ))
    }
    r
  }
  
  # Stable, central-difference Jacobian (O(h^2))
  jac_fd <- function(z, fz = NULL) {
    n <- length(z); J <- matrix(0.0, n, n)
    if (is.null(fz)) fz <- fn(z)
    for (j in 1:n) {
      h <- fd_step * (1 + abs(z[j]))
      zp <- z; zm <- z
      zp[j] <- z[j] + h;  zm[j] <- z[j] - h
      fp <- fn(zp);       fm <- fn(zm)
      J[, j] <- (fp - fm) / (2*h)
    }
    J
  }
  
  # -------------------- Phase 1: robust Broyden --------------------
  ans_broyden <- nleqslv::nleqslv(
    x = z0, fn = fn, method = method,
    global = "dbldog", xscalm = "auto", control = control
  )
  
  # Decide whether to switch to Newton
  do_switch <- isTRUE(switch_to_newton)
  if (do_switch) {
    fmax   <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
    iters  <- ans_broyden$iter %||% 0L
    do_switch <- is.finite(fmax) && (fmax <= switch_ftol) && (iters >= switch_iter_min)
  }
  
  # -------------------- Phase 2: pure Newton (optional) --------------------
  if (do_switch) {
    z_start <- ans_broyden$x
    jf <- function(z) jac_fd(z)  # supply Jacobian explicitly
    ans_newton <- tryCatch(
      nleqslv::nleqslv(
        x = z_start, fn = fn, jac = jf, method = "Newton",
        global = "none", xscalm = "auto", control = newton_control
      ),
      error = function(e) e
    )
    
    # If Newton succeeds and improves max|f|, take it; otherwise keep Broyden result
    if (!inherits(ans_newton, "error")) {
      f_b <- max(abs(ans_broyden$fvec %||% fn(ans_broyden$x)))
      f_n <- max(abs(ans_newton$fvec  %||% fn(ans_newton$x)))
      final_ans <- if (is.finite(f_n) && (f_n <= f_b)) ans_newton else ans_broyden
      chosen <- if (identical(final_ans, ans_newton)) "newton" else "broyden"
    } else {
      final_ans <- ans_broyden
      chosen <- "broyden"
      warning(sprintf("Newton phase skipped (error: %s). Keeping Broyden solution.", conditionMessage(ans_newton)))
    }
  } else {
    final_ans <- ans_broyden
    chosen <- "broyden"
  }
  
  # -------------------- Build final suboptimal H and return --------------------
  xs_hat <- x_from_z(final_ans$x)
  sol_subopt <- build_suboptimal_H(
    p, xs_hat[["xL"]], xs_hat[["xk"]], xs_hat[["xl"]], xs_hat[["xU"]]
  )
  
  list(
    H      = sol_subopt$H,
    Hp     = sol_subopt$Hp,
    gamma  = sol_subopt$gamma,
    pieces = sol_subopt$pieces,
    params = sol_subopt$params,
    x      = sol_subopt$x,
    nleqslv_primary = ans_broyden,
    nleqslv_refined = if (exists("ans_newton")) ans_newton else NULL,
    chosen_solver   = chosen
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ================================ Figures =======================================
.with_plot_margins <- function(expr) { 
  op <- par(mar=c(2,2,0.5,0.5)+0.2)
  on.exit(par(op), add=TRUE); force(expr)
}

render_and_save <- function(
    fname_base, plotfun, show=interactive(), save=TRUE,
    dir="figures",
    width_in = 6.67, height_in = 4.67, dpi = 300,
    pointsize = NULL, family = NULL,
    match_current = TRUE, use_cairo_pdf = TRUE
) {
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  # Derive device settings (size, ps, family) from current device when available
  if (match_current && dev.cur() != 1L) {
    sz <- dev.size("in")
    width_in  <- sz[1]; height_in <- sz[2]
    if (is.null(pointsize)) pointsize <- par("ps")
    if (is.null(family))    family    <- par("family")
  }
  if (is.null(pointsize)) pointsize <- 12
  fam <- if (is.null(family) || !nzchar(family)) "sans" else family
  
  # SHOW: draw on the current device but force the same family for consistency
  if (isTRUE(show)) {
    op <- par(family = fam)               # temporary — reverts on exit
    on.exit(par(op), add = TRUE)
    .with_plot_margins(plotfun())
  }
  
  if (!isTRUE(save)) return(invisible())
  
  # PNG: exact physical size + pointsize
  png(file.path(dir, paste0(fname_base, ".png")),
      width = width_in, height = height_in, units = "in",
      res = dpi, pointsize = pointsize,
      type = getOption("bitmapType", "cairo"), antialias = "subpixel")
  .with_plot_margins(plotfun()); dev.off()
  
  # PDF: same physical size + pointsize + family
  if (isTRUE(use_cairo_pdf) && capabilities("cairo")) {
    cairo_pdf(file.path(dir, paste0(fname_base, ".pdf")),
              width = width_in, height = height_in,
              pointsize = pointsize, family = fam)
  } else {
    pdf(file.path(dir, paste0(fname_base, ".pdf")),
        width = width_in, height = height_in,
        pointsize = pointsize, family = fam, useDingbats = FALSE)
  }
  .with_plot_margins(plotfun()); dev.off()
}

plot_H <- function(sol, params, show=interactive(), save=TRUE, n = 600,
                   pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                   out_base="H_thresholds", dir="figures",
                   margins = c(3.2, 3.25, 1.05, 0.9),
                   axis_mgp = c(2.2, 0.7, 0),
                   tcl = -0.25,
                   show_axis_titles = TRUE,
                   show_tick_labels = TRUE,
                   tick_cex = 1.7,
                   base_cex = 1, mex = 1) {
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22"); ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(sol$H(seq(x_from, x_to, l = n)))
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(xaxs="i",
        oma = c(0,0,0,0),
        cex = base_cex,
        mex = mex,
        mar = margins,
        mgp = axis_mgp,
        tcl = tcl,
        cex.axis = tick_cex)
    
    curve(sol$H(x), from=x_from, to=x_to, n=n,
          xlab="",
          ylab=if (show_axis_titles) expression(H(x) == V*minute*(x)) else "",
          col="black", lwd=1.4, ylim=ylim_y,
          xaxt="n", yaxt="n", bty="n")
    
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    
    segments(x0=xs, y0=rep(-params$u, length(xs)), x1=xs, y1=rep(params$l, length(xs)),
             lty=ltys, col=cols, lwd=2)
    abline(h=-params$u, lty=3, col="#777777")
    abline(h=0,       lty=3, col="#BBBBBB")
    abline(h=params$l, lty=4, col="#777777")
    
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 5), 0.38 * inch_to_user)
    
    legend("topleft",
           legend=expression(H(x), underline(x), x^kappa, x^lambda, bar(x), -u, l),
           lty=c(1,1,2,2,1,3,4),
           col=c("black", cols, "#777777", "#777777"),
           lwd=c(1.4, rep(2,6)),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}


plot_H_prime <- function(sol, params, show=interactive(), save=TRUE, n = 800,
                         pad_x_frac = 0.10, top_blank = 0.10, bottom_blank = 0.05,
                         out_base="Hprime_thresholds", dir="figures") {
  stopifnot(!is.null(sol$Hp))
  
  xs   <- unlist(sol$x)[c("xL","xk","xl","xU")]
  xL <- xs[1]; xk <- xs[2]; xl <- xs[3]; xU <- xs[4]
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22")  # xL,xk,xl,xU
  ltys <- c(1,2,2,1)
  
  x_max <- max(xU, xl)
  W <- x_max - xL
  x_from <- xL - pad_x_frac*W
  x_to   <- x_max + pad_x_frac*W
  
  # sample grid and break the line at thresholds to avoid artificial connections
  xg <- seq(x_from, x_to, length.out = n)
  yg <- sol$Hp(xg)
  
  # Force line breaks around the 4 special points (xL,xk,xl,xU)
  for (b in xs) {
    j <- which.min(abs(xg - b))
    for (k in c(j-1, j, j+1)) {
      if (k >= 1 && k <= length(yg)) yg[k] <- NA_real_
    }
  }
  
  # y-limits with controlled blank margins (like plot_H)
  frac_band <- 1 - top_blank - bottom_blank
  ylim_y <- range(yg, finite = TRUE)
  if (!all(is.finite(ylim_y))) ylim_y <- c(-1, 1)
  if (diff(ylim_y) == 0) {
    bump <- max(1, abs(ylim_y[1]))
    ylim_y <- ylim_y + c(-1, 1) * 0.05 * bump
  }
  ylim_y <- ylim_y + diff(ylim_y)/frac_band * c(-bottom_blank, top_blank)
  
  plotfun <- function() {
    op <- par(xaxs="i"); on.exit(par(op), add=TRUE)
    
    plot(xg, yg, type="l",
         xlab="", ylab=expression(H*minute*(x) == W(x)),
         col="black", lwd=1.4, ylim=ylim_y)
    
    # thresholds as vertical segments spanning the panel
    segments(x0=xs, y0=rep(ylim_y[1], length(xs)),
             x1=xs, y1=rep(ylim_y[2], length(xs)),
             lty=ltys, col=cols, lwd=2)
    
    # reference line at 0
    abline(h=0, lty=3, col="#BBBBBB")
    
    # legend (same style as plot_H)
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user,
            rep(col_w, 5),
            0.20 * inch_to_user)
    
    legend("topleft",
           legend = expression(H*minute*(x), underline(x), x^kappa, x^lambda, bar(x), 0),
           lty    = c(1, 1, 2, 2, 1, 3),
           col    = c("black", cols, "#BBBBBB"),
           lwd    = c(1.4, rep(2,4), 1.2),
           bty="n", horiz=TRUE, seg.len=1.5, x.intersp=0.5,
           text.width=tw)
    
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}

# --------------------------- Reflected process sim ------------------------------
simulate_reflected_jd <- function(T=8, dt=0.001, x0=NULL, params, thresholds,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(round(T/dt)); tt <- seq(0, T, length.out=n+1)
  xL <- thresholds$xL; xU <- thresholds$xU
  if (is.null(x0)) x0 <- (thresholds$xL + thresholds$xU)/2
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

# --------------------------- Small reusable drawer for the middle panel --------
# Draw the reflected trajectory with thresholds; returns the ylim used (so panels
# can be tuned consistently if needed).
.draw_reflected_panel <- function(t, sim, xL, xU, xk, xl,
                                  top_blank=0.075, bottom_blank=0.05,
                                  draw_legend=TRUE) {
  Hcorr <- xU - xL
  frac  <- 1 - top_blank - bottom_blank
  Htot  <- Hcorr/frac
  ylim  <- c(xL - Htot*bottom_blank, xU + Htot*top_blank)
  
  cols <- c("#B22222","#1E90FF","#8A2BE2","#228B22") # xL, xk, xl, xU
  ltys <- c(1,2,2,1)
  
  plot(t, sim, type="n", xlab="", ylab="", ylim=ylim, xaxt="n")
  usr <- par("usr")
  rect(usr[1], xL, usr[2], xU, col=adjustcolor("gray85", 0.6), border=NA)
  
  lines(t, sim, lwd=1.2)
  abline(h=c(xL, xU), col=c(cols[1], cols[4]), lty=c(ltys[1], ltys[4]), lwd=2)
  abline(h=c(xk, xl), col=c(cols[2], cols[3]), lty=c(ltys[2], ltys[3]), lwd=2)
  
  if (draw_legend) {
    inch_to_user <- diff(grconvertX(c(0,1), from="in", to="user"))
    col_w <- 0.23 * inch_to_user
    tw <- c(0.36 * inch_to_user, rep(col_w, 4))  # 5 entries total
    
    legend("topleft",
           legend = expression(bar(X)[t], underline(x), x^kappa, x^lambda, bar(x)),
           lty    = c(1, 1, 2, 2, 1),                # 5 entries
           col    = c("black", cols),                # 1 + 4 entries
           lwd    = c(1.2, rep(2,4)),
           bty    = "n", horiz = TRUE, x.intersp = 0.5, seg.len = 2,
           text.width = tw
    )
  }
  box()
  invisible(ylim)
}

# --------------------------- Original single plots (kept, just call helper) ----
plot_reflected_jd <- function(sol, params, sim=NULL, seed=123, show=interactive(), 
                              save=TRUE, out_base="reflected_jd", dir="figures",
                              top_blank=0.075, bottom_blank=0.05) {
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  plotfun <- function() {
    .draw_reflected_panel(t, X, xL, xU, xk, xl, top_blank=top_blank, bottom_blank=bottom_blank)
    axis(1); mtext("time", side=1, line=2)
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}

plot_controls <- function(sim, show=interactive(), save=TRUE,
                          out_base="singular_controls", dir="figures") {
  plotfun <- function() {
    rng <- range(sim$U, sim$L)
    plot(sim$time, sim$L, type="s", xlab="time", ylab="cumulative push",
         lwd=1.6, ylim=rng, lty=2)
    lines(sim$time, sim$U, type="s", lwd=1.6)
    legend("topleft",
           legend=c(expression(L[t]~"(pushes down)"), expression(U[t]~"(pushes up)")),
           lty=c(2,1), lwd=2, bty="n")
    box()
  }
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}

# --------------------------- NEW: merged 3-panel figure -------------------------
# Top:   L_t (pushes down)
# Middle: reflected process + thresholds
# Bottom: U_t (pushes up)
plot_reflected_with_controls <- function(sol, params, sim=NULL, seed=123,
                                         show=interactive(), save=TRUE,
                                         out_base="reflected_with_controls",
                                         dir="figures", draw_legend = TRUE,
                                         heights = c(0.7, 1.7, 0.7),
                                         top_blank=0.075, bottom_blank=0.05,
                                         margins = c(3.2, 3.25, 1.05, 0.9),
                                         axis_mgp = c(2.2, 0.7, 0),
                                         tcl = -0.25,
                                         show_axis_titles = TRUE,
                                         show_tick_labels = TRUE,
                                         tick_cex = 1.7,
                                         base_cex = 1, mex = 1) { 
  xs <- sol$x
  if (is.null(sim)) { set.seed(seed); sim <- simulate_reflected_jd(params=params, thresholds=xs) }
  t <- sim$time; X <- sim$X; xL <- xs$xL; xU <- xs$xU; xk <- xs$xk; xl <- xs$xl
  
  bottom <- margins[1]; left <- margins[2]; top <- margins[3]; right <- margins[4]
  
  plotfun <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add=TRUE)
    
    par(oma = c(0,0,0,0), cex = base_cex, mex = mex)
    
    layout(matrix(1:3, nrow=3), heights = heights)
    
    # --- TOP: L_t
    par(mar = c(0.4, left, top, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex)
    plot(t, sim$L, type="s",
         xlab="", ylab=if (show_axis_titles) expression(L[t]) else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#228B22")
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(L[t]), cex = 1.5,
           lty = 1, lwd = 2, col = "#228B22", bty = "n")
    box()
    
    # --- MIDDLE
    par(mar = c(0.2, left, 0.2, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex)
    .draw_reflected_panel(t, X, xL, xU, xk, xl,
                          top_blank=top_blank, bottom_blank=bottom_blank,
                          draw_legend=draw_legend)
    
    # --- BOTTOM: U_t
    par(mar = c(bottom, left, 0.4, right),
        mgp=axis_mgp, tcl=tcl,
        cex = base_cex, mex = mex,
        cex.axis=tick_cex)
    plot(t, sim$U, type="s",
         xlab=if (show_axis_titles) "time" else "",
         ylab=if (show_axis_titles) expression(U[t]) else "",
         xaxt="n", yaxt="n",
         lwd=2, col = "#B22222")
    axis(1, labels = show_tick_labels)
    axis(2, labels = show_tick_labels)
    legend("topleft", legend = expression(U[t]), cex = 1.7,
           lty = 1, lwd = 2, col = "#B22222", bty = "n")
    box()
  }
  
  render_and_save(out_base, plotfun, show=show, save=save, dir=dir)
}


# ============================= Comparative sweeps ================================
comparative_sweeper <- function(
    sweep_param = "b",
    sweep_values = seq(-10, 10, by = 0.1),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1.0, l = 1.0,
    xL0 = -0.5, xk0 = -0.1, xl0 = 0.3, xU0 = 1.0,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    save = FALSE
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sweep_param <- as.character(sweep_param)
  # NEW: allow "inv_mu" as synthetic parameter = 1/mu
  known_params <- c("b","delta","r","eps","sigma","mu","inv_mu","u","l")
  if (!sweep_param %in% known_params) {
    stop(sprintf("sweep_param must be one of: %s", paste(known_params, collapse = ", ")))
  }
  fixed_params <- list(b = b, delta = delta, r = r, eps = eps, sigma = sigma, mu = mu, u = u, l = l)
  
  # Metrics: thresholds + ergodic value γ
  res <- data.frame(
    sweep_param = sweep_param,
    sweep_value = as.numeric(sweep_values),
    xL = NA_real_, xk = NA_real_, xl = NA_real_, xU = NA_real_,
    gamma = NA_real_,
    converged = FALSE, error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  last_sol <- NULL
  
  for (i in seq_along(sweep_values)) {
    val <- sweep_values[i]
    message(sprintf("[sweep %s] %2d/%d  %s = % .6g ...",
                    sweep_param, i, length(sweep_values), sweep_param, val)); flush.console()
    
    p <- fixed_params
    
    # NEW: map 1/mu ("inv_mu") to actual mu > 0
    if (sweep_param == "inv_mu") {
      if (val <= 0) {
        res$converged[i] <- FALSE
        res$error_msg[i] <- "inv_mu (1/mu) must be > 0"
        next
      }
      p$mu <- 1/val
    } else {
      p[[sweep_param]] <- val
    }
    
    params <- do.call(make_params, p)
    
    # warm-start from last solution if available
    xL_init <- xL0; xk_init <- xk0; xl_init <- xl0; xU_init <- xU0
    if (!is.null(last_sol)) {
      xL_init <- last_sol$x$xL; xk_init <- last_sol$x$xk
      xl_init <- last_sol$x$xl; xU_init <- last_sol$x$xU
    }
    
    sol_opt <- tryCatch(
      solve_optimal_barriers(params, xL_init, xk_init, xl_init, xU_init, verbose = FALSE),
      error = function(e) e
    )
    
    if (inherits(sol_opt, "error")) {
      res$converged[i] <- FALSE
      res$error_msg[i] <- conditionMessage(sol_opt)
    } else {
      th <- sol_opt$x
      res[i, c("xL","xk","xl","xU")] <- unlist(th[c("xL","xk","xl","xU")])
      res$gamma[i] <- sol_opt$gamma
      res$converged[i] <- TRUE
      res$error_msg[i] <- NA_character_
      last_sol <- sol_opt
    }
  }
  
  # Default name (used only for CSV; plotting can choose its own if desired)
  if (is.null(out_name)) {
    step_str <- if (length(sweep_values) > 1) {
      diffs <- unique(round(diff(sweep_values), 10))
      if (length(diffs) == 1) sprintf("by_%g", diffs) else "custom_steps"
    } else "single_value"
    out_name <- sprintf("thresholds_vs_%s_%g_to_%g_%s",
                        sweep_param, min(sweep_values), max(sweep_values), step_str)
  }
  
  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # NEW: build a filename suffix with all NON-swept parameter values
    fixed_for_name <- fixed_params
    if (sweep_param %in% names(fixed_for_name)) fixed_for_name[[sweep_param]] <- NULL
    if (sweep_param == "inv_mu") fixed_for_name$mu <- NULL  # mu is being varied via inv_mu
    
    fixed_tag <- paste(
      sprintf("%s_%s",
              names(fixed_for_name),
              vapply(fixed_for_name, function(x) sprintf("%.6g", x), "")
              
      ),
      collapse = "__"
    )
    
    csv_path <- file.path(out_dir, paste0(out_name, "__", fixed_tag, ".csv"))
    write.csv(res, csv_path, row.names = FALSE)
    message(sprintf("Saved sweep metrics to: %s", normalizePath(csv_path)))
  }
  
  # Return metrics object (no plotting)
  list(
    results     = res,
    out_dir     = out_dir,
    out_name    = out_name,
    fixed_params = fixed_params
  )
}


plot_sweep <- function(
    sweep_obj,
    out_dir  = file.path("figures", "sweeps"),
    out_name = NULL,
    show = interactive(),
    save = FALSE,
    cols = c("#B22222", "#1E90FF", "#8A2BE2", "#228B22"),
    ltys = c(1, 1, 1, 1),
    title = TRUE, axis_labs = TRUE,
    plot_gamma = FALSE, 
    gamma_layout = c("stacked","separate")
) {
  gamma_layout <- match.arg(gamma_layout)
  
  # Accept either the list returned by comparative_sweeper() or a bare data.frame
  res <- if (is.data.frame(sweep_obj)) {
    sweep_obj
  } else {
    sweep_obj$results %||% sweep_obj
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ok <- res$converged
  if (!any(ok)) {
    warning("No successful solutions; nothing to plot.")
    return(invisible(NULL))
  }
  
  sweep_param <- as.character(unique(res$sweep_param))
  if (length(sweep_param) != 1L) {
    warning("Multiple sweep_param values detected; using the first.")
    sweep_param <- sweep_param[1L]
  }
  
  sweep_values <- res$sweep_value
  
  # Default fname base if not provided (use raw sweep_param, not pretty label)
  if (is.null(out_name)) {
    step_str <- if (length(sweep_values) > 1) {
      diffs <- unique(round(diff(sweep_values), 10))
      if (length(diffs) == 1) sprintf("by_%g", diffs) else "custom_steps"
    } else "single_value"
    out_name <- sprintf("thresholds_vs_%s_%g_to_%g_%s",
                        sweep_param, min(sweep_values), max(sweep_values), step_str)
  }
  
  # append non-swept parameter values to filename (only if available)
  fixed_for_name <- NULL
  if (!is.data.frame(sweep_obj) && !is.null(sweep_obj$fixed_params)) {
    fixed_for_name <- sweep_obj$fixed_params
  }
  
  if (!is.null(fixed_for_name)) {
    if (sweep_param %in% names(fixed_for_name)) fixed_for_name[[sweep_param]] <- NULL
    if (sweep_param == "inv_mu") fixed_for_name$mu <- NULL
    
    fixed_tag <- paste(
      sprintf("%s_%s",
              names(fixed_for_name),
              vapply(fixed_for_name, function(x) sprintf("%.6g", x), "")
      ),
      collapse = "__"
    )
    fixed_tag <- gsub("\\s+", "", fixed_tag)
    out_name <- paste0(out_name, "__", fixed_tag)
  }
  
  # nice label for titles (e.g. "1/mu" instead of "inv_mu")
  param_label_str <- if (sweep_param == "inv_mu") "1/mu" else sweep_param
  
  # y-range for thresholds
  ylm_th <- range(as.numeric(unlist(res[ok, c("xL","xk","xl","xU")])), na.rm = TRUE)
  
  # y-range for gamma (if present)
  has_gamma <- "gamma" %in% names(res) && any(is.finite(res$gamma[ok]))
  if (!has_gamma) plot_gamma <- FALSE
  if (plot_gamma) {
    ylm_g <- range(res$gamma[ok], na.rm = TRUE)
  }
  
  xlab_expr <- if (axis_labs) {
    switch(sweep_param,
           "b"      = expression(b),
           "delta"  = expression(delta),
           "r"      = expression(r),
           "eps"    = expression(epsilon),
           "sigma"  = expression(sigma),
           "mu"     = expression(mu),
           "inv_mu" = expression(1/mu),
           "u"      = expression(u),
           "l"      = expression(l),
           sweep_param)
  } else ""
  
  title_str   <- if (title)
    sprintf("Ambiguity thresholds & barriers vs %s", param_label_str) else ""
  title_gamma <- if (title)
    sprintf("Ergodic value %s vs %s", "\u03b3", param_label_str) else ""
  
  ## ---------- Helper to draw thresholds-only panel ----------
  thresholds_panel <- function(xlab_bottom = TRUE, mar_override = NULL) {
    mar <- if (is.null(mar_override)) c(2.2, 2.2, 1.5, 1.2) else mar_override
    op <- par(xaxs = "i", mar = mar, cex.axis=1.4); on.exit(par(op), add = TRUE)
    
    xlab <- if (xlab_bottom) xlab_expr else ""
    plot(res$sweep_value, res$xL, type = "n", xlab = xlab,
         ylab = ifelse(axis_labs, "threshold value", ""), ylim = ylm_th)
    grid()
    lines(res$sweep_value, res$xL, lwd = 2, lty = ltys[1], col = cols[1])
    lines(res$sweep_value, res$xk, lwd = 2, lty = ltys[2], col = cols[2])
    lines(res$sweep_value, res$xl, lwd = 2, lty = ltys[3], col = cols[3])
    lines(res$sweep_value, res$xU, lwd = 2, lty = ltys[4], col = cols[4])
    legend("topleft",
           legend = expression(underline(x), x^kappa, x^lambda, bar(x)),
           col = cols, lty = ltys, lwd = 2, bty = "n", horiz = TRUE,
           x.intersp = 0.5, seg.len = 2)
    title(title_str)
  }
  
  ## ---------- Helper to draw gamma-only panel ----------
  gamma_panel <- function(xlab_bottom = TRUE, mar_override = NULL) {
    mar  <- if (is.null(mar_override)) c(2.2, 2.2, 1.5, 1.2) else mar_override
    op   <- par(xaxs = "i", mar = mar, cex.axis=1.4); on.exit(par(op), add = TRUE)
    
    xlab <- if (xlab_bottom) xlab_expr else ""
    xaxt <- if (xlab_bottom) "s" else "n"   # hide x-axis when xlab_bottom = FALSE
    
    plot(res$sweep_value, res$gamma, type = "l", lwd = 2,
         xlab = xlab, ylab = expression(gamma),
         ylim = ylm_g, xaxt = xaxt)
    grid()
    legend("topleft", legend = expression(gamma), lwd = 2, bty = "n")
    title(title_gamma)
  }
  
  ## ======================= CASE 1: thresholds only =======================
  if (!plot_gamma) {
    plotfun <- function() {
      thresholds_panel()
    }
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name, plotfun = plotfun,
                      show = show, save = save, dir = out_dir)
    } else {
      if (isTRUE(show)) plotfun()
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1200, height = 800, res = 150); plotfun(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 9, height = 6); plotfun(); dev.off()
      }
    }
    
    ## ======================= CASE 2: gamma + thresholds, separate =======================
  } else if (gamma_layout == "separate") {
    plotfun_th <- function() thresholds_panel()
    plotfun_g  <- function() gamma_panel()
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name,
                      plotfun = plotfun_th, show = show, save = save, dir = out_dir)
      render_and_save(fname_base = paste0(out_name, "_gamma"),
                      plotfun = plotfun_g, show = show, save = save, dir = out_dir)
    } else {
      if (isTRUE(show)) {
        plotfun_th()
        plotfun_g()
      }
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1200, height = 800, res = 150); plotfun_th(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 9, height = 6); plotfun_th(); dev.off()
        png(file.path(out_dir, paste0(out_name, "_gamma.png")),
            width = 1200, height = 800, res = 150); plotfun_g(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, "_gamma.pdf")),
            width = 9, height = 6); plotfun_g(); dev.off()
      }
    }
    
    ## ======================= CASE 3: gamma + thresholds, stacked 2×1 =======================
  } else {  # gamma_layout == "stacked"
    plotfun <- function() {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      layout(matrix(1:2, nrow = 2), heights = c(0.35, 0.65))  # (% gamma, % thresholds)
      
      # Top: gamma (shared x, no axis)
      gamma_panel(xlab_bottom = FALSE, mar_override = c(0.4, 3, 1.05, 0.9))
      
      # Bottom: thresholds (x-axis shown)
      thresholds_panel(xlab_bottom = TRUE, mar_override = c(3, 3, 0.4, 0.9))
    }
    
    if (exists("render_and_save")) {
      render_and_save(fname_base = out_name, plotfun = plotfun,
                      show = show, save = save, dir = out_dir)
    } else {
      if (isTRUE(show)) plotfun()
      if (isTRUE(save)) {
        png(file.path(out_dir, paste0(out_name, ".png")),
            width = 1600, height = 1600, res = 150); plotfun(); dev.off()
        pdf(file.path(out_dir, paste0(out_name, ".pdf")),
            width = 12, height = 12); plotfun(); dev.off()
      }
    }
  }
  
  invisible(NULL)
}

fit_best_family <- function(
    x, y,
    families  = c("monomial", "exponential"),
    criterion = c("AIC", "RSS"),
    verbose   = TRUE,
    max_tries = 6              # how many jittered tries for each nonlinear family
) {
  criterion <- match.arg(criterion)
  families  <- intersect(families, c("monomial", "exponential"))
  if (length(families) == 0L) stop("No valid families given.")
  
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!all(ok)) {
    if (verbose) message("Dropping ", sum(!ok), " non-finite (x,y) pairs.")
    x <- x[ok]; y <- y[ok]
  }
  if (length(x) < 3L) stop("Need at least 3 data points.")
  
  df <- data.frame(x = x, y = y)
  xrange <- range(df$x)
  xspan  <- diff(xrange)
  if (xspan <= 0) stop("x must not be constant.")
  
  # -------- safe building blocks to avoid NaN/Inf -----------------
  safe_pow <- function(x, b, c, eps = 1e-8, max_power = 50) {
    # base >= eps so no negative or zero to weird powers
    base <- pmax(x - b, eps)
    val  <- base^c
    val[!is.finite(val)] <- max_power
    pmin(val, max_power)
  }
  
  safe_exp_pow <- function(x, b, c, eps = 1e-8, max_exponent = 700) {
    base <- pmax(x - b, eps)
    power <- base^c
    power[!is.finite(power)] <- max_exponent
    exp(pmin(power, max_exponent))
  }
  
  fits <- list()
  crit <- data.frame(
    family = character(),
    RSS    = numeric(),
    AIC    = numeric(),
    stringsAsFactors = FALSE
  )
  
  use_nlsLM <- requireNamespace("minpack.lm", quietly = TRUE)
  if (verbose && use_nlsLM) {
    message("Using minpack.lm::nlsLM for nonlinear fits.")
  }
  
  # ---------- helper: robust nonlinear fitter with jittered starts ----------
  fit_nonlinear_family <- function(formula, start, lower, upper, family_label) {
    best_fit <- NULL
    best_RSS <- Inf
    last_err <- NULL
    
    for (k in seq_len(max_tries)) {
      st <- start
      
      # Jitter after the first attempt
      if (k > 1) {
        for (nm in names(st)) {
          v <- st[[nm]]
          if (!is.finite(v)) next
          j  <- rnorm(1, mean = 0, sd = 0.3)
          st[[nm]] <- v * (1 + j)
        }
        # enforce bounds
        for (nm in names(st)) {
          if (!is.null(lower) && nm %in% names(lower)) {
            st[[nm]] <- max(st[[nm]], lower[[nm]])
          }
          if (!is.null(upper) && nm %in% names(upper)) {
            st[[nm]] <- min(st[[nm]], upper[[nm]])
          }
        }
      }
      
      fit_try <- try(
        if (use_nlsLM) {
          minpack.lm::nlsLM(
            formula, data = df, start = st,
            lower = lower, upper = upper,
            control = minpack.lm::nls.lm.control(maxiter = 1000)
          )
        } else {
          nls(
            formula, data = df, start = st,
            algorithm = "port",
            lower = lower, upper = upper,
            control = list(maxiter = 500)
          )
        },
        silent = TRUE
      )
      
      if (!inherits(fit_try, "try-error")) {
        RSS <- sum(resid(fit_try)^2)
        if (RSS < best_RSS) {
          best_RSS <- RSS
          best_fit <- fit_try
        }
      } else {
        last_err <- fit_try
      }
    }
    
    if (is.null(best_fit) && verbose && !is.null(last_err)) {
      cond <- attr(last_err, "condition")
      msg  <- if (inherits(cond, "condition")) conditionMessage(cond) else as.character(last_err)
      message("  ", family_label, " fit failed after ", max_tries, " tries: ", msg)
    }
    best_fit
  }
  
  ## ------------- 2) Monomial: y = d + a * safe_pow(x, b, c) -------------
  if ("monomial" %in% families) {
    if (verbose) message("Fitting monomial model y = d + a*safe_pow(x, b, c) ...")
    
    # base starting values
    b0 <- xrange[1] - 0.1 * xspan      # somewhere to the left
    t0 <- safe_pow(df$x, b0, 1)        # essentially x-b0
    d0 <- min(df$y)
    
    # crude a0, c0 via log–log if y-d0>0 and t0>0
    if (all(t0 > 0)) {
      y_shift <- df$y - d0
      if (all(y_shift > 0)) {
        llm <- try(lm(log(y_shift) ~ log(t0), data = df), silent = TRUE)
        if (!inherits(llm, "try-error")) {
          c0  <- coef(llm)[2]
          a0  <- exp(coef(llm)[1])
        } else {
          a0 <- mean(y_shift)
          c0 <- 1
        }
      } else {
        a0 <- mean(df$y)
        c0 <- 1
      }
    } else {
      a0 <- mean(df$y)
      c0 <- 1
    }
    
    # bounds: c in [0, 5] to avoid wild curvature; b within a generous x-range; d free
    lower_m <- c(a = -Inf, b = xrange[1] - 5*xspan, c = -5,   d = -Inf)
    upper_m <- c(a =  Inf, b = xrange[2] + 5*xspan, c = 5, d =  Inf)
    start_m <- list(a = a0, b = b0, c = c0, d = d0)
    
    m_mono <- fit_nonlinear_family(
      y ~ d + a * safe_pow(x, b, c),
      start = start_m,
      lower = lower_m,
      upper = upper_m,
      family_label = "Monomial"
    )
    
    if (!is.null(m_mono)) {
      RSS_m <- sum(resid(m_mono)^2)
      AIC_m <- AIC(m_mono)
      fits$monomial <- m_mono
      crit <- rbind(crit,
                    data.frame(family = "monomial", RSS = RSS_m, AIC = AIC_m))
    }
  }
  
  ## --- 3) Exponential: y = d + a * safe_exp_pow(x, b, c) (or decay version) ---
  if ("exponential" %in% families) {
    if (verbose) message("Fitting exponential model y = d + a*safe_exp_pow(x, b, c) ...")
    
    b0 <- xrange[1] - 0.1 * xspan
    d0 <- min(df$y)
    a0 <- max(df$y) - d0
    if (!is.finite(a0) || a0 == 0) a0 <- mean(df$y) - d0
    c0 <- 1
    
    lower_e <- c(a = -Inf, b = xrange[1] - 5*xspan, c = -Inf,   d = -Inf)
    upper_e <- c(a =  Inf, b = xrange[2] + 5*xspan, c = 5.0, d =  Inf)
    start_e <- list(a = a0, b = b0, c = c0, d = d0)
    
    m_exp <- fit_nonlinear_family(
      # If you want explicit decay, switch the model to:
      # y ~ d + a * (1 / safe_exp_pow(x, b, c))
      y ~ d + a * safe_exp_pow(x, b, c),
      start = start_e,
      lower = lower_e,
      upper = upper_e,
      family_label = "Exponential"
    )
    
    if (!is.null(m_exp)) {
      RSS_e <- sum(resid(m_exp)^2)
      AIC_e <- AIC(m_exp)
      fits$exponential <- m_exp
      crit <- rbind(crit,
                    data.frame(family = "exponential", RSS = RSS_e, AIC = AIC_e))
    }
  }
  
  if (nrow(crit) == 0L) stop("No model could be fitted.")
  
  ## ------------------------- Choose best model -------------------------
  best_row <- if (criterion == "AIC") which.min(crit$AIC) else which.min(crit$RSS)
  best_family <- crit$family[best_row]
  best_fit    <- fits[[best_family]]
  
  if (verbose) {
    message(sprintf("Best family: %s (criterion = %s)", best_family, criterion))
    print(crit)
  }
  
  list(
    best_family = best_family,
    best_fit    = best_fit,
    all_fits    = fits,
    criteria    = crit
  )
}




# # ---------------------- EXAMPLE -----------------------------------------------
# Parameters
p <- make_params(b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1, l = 1)
xL0 <- -0.5; xk0 <- -0.1; xl0 <- 0.5; xU0 <- 1

# Suboptimal H (inner solver)
subopt_H <- build_suboptimal_H(p, xL0, xk0, xl0, xU0)
diagnose(subopt_H)
# Plot H
plot_H(subopt_H, p, show = TRUE, save = FALSE)
plot_H_prime(subopt_H, p, show = TRUE, save = FALSE)

# Optimal H (outer solver)
sol_opt <- solve_optimal_barriers(p, xL0, xk0, xl0, xU0, verbose=TRUE, 
                                  switch_to_newton = FALSE, tol_regime_switch = 1e-3)
diagnose(sol_opt)
# Plot H and simulate the reflected process
plot_H(sol_opt, p, show = TRUE, save = TRUE, 
       top_blank = 0.075, bottom_blank = 0.05, 
       margins = c(2, 2, 1, 1), show_axis_titles = FALSE, tick_cex = 1)
plot_H_prime(sol_opt, p, show = TRUE, save = FALSE)
sim <- simulate_reflected_jd(params=p, thresholds=sol_opt$x, seed = 123)
plot_reflected_jd(sol_opt, p, sim=sim, show=TRUE, save=TRUE)
plot_controls(sim, show=TRUE, save=TRUE)
plot_reflected_with_controls(sol_opt, p, sim, top_blank = 0, bottom_blank = 0, 
                             heights = c(0.7, 1.7, 0.85), draw_legend = FALSE,
                             margins = c(2, 2, 1, 1), axis_mgp = c(0.3, 0.7, 0),
                             show_axis_titles = FALSE, tick_cex = 1,
                             show=TRUE, save=TRUE)

# ---------------------- EXAMPLE SWEEPER ---------------------------------------
cost_matrix <- matrix(c(1, 1,
                        2, 1,
                        1, 2), nrow = 3, byrow = T)

## 1) Sweep b
sweep_b <- comparative_sweeper(
  sweep_param  = "b",
  sweep_values = seq(-10, 10, by = 0.01),
  delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1.5, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_b,
  title = FALSE, axis_labs = FALSE,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE, 
)

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 1) Sweep b
  sweep_b <- comparative_sweeper(
    sweep_param  = "b",
    sweep_values = seq(-20, 20, by = 0.01),
    delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_b,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 2) Sweep delta
  sweep_delta <- comparative_sweeper(
    sweep_param  = "delta",
    sweep_values = seq(0, 100, by = 0.01),
    b = 0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_delta,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_delta$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 3) Sweep r
  sweep_r <- comparative_sweeper(
    sweep_param  = "r",
    sweep_values = seq(0.05, 100, by = 0.05),
    b = 0, delta = 1.0, eps = 0.5, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_r,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_r$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 4) Sweep eps (epsilon)
  sweep_eps <- comparative_sweeper(
    sweep_param  = "eps",
    sweep_values = seq(0, 1, by = 0.001),
    b = 0, delta = 1.0, r = 1, sigma = 1, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_eps,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_eps$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 5) Sweep sigma
  sweep_sigma <- comparative_sweeper(
    sweep_param  = "sigma",
    sweep_values = seq(0.1, 50, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, mu = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_sigma,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_sigma$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 6) Sweep mu
  sweep_mu <- comparative_sweeper(
    sweep_param  = "mu",
    sweep_values = 1/seq(0.51, 5, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_mu,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_mu$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

for (i in 1:nrow(cost_matrix)) {
  
  cat("u = ", cost_matrix[i,1], ", l = ", "u = ", cost_matrix[i,2], ". \n", sep = "")
  
  ## 6.1) Sweep 1/mu (E[Y])
  sweep_inv_mu <- comparative_sweeper(
    sweep_param  = "inv_mu",
    sweep_values = seq(0, 3, by = 0.01),
    b = 0, delta = 1.0, r = 1, eps = 0.5, sigma = 1, u = cost_matrix[i,1], l = cost_matrix[i,2],,,
    save = TRUE
  )
  plot_sweep(
    sweep_obj = sweep_inv_mu,
    title = FALSE, axis_labs = FALSE,
    plot_gamma = TRUE, gamma_layout = "stacked",
    save = TRUE, 
  )
  # Fitting from family of functions
  res <- sweep_inv_mu$results
  ok  <- res$converged
  x <- res$sweep_value[ok]
  y <- res$gamma[ok]
  fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
  fit_obj$best_family
  coef(fit_obj$best_fit)
  sum(residuals(fit_obj$best_fit)^2)
  # Make a smooth curve over the same range
  yy <- fit_obj$best_fit$m$fitted()
  plot(x, y, pch = 16)
  lines(x, yy, lwd = 2, col = "red")
  
}

## 7) Sweep u
# b = 3
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 5, by = 0.01),
  b = 3, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_u,
  title = FALSE, axis_labs = FALSE,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE, 
)
# b = -3
sweep_u <- comparative_sweeper(
  sweep_param  = "u",
  sweep_values = seq(0.01, 5, by = 0.01),
  b = -3, delta = 1, r = 1, eps = 0.5, sigma = 1, mu = 1, l = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_u,
  title = FALSE, axis_labs = FALSE,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE, 
)
# Fitting from family of functions
res <- sweep_u$results
ok  <- res$converged
x <- res$sweep_value[ok]
y <- res$gamma[ok]
fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
fit_obj$best_family
coef(fit_obj$best_fit)
sum(residuals(fit_obj$best_fit)^2)
# Make a smooth curve over the same range
yy <- fit_obj$best_fit$m$fitted()
plot(x, y, pch = 16)
lines(x, yy, lwd = 2, col = "red")

## 8) Sweep l
# b = 3
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 5, by = 0.01),
  b = 3, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_l,
  title = FALSE, axis_labs = FALSE,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE, 
)
# b = -3
sweep_l <- comparative_sweeper(
  sweep_param  = "l",
  sweep_values = seq(1, 5, by = 0.01),
  b = -3, delta = 1.0, r = 1, eps = 0.5, sigma = 1, mu = 1, u = 1,
  save = TRUE
)
plot_sweep(
  sweep_obj = sweep_l,
  title = FALSE, axis_labs = FALSE,
  plot_gamma = TRUE, gamma_layout = "stacked",
  save = TRUE, 
)
# Fitting from family of functions
res <- sweep_l$results
ok  <- res$converged
x <- res$sweep_value[ok]
y <- res$gamma[ok]
fit_obj <- fit_best_family(x, y, criterion = "RSS", verbose = TRUE)
fit_obj$best_family
coef(fit_obj$best_fit)
sum(residuals(fit_obj$best_fit)^2)
# Make a smooth curve over the same range
yy <- fit_obj$best_fit$m$fitted()
plot(x, y, pch = 16)
lines(x, yy, lwd = 2, col = "red")

# ---------------------- RESIDUALS ---------------------------------------------
# Parameters
p  <- make_params(b=100, delta=0.4, r=1.0, eps=0.2, sigma=0.6, mu=1.5, u=1.0, l=1.2)
xL0 <- -1; xk0 <- -0.5; xl0 <- 0.5; xU0 <- 1.0
base_x <- list(xL=xL0, xk=xk0, xl=xl0, xU=xU0)

# Threshold ranges (respect ordering gaps)
xL_seq <- seq(xk0 - 20, xk0 - 1e-2,  length.out = 1000)       # xL < xk
xk_seq <- seq(xL0 + 2e-2, xl0 - 2e-2,  length.out = 1000)       # xL < xk < xl
xl_seq <- seq(xk0 + 2e-2, xU0 - 2e-2,  length.out = 1000)       # xk < xl < xU
xU_seq <- seq(xl0 + 2e-2, xl0 + 20,    length.out = 1000)       # xU > xl
seqs   <- list(xL=xL_seq, xk=xk_seq, xl=xl_seq, xU=xU_seq)

# Safe residual evaluator (returns c(r1,r2,r3,r4) or NA on failure)
residuals_safe <- function(x, p, tol_build=1e-10) {
  tryCatch(
    opt_conditions_residuals(
      z_from_x(x$xL, x$xk, x$xl, x$xU), p, tol_build = tol_build
    ),
    error = function(e) rep(NA_real_, 4)
  )
}

# Scan all four residuals while varying one variable
scan_all_residuals_vs <- function(var, grid, base_x, p) {
  R <- matrix(NA_real_, nrow = length(grid), ncol = 4)
  x <- base_x
  for (i in seq_along(grid)) {
    x[[var]] <- grid[i]
    R[i, ] <- residuals_safe(x, p)
  }
  list(x = grid, R = R)
}

# Run scans for each variable
scans <- lapply(names(seqs), function(v) scan_all_residuals_vs(v, seqs[[v]], base_x, p))
names(scans) <- names(seqs)

# 4×4 plots: each residual (rows) vs each variable (cols)
res_names <- paste0("Residual ", 1:4)
par(mfrow = c(4,4), mar = c(3,4,1.2,0.5) + 0.4)
for (ri in 1:4) {
  for (vj in names(scans)) {
    x <- scans[[vj]]$x
    y <- scans[[vj]]$R[, ri]
    plot(x, y, type = "l", xlab = vj, ylab = res_names[ri],
         main = paste(res_names[ri], "vs", vj))
    abline(h = 0, lty = 3, col = "#BBBBBB")
  }
}
