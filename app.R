# app.R — Shiny UI for ergodic singular control with ambiguous jump diffusion
# Uses your functions_solver.R and functions_images.R

suppressPackageStartupMessages({
  library(shiny)
})

# ---- source your scripts ----
# The app expects these two files to live in the same folder as app.R
source("functions_solver.R", local = TRUE)
source("functions_images.R", local = TRUE)

# ---- helpers ----
make_param_list <- function(input) {
  make_params(
    b     = input$b,
    delta = input$delta,
    r     = input$r,
    eps   = input$eps,
    sigma = input$sigma,
    mu    = input$mu,
    u     = input$u,
    l     = input$l
  )
}

solve_model_safe <- function(params, xL, gapU, use_nleqslv) {
  tryCatch({
    solve_two_stage(params, init = list(xL = xL, gapU = gapU), use_nleqslv = use_nleqslv)
  }, error = function(e) {
    structure(list(error = TRUE, message = conditionMessage(e)), class = "solve_error")
  })
}

# ---- UI ----
ui <- fluidPage(
  title = "Ergodic singular control — ambiguous jump diffusion",
  tags$head(tags$style(HTML(".small-note{color:#666;font-size:12px;} .ok{color:#2e7d32} .bad{color:#c62828}"))),
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("Model parameters"),
                 fluidRow(
                   column(6, numericInput("b",     HTML("b"),     value = -1.8, step = 0.05)),
                   column(6, numericInput("delta", HTML("&delta;"), value = 0.20, step = 0.01))
                 ),
                 fluidRow(
                   column(6, numericInput("r",     HTML("r"),     value = 1.00, step = 0.05)),
                   column(6, numericInput("eps",   HTML("&epsilon;"), value = 0.30, step = 0.01, min = 0, max = 0.999))
                 ),
                 fluidRow(
                   column(6, numericInput("sigma", HTML("&sigma;"), value = 0.60, step = 0.05, min = 1e-6)),
                   column(6, numericInput("mu",    HTML("&mu;"),    value = 1.50, step = 0.05, min = 1e-6))
                 ),
                 fluidRow(
                   column(6, numericInput("u", value = 1.00, label = HTML("u (lower cost)"), step = 0.1, min = 0)),
                   column(6, numericInput("l", value = 1.00, label = HTML("l (upper cost)"), step = 0.1, min = 0))
                 ),
                 tags$hr(),
                 h5("Solver options"),
                 fluidRow(
                   column(6, numericInput("xL0",  "Initial xL", value = -2.0, step = 0.1)),
                   column(6, numericInput("gapU", "Initial gapU", value = 1.0, step = 0.05, min = 1e-6))
                 ),
                 checkboxInput("use_nleq", "Use nleqslv (Broyden) if available", value = TRUE),
                 actionButton("solve", "Solve / Recompute", class = "btn-primary"),
                 tags$div(class = "small-note", "Tip: if the solver struggles, try nudging xL or gapU, or uncheck 'Use nleqslv'.")
    ),
    mainPanel(width = 8,
              tabsetPanel(id = "tabs", type = "tabs",
                          tabPanel("Solution & plots",
                                   br(),
                                   fluidRow(
                                     column(6,
                                            h4("Thresholds & status"),
                                            tableOutput("thresholds_tbl"),
                                            verbatimTextOutput("status_txt")
                                     ),
                                     column(6,
                                            h4("Gamma"),
                                            verbatimTextOutput("gamma_txt"),
                                            h5("Outer residuals at optimum (r3, r4)"),
                                            verbatimTextOutput("resid_txt")
                                     )
                                   ),
                                   tags$hr(),
                                   h4("Fig 1 — H(x) with thresholds"),
                                   plotOutput("fig1", height = 420)
                          ),
                          tabPanel("Simulation",
                                   br(),
                                   fluidRow(
                                     column(3, numericInput("T",  "Horizon T", value = 8, step = 0.5, min = 0.5)),
                                     column(3, numericInput("dt", "Step dt",  value = 0.001, step = 0.001, min = 1e-4)),
                                     column(3, numericInput("seed", "Seed", value = 123, step = 1)),
                                     column(3, actionButton("resim", "Run simulation", class = "btn-secondary"))
                                   ),
                                   tags$hr(),
                                   h4("Fig 2 — Reflected path of X"),
                                   plotOutput("fig2", height = 360),
                                   h4("Fig 3 — Cumulative controls"),
                                   plotOutput("fig3", height = 320)
                          ),
                          tabPanel("Diagnostics",
                                   br(),
                                   fluidRow(
                                     column(4, numericInput("tol", "Pass/fail tolerance (abs)", value = 1e-6, min = 0, step = 1e-6)),
                                     column(8, tags$div(class = "small-note", "Green = within tolerance; red = outside tolerance."))
                                   ),
                                   tags$hr(),
                                   h4("Printed diagnostic (from diagnose_solution)"),
                                   verbatimTextOutput("diag_print"),
                                   h4("11 checks"),
                                   tableOutput("diag_table")
                          )
              )
    )
  )
)

# ---- server ----
server <- function(input, output, session) {
  # compute solution on demand (and once at startup)
  sol_rv <- reactiveVal(NULL)
  
  compute_and_store <- function() {
    # Build params directly from current inputs (no reactive needed here)
    p <- make_params(
      b = input$b, delta = input$delta, r = input$r, eps = input$eps,
      sigma = input$sigma, mu = input$mu, u = input$u, l = input$l
    )
    withProgress(message = "Solving two-stage system...", value = 0.3, {
      sol <- solve_model_safe(p, xL = input$xL0, gapU = input$gapU, use_nleqslv = isTRUE(input$use_nleq))
      incProgress(0.7)
      if (!inherits(sol, "solve_error")) sol$params <- p
      sol_rv(sol)
    })
  }
  
  # Initial solve (run once when app starts)
  observeEvent(TRUE, { compute_and_store() }, once = TRUE)
  
  # Recompute on button click
  observeEvent(input$solve, { compute_and_store() })
  
  # convenience getter with guards
  sol_ok <- reactive({
    s <- sol_rv()
    req(!inherits(s, "solve_error"))
    req(!is.null(s$thresholds))
    xs <- unlist(s$thresholds, use.names = FALSE)
    req(all(is.finite(xs)))
    s
  })
  
  # status & scalars
  output$status_txt <- renderPrint({
    s <- sol_rv()
    if (inherits(s, "solve_error")) {
      cat("❌ ", s$message, sep = "")
    } else if (is.null(s)) {
      cat("(no solution yet)")
    } else {
      cat("Converged:", s$status$converged, "-", s$status$message, "
")
    }
  })
  
  output$gamma_txt <- renderPrint({
    s <- sol_ok()
    cat(sprintf("gamma = %.6f
", s$gamma))
  })
  
  output$resid_txt <- renderPrint({
    s <- sol_ok()
    r <- s$diagnostics$outer_residuals
    cat(sprintf("r3 = %.3e, r4 = %.3e
", r[1], r[2]))
  })
  
  output$thresholds_tbl <- renderTable({
    s <- sol_ok()
    data.frame(Threshold = names(s$thresholds), Value = unlist(s$thresholds), row.names = NULL)
  }, digits = 6)
  
  # --- plots ---
  output$fig1 <- renderPlot({
    s <- sol_ok(); p <- s$params
    req(all(is.finite(unlist(s$thresholds, use.names = FALSE))))
    draw_fig1(sol = s, params = p, show = TRUE, save = FALSE, dir = tempdir())
  })
  
  # simulation (depends on solution)
  sim_react <- eventReactive(list(input$resim, sol_ok(), input$T, input$dt, input$seed), {
    s <- sol_ok(); p <- s$params
    set.seed(input$seed)
    simulate_reflected_jd(T = input$T, dt = input$dt, x0 = NULL, params = p, thresholds = s$thresholds)
  }, ignoreInit = FALSE)
  
  output$fig2 <- renderPlot({
    s <- sol_ok(); p <- s$params; sim <- sim_react()
    draw_fig2(sol = s, params = p, sim = sim, show = TRUE, save = FALSE, dir = tempdir())
  })
  
  output$fig3 <- renderPlot({
    sim <- sim_react()
    draw_fig3(sim = sim, show = TRUE, save = FALSE, dir = tempdir())
  })
  
  # diagnostics (printed text + table with pass/fail)
  diag <- reactive({
    s <- sol_ok(); p <- s$params
    txt <- capture.output(d <- diagnose_solution(s, p))
    list(text = paste(txt, collapse = "
"), list = d)
  })
  
  output$diag_print <- renderPrint({
    cat(diag()$text)
  })
  
  output$diag_table <- renderTable({
    d <- diag()$list
    keep <- c(
      "H_xk",
      "IH_scaled_at_xlambda",
      "IH_unscaled_at_xlambda",
      "H_xL_plus_u",
      "Hp_xL",
      "val_jump_at_xk",
      "slope_jump_at_xk",
      "val_jump_at_xlambda",
      "slope_jump_at_xlambda",
      "H_xU_minus_l",
      "Hp_xU"
    )
    vals <- unlist(d)[keep]
    df <- data.frame(
      check = keep,
      value = as.numeric(vals),
      pass  = ifelse(abs(as.numeric(vals)) <= input$tol, "✅", "❌"),
      stringsAsFactors = FALSE
    )
    df
  }, digits = 6)
}

# ---- run ----
shinyApp(ui, server)
