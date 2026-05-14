library(shiny)
library(bslib)
library(DT)

results <- read.csv("results_tibble.csv", stringsAsFactors = FALSE)

results$type_I  <- round(results$type_I,  3)
results$type_II <- round(results$type_II, 3)
results$type_I  <- ifelse(is.na(results$type_I),  "—", results$type_I)
results$type_II <- ifelse(is.na(results$type_II), "—", results$type_II)

method_labels <- c(
  "No Covariates (NoC)"                              = "no_covs",
  "All Covariates (AllC)"                            = "all_covs",
  "P-Hacking (Hack)"                                 = "p_hacked",
  "Single Covariate Linear Model without X (LM 1C)"  = "r",
  "Single Covariate Linear Model with X (LM 1C + X)" = "partial_r",
  "All Covariates Linear Model without X (LM AllC)"  = "full_lm_wo_x",
  "All Covariates Linear Model with X (LM AllC + X)" = "full_lm",
  "All Covariates LASSO without X (LASSO AllC)"      = "lasso_wo_x",
  "All Covariates LASSO with X (LASSO AllC + X)"     = "lasso"
)

method_abbrevs <- c(
  "no_covs"      = "NoC",
  "all_covs"     = "AllC",
  "p_hacked"     = "Hack",
  "r"            = "LM 1C",
  "partial_r"    = "LM 1C + X",
  "full_lm_wo_x" = "LM AllC",
  "full_lm"      = "LM AllC + X",
  "lasso_wo_x"   = "LASSO AllC",
  "lasso"        = "LASSO AllC + X"
)

ui <- page_navbar(
  title = tags$span("Type I and Type II Error Rates by Covariate Selection Method and Research Context", style = "font-size: 1.4rem; font-weight: 600;"),
  theme = bs_theme(
    bg        = "#ffffff",
    fg        = "#343a40",
    primary   = "#C5050C",
    base_font = font_google("Source Sans Pro"),
    "navbar-bg" = "#C5050C"
  ),
  nav_panel(NULL,
  layout_sidebar(
  sidebar = sidebar(width = 400,
    accordion(
      open = FALSE,
      accordion_panel("Method",
        checkboxGroupInput("method", label = NULL,
                choices  = method_labels,
                selected = method_labels)
      ),
      accordion_panel("Sample Size (n)",
        checkboxGroupInput("n_obs", label = NULL,
                choices  = sort(unique(results$n_obs)),
                selected = sort(unique(results$n_obs)))
      ),
      accordion_panel("No. of Available Covariates (nC)",
        checkboxGroupInput("n_covs", label = NULL,
                choices  = sort(unique(results$n_covs)),
                selected = sort(unique(results$n_covs)))
      ),
      accordion_panel("Covariate-Outcome Strength (rCY)",
        checkboxGroupInput("r_ycov", label = NULL,
                choices  = sort(unique(results$r_ycov)),
                selected = sort(unique(results$r_ycov)))
      ),
      accordion_panel("% Covariates with Non-Zero Effects (%C)",
        checkboxGroupInput("p_good_covs", label = NULL,
                choices  = sort(unique(results$p_good_covs)),
                selected = sort(unique(results$p_good_covs)))
      ),
      accordion_panel("True IV Effect (BX)",
        checkboxGroupInput("b_x", label = NULL,
                choices  = sort(unique(results$b_x)),
                selected = sort(unique(results$b_x)))
      )
    ),
    hr(),
    actionButton("reset", "Reset All Filters", width = "100%")
  ),
  card(
    full_screen = TRUE,
    card_header("Error rates estimated across 40,000 simulations per condition"),
    DTOutput("table")
  )
  ))
)

server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateCheckboxGroupInput(session, "method",      selected = method_labels)
    updateCheckboxGroupInput(session, "n_obs",       selected = sort(unique(results$n_obs)))
    updateCheckboxGroupInput(session, "n_covs",      selected = sort(unique(results$n_covs)))
    updateCheckboxGroupInput(session, "r_ycov",      selected = sort(unique(results$r_ycov)))
    updateCheckboxGroupInput(session, "p_good_covs", selected = sort(unique(results$p_good_covs)))
    updateCheckboxGroupInput(session, "b_x",         selected = sort(unique(results$b_x)))
  })

  filtered <- reactive({
    df <- results
    df <- df[df$method      %in% input$method,                  ]
    df <- df[df$n_obs       %in% as.numeric(input$n_obs),       ]
    df <- df[df$n_covs      %in% as.numeric(input$n_covs),      ]
    df <- df[df$r_ycov      %in% as.numeric(input$r_ycov),      ]
    df <- df[df$p_good_covs %in% as.numeric(input$p_good_covs), ]
    df <- df[df$b_x         %in% as.numeric(input$b_x),         ]
    df$method <- method_abbrevs[df$method]
    df
  })

  output$table <- renderDT({
    df <- filtered()
    names(df) <- c("Method", "n", "nC", "rCY", "%C", "BX", "Type I", "Type II")
    datatable(
      df,
      rownames = FALSE,
      filter   = "top",
      options  = list(
        pageLength      = 15,
        searchHighlight = TRUE,
        dom             = "tip",
        scrollX         = TRUE,
        autoWidth  = TRUE,
        columnDefs = list(
          list(searchable = FALSE, targets = c(6, 7)),
          list(width = "110px", targets = 0),
          list(width = "45px",  targets = c(1, 2, 3, 4, 5)),
          list(width = "65px",  targets = c(6, 7))
        )
      )
    )
  })
}

shinyApp(ui, server)
