# app.R
# Minimal Shiny dashboard for massiveNGSpipe runtime monitoring (full-width controls row)

# ---- load massiveNGSpipe (dev path first if present) ----

#'
#' @import shiny
mNGSp_app <- function(config = pipeline_config(),
                      pipelines_all = pipeline_init_all(config,
      only_complete_genomes = TRUE, gene_symbols = FALSE, show_status_per_exp = FALSE)) {
  pipeline_sessions <- session_info_table(config)
  default_show_done  <- TRUE
  default_show_stats <- FALSE

  # ---- ui ----
  ui <- fluidPage(
    titlePanel("massiveNGSpipe monitor"),

    # --- FULL-WIDTH SETTINGS ROW (replaces sidebar) ---
    fluidRow(
      column(
        width = 12,
        wellPanel(
          fluidRow(
            column(
              width = 3,
              # Project (single option)
              selectizeInput(
                "project", "Project",
                choices = c("NGS_pipeline"),
                selected = "NGS_pipeline",
                options = list(create = FALSE)
              )
            ),
            column(
              width = 3,
              # Session selector
              selectInput(
                "session_selected",
                "Select pipeline session",
                choices = setNames(seq_len(nrow(pipeline_sessions)), basename(pipeline_sessions$path)),
                selected = 1
              )
            ),
            column(
              width = 3,
              # Toggles
              checkboxInput("show_done",  "Include done steps", value = default_show_done),
              checkboxInput("show_stats", "Show stats",          value = default_show_stats)
            ),
            column(
              width = 3,
              # Refresh controls
              sliderInput("refresh_ms", "Auto-refresh interval (s)",
                          min = 1, max = 30, value = 5, step = 1),

            )
          ),
          # --- INFO ROW: overall/system + this-session usage ---
          fluidRow(
            column(
              width = 8,
              h4("System usage (overall)"),
              verbatimTextOutput("system_usage")
            ),
            column(
              width = 4,
              h4("This session usage"),
              tableOutput("session_usage")
            )
          )
        )
      )
    ),

    # --- FIRST CONTENT: PLOTLY STATUS PLOT ---
    fluidRow(
      column(
        width = 12,
        plotlyOutput("status_plot", height = "600px")
      )
    ),

    # --- PROCESSES TABLE ---
    fluidRow(
      column(
        width = 12,
        h4("Running processes for selected session"),
        tableOutput("proc_table")
      )
    )
  )

  # ---- server ----
  server <- function(input, output, session) {

    # helper: manual or timer invalidation
    invalidate_tick <- reactive({
      invalidateLater(input$refresh_ms * 1e3, session)
      Sys.time()
    })

    # steps (static from config)
    output$pipeline_steps <- renderText({
      paste(names(config$flag), collapse = ", ")
    })

    # overall system usage (one-liner string)
    output$system_usage <- renderText({
      invalidate_tick()
      tryCatch(
        capture.output(get_system_usage(one_liner = TRUE))[1],
        error = function(e) paste("System usage unavailable:", conditionMessage(e))
      )
    })

    # selected session info -> pipelines, status, processes
    sel <- reactive({
      req(nrow(pipeline_sessions) > 0)
      idx <- as.integer(input$session_selected)
      validate(need(idx >= 1 && idx <= nrow(pipeline_sessions), "Invalid session selection"))
      list(
        idx  = idx,
        path = pipeline_sessions$path[idx],
        name = pipeline_sessions$name[idx]
      )
    })

    # read session + subset pipelines
    session_obj <- reactive({
      invalidate_tick()
      s <- session_info_read(list(session_dir = sel()$path))
      keep <- intersect(s$pipeline_names, names(pipelines_all))
      pl   <- pipelines_all[keep]
      list(s = s, pipelines = pl)
    })

    # compute status table and plot
    status_reactive <- reactive({
      invalidate_tick()
      s <- session_obj()
      req(length(s$pipelines) > 0)
      status <- status_per_study(
        s$pipelines, config,
        show_status_per_exp = FALSE,
        show_done           = isTRUE(input$show_done),
        show_stats          = isTRUE(input$show_stats)
      )
      status
    })

    output$status_plot <- renderPlotly({
      st <- status_reactive()
      validate(need(NROW(st) > 0, "No status to plot for this session"))
      status_plot(st)  # returns plotly object
    })

    # per-session processes + summed CPU/MEM
    session_ps <- reactive({
      invalidate_tick()
      s <- session_obj()$s
      tryCatch(
        get_running_processes(s$pgid),
        error = function(e) data.table(message = conditionMessage(e))
      )
    })

    output$proc_table <- renderTable({
      session_ps()
    })

    output$session_usage <- renderTable({
      ps <- session_ps()
      if (!all(c("%CPU", "%MEM") %in% names(ps))) {
        return(data.frame(Metric = c("CPU", "MEM"), Value = c("NA", "NA")))
      }
      data.frame(
        Metric = c("CPU % (sum)", "MEM % (sum)"),
        Value  = c(sum(ps$`%CPU`, na.rm = TRUE), sum(ps$`%MEM`, na.rm = TRUE)),
        check.names = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s")
  }

  shinyApp(ui, server)
}
