#' For each study that is done, validate / fix shifts
#' @import shinycssloaders shinyjqui shinyWidgets base64enc
#' @importFrom DT datatable renderDT DTOutput dataTableProxy
#' @export
psite_reshift_app <- function(config, exp_names = reshift_app_exp_names(config),
                              options = list("launch.browser" = TRUE)) {
  stopifnot(is.null(config) || is.list(config))
  stopifnot(is.null(options) || is.list(options))
  ui <- fluidPage(
    fluidRow(h3("Ribo-seq P-site Shift Table Editor")),
    fluidRow(
      column(1, actionButton("save", "💾 Save", class = "btn-primary")),
      column(1, actionButton("verify", "✔️ Verify only", class = "btn-secondary")),
      column(1, actionButton("undo", "↶ Undo shift", class = "btn-default")), # NEW BUTTON
      column(3, fluidRow(selectizeInput(
        inputId = "study",
        label = NULL,
        choices = exp_names,
        multiple = FALSE,
        width = "100%",
        options = list(placeholder = "Select a study")
      )), uiOutput("shift_verification")),
      column(2, selectizeInput("libs", "Libs (subset)", choices = NULL, multiple = TRUE, width = "100%")),
      column(2, shinyWidgets::pickerInput(
        inputId = "fractions_to_shift",
        label = "Select fractions to shift",
        choices = c(20,21, 25:33),
        multiple = TRUE,
        selected = 26,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          `selected-text-format` = "count > 15"
        ),
        width = "100%"
      )),
      column(1, numericInput("relative_shift", "Shift (Rel)", value = 3, step = 1, width = "100%")),
      column(1, fluidRow(actionButton("apply_shift", "↻ Apply Shift", class = "btn-warning")))
    ),
    fluidRow(
      column(6,
             div(
               style = "height: 800px; overflow-y: auto; border: 1px solid #ddd; padding: 5px;",
               uiOutput("image_ui")
             ) %>% shinycssloaders::withSpinner(color = "#0dc5c1")),
      column(6, DTOutput("table_ui"))
    )
  )

  server <- function(input, output, session) {
    rv <- reactiveValues(
      study = NULL,
      lib_paths = NULL,
      lib_names = NULL,
      current_table_path = NULL,
      flat_table = NULL,
      history = list() # NEW: stores past table states
    )
    verified_status <- reactiveVal(FALSE)

    check_verified_file <- function() {
      req(input$study)
      manual_shifts <- file.path(dirname(input$study), "manually_checked_shifts.rds")
      verified_status(file.exists(manual_shifts))
    }

    observeEvent(input$study, {
      req(input$study)
      message("Selected study: ", input$study)
      check_verified_file()
      rv$study <- basename(dirname(dirname(dirname(input$study))))
      shift_table_path <- file.path(dirname(dirname(input$study)),
                                    "pshifted", "shifting_table.rds")
      rv$current_table_path <- shift_table_path

      if (!file.exists(shift_table_path)) {
        showNotification(paste("Missing shift table:", shift_table_path), type = "error")
        rv$flat_table <- NULL
        return()
      }
      shift_list <- readRDS(shift_table_path)
      shift_list_ids <- sub("\\.ofst$", "", basename(names(shift_list)))
      df <- read.experiment(rv$study, validate = FALSE)
      run_ids <- runIDs(df)
      shift_list_corrected_order <- match(run_ids, shift_list_ids)
      stopifnot(all(!is.na(shift_list_corrected_order)))
      shift_list <- shift_list[shift_list_corrected_order]
      shift_table <- rbindlist(shift_list, idcol = "lib")
      rv$lib_paths <- unique(shift_table$lib)
      rv$lib_names <- bamVarName(df)
      updateSelectizeInput(inputId = "libs", choices = rv$lib_names, selected = NULL)
      shift_table <- cbind(libname = rv$lib_names[as.numeric(factor(shift_table$lib, levels = unique(shift_table$lib)))],
                           shift_table)
      shift_table[, offset_original := offsets_start]
      shift_table[, lib := NULL]

      rv$flat_table <- shift_table
      rv$history <- list(copy(shift_table)) # store initial state
    }, ignoreInit = FALSE)

    output$image_ui <- renderUI({
      req(input$study)
      img_path <- input$study
      format <- ifelse(tools::file_ext(img_path) == "fst", "plotly", "image")
      if (format == "plotly") {
        message("Using plotly, not implemented yet!")
        img_path <- sub("\\.fst$", "", img_path)
        if (!file.exists(img_path)) return(tags$p("Image not found"))
        tags$img(src = base64enc::dataURI(file = img_path, mime = "image/png"), width = "100%")
      } else {
        if (!file.exists(img_path)) return(tags$p("Image not found"))
        tags$img(src = base64enc::dataURI(file = img_path, mime = "image/png"), width = "100%")
      }
    })

    output$table_ui <- renderDT({
      req(rv$flat_table)
      datatable(rv$flat_table, editable = TRUE, options = list(pageLength = 15))
    })

    proxy <- dataTableProxy("table_ui")

    observeEvent(input$table_ui_cell_edit, {
      info <- input$table_ui_cell_edit
      rv$flat_table[info$row, info$col] <- info$value
      rv$history <- c(rv$history, list(copy(rv$flat_table))) # save after edit
      replaceData(proxy, rv$flat_table, resetPaging = FALSE)
    })

    observeEvent(input$undo, {
      if (length(rv$history) > 1) {
        rv$history <- rv$history[-length(rv$history)] # drop latest
        rv$flat_table <- copy(rv$history[[length(rv$history)]])
        replaceData(proxy, rv$flat_table, resetPaging = FALSE)
        showNotification("Undid last change", type = "message")
      } else {
        showNotification("Nothing to undo", type = "warning")
      }
    })

    observeEvent(input$save, {
      req(rv$current_table_path, rv$flat_table)
      tab <- copy(rv$flat_table)
      if (any(tab$offsets_start > 0)) {
        shiny::showModal(modalDialog("You can not have shifts > 0, please fix!"))
        return()
      }
      tab[, lib := rv$lib_paths[as.numeric(factor(libname, levels = unique(libname)))]]
      tab[, libname := NULL]
      tab[, offset_original := NULL]
      split_list <- split(tab, by = "lib", keep.by = FALSE)
      split_list <- lapply(split_list, setDT)

      saveRDS(split_list, rv$current_table_path)
      saveRDS(TRUE, file.path(dirname(input$study), "manually_checked_shifts.rds"))
      check_verified_file()
      steps_to_remove <- names(config$flag)[-seq(8)]
      stopifnot(steps_to_remove[1] == "pshifted")
      remove_flag_all_exp(config, steps = steps_to_remove, exps = isolate(rv$study))
      showNotification(paste("Saved and reset flags for:", rv$study), type = "message")
    })

    observeEvent(input$verify, {
      req(rv$current_table_path, rv$flat_table)
      saveRDS(TRUE, file.path(dirname(input$study), "manually_checked_shifts.rds"))
      check_verified_file()
      showNotification(paste("Verified shifts for: ", rv$study), type = "message")
    })

    output$shift_verification <- renderUI({
      if (verified_status()) {
        div("✔ Verified", style = "color: white; background-color: #28a745;
         padding: 6px 12px; border-radius: 5px; font-weight: bold; text-align: center;")
      } else {
        div("✘ Not Verified", style = "color: white; background-color: #dc3545;
         padding: 6px 12px; border-radius: 5px; font-weight: bold; text-align: center;")
      }
    })

    observeEvent(input$apply_shift, {
      req(rv$flat_table, input$fractions_to_shift, input$relative_shift)
      rows_to_shift <- rv$flat_table$fraction %in% input$fractions_to_shift
      if (isTruthy(input$libs))
        rows_to_shift <- rows_to_shift & (rv$flat_table$libname %in% input$libs)
      n <- sum(rows_to_shift)

      if (n == 0) {
        showNotification(paste("No rows found with fraction =",
                               paste(input$fractions_to_shift, collapse = ", ")), type = "warning")
        return()
      }

      rv$history <- c(rv$history, list(copy(rv$flat_table))) # save before shift
      rv$flat_table[rows_to_shift, offsets_start := offsets_start + input$relative_shift]
      replaceData(proxy, rv$flat_table, resetPaging = FALSE)

      showNotification(paste("Shifted", n, "rows for fraction",
                             paste(input$fractions_to_shift, collapse = ", "),
                             "by", input$relative_shift), type = "message")
    })
  }

  shinyApp(ui, server, options = options)
}

#' @export
reshift_app_exp_names <- function(config, pipelines = pipeline_init_all(config, only_complete_genomes = TRUE, gene_symbols = FALSE,
                                                                        verbose_load_annotation = FALSE, show_status_per_exp = FALSE),
                                  species_subset = NULL, study_subset = NULL) {
  if (!is.null(study_subset)) {
    if (all(study_subset) %in% names(pipelines)) stop("When study subset is not NULL,
                                                      it must be a subset of the names of 'pipelines'")
    pipelines <- pipelines[subset]
  }
  vec <- progress_report(pipelines, config, return_progress_vector = TRUE,
                         named_progress_vector = TRUE, show_status_per_exp = FALSE)
  # Only use subset that is done reprocessing
  vec <- vec[vec == 14]
  shift_plot_paths <- file.path(ORFik::config()["bam"], names(vec), "aligned", "QC_STATS", "pshifts_barplots.png")
  shift_plot_paths_fst <- ((paste0(shift_plot_paths, ".fst")))
  shift_plot_paths[file.exists(shift_plot_paths_fst)] <- shift_plot_paths_fst[file.exists(shift_plot_paths_fst)]
  stopifnot(all(file.exists(shift_plot_paths)))
  # Ignore manually checked ones
  checked_shifts_rds <- file.path(dirname(shift_plot_paths), "manually_checked_shifts.rds")

  exp_names <- shift_plot_paths
  names(exp_names) <- names(vec) # Display names to show (i.e. shorter)
  # Remove ones that are manually checked already
  message("Studies with validated pshifts: ", sum(file.exists(checked_shifts_rds)), " / ", length(exp_names))
  exp_names <- exp_names[!file.exists(checked_shifts_rds)]
  exp_names <- exp_names[order(sub(".*-", "", names(exp_names)))]
  if (!is.null(species_subset)) {
    species_subset <- gsub(" ", "_", tolower(species_subset))
    species <- sub(".*-", "", names(exp_names))
    stopifnot(all(species_subset %in% species))
    exp_names <- exp_names[species %in% species_subset]
  }
  message("Number of final studies: ", length(exp_names))
  message("Studies per organism")
  print(table(sub(".*-", "", names(exp_names))))
  return(exp_names)
}
