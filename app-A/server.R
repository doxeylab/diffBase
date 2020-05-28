#-------------------------------------------------------------------------------
# server.R
# Last modified: 2020-05-28 23:15:37 (CEST)
# BJM Tremblay

msg("Loading server")
server <- function(input, output, session) {

  msg("Session start:", session$token)
  onSessionEnded({function() msg("Session stop:", session$token)})

  CURRENT_PAGE <- reactiveValues(
    WHICH = "WELCOME"
  )

  output$CURRENT_PAGE <- eventReactive(CURRENT_PAGE$WHICH, {
    CURRENT_PAGE$WHICH
  })
  outputOptions(output, "CURRENT_PAGE", suspendWhenHidden = FALSE)

  SELECTED_TYPE <- reactiveValues(
    WHICH = "A1"
  )

  SELECTED_SUBTYPE <- reactiveValues(
    "A1" = "A1.1",
    "A2" = "A2.1",
    "A3" = "A3.1",
    "A4" = "A4.1",
    "A5" = "A5.1",
    "A6" = "A6.1",
    "A7" = "A7.1",
    "sordellii_TcsH" = "sordellii_TcsH.1"
  )

  QUERY <- reactiveValues(STRING = NULL)
  observe({
    QUERY$STRING <- parseQueryString(session$clientData$url_search)[["goto"]]
  })

  observeEvent(QUERY$STRING, {
    req(QUERY$STRING)
    if (!is.null(QUERY$STRING)) {
      SType <- QUERY$STRING
      print(SType)
      QType <- strsplit(SType, ".", fixed = TRUE)[[1]][1]
      if (!QType %in% ALL_TYPES || !SType %in% ALL_TYPES_SUBTYPES[[QType]]) {
        showModal(modalDialog(
          title = "Error", paste("Incorrect query.", SType, "does not exist.")
        ))
      } else {
        SELECTED_TYPE$WHICH <- QType
        SELECTED_SUBTYPE[[QType]] <- SType
        CURRENT_PAGE$WHICH <- "INFO"
        updateSelectInput(
          session, "SUBTYPE_SELECTOR",
          label = "",
          choices = names(SEQ_NAMES_LIST[[SELECTED_TYPE$WHICH]]),
          selected = SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
        )
      }
    }
  })

  output$PANEL_LEFT_CURRENT_TYPE <- renderText({
    paste("<b>Type:</b>", SELECTED_TYPE$WHICH)
  })

  output$PANEL_LEFT_CURRENT_SUBTYPE <- renderText({
    paste(
      "<b>Currently selected subtype:</b>",
      SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
    )
  })

  output$PANEL_TOP_RIGHT_CURRENT_SUBTYPE <- renderText({
    if (CURRENT_PAGE$WHICH != "INFO") return()
    paste(
      "<b>Currently selected subtype sequence:</b>",
      SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
    )
  })

  output$PANEL_TOP_RIGHT_METADATA <- DT::renderDataTable({
    if (CURRENT_PAGE$WHICH != "INFO") return()
    out <- show_metadata(
      SEQ_NAMES_ALL[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]]
    )
    if (is.null(out)) return()
    DT::datatable(
      out,
      selection = "none",
      options = list(
        pageLength = 10,
        ordering = FALSE,
        scrollX = TRUE,
        dom = "tip"
      ),
      rownames = FALSE
    )
  })

  observe({
    lapply(
      ALL_TYPES,
      function(x) {
        observeEvent(input[[paste0("BUTTON_", x)]], {
          CURRENT_PAGE$WHICH <- "INFO"
          SELECTED_TYPE$WHICH <- x
        })
      }
    )
  })

  observeEvent(input$BUTTON_GO_BACK_TO_WELCOME, {
    CURRENT_PAGE$WHICH <- "WELCOME"
  })

  observeEvent(input$BUTTON_GO_BACK_TO_WELCOME_BLAST, {
    CURRENT_PAGE$WHICH <- "WELCOME"
  })

  observeEvent(input$SUBTYPE_SELECTOR, {
    SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]] <- input$SUBTYPE_SELECTOR
  })

  output$PANEL_TOP_RIGHT_CURRENT_SUBTYPE_SEQUENCE <- renderText({
    if (CURRENT_PAGE$WHICH != "INFO") return()
    as.character(SEQS_ALL[[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]]])
  })

  output$CURRENT_ACCESSION <- renderText({
    paste(
      "<b>Representative sequence:</b>",
      names(SEQ_NAMES_ALL[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]])
    )
  })

  observeEvent(input$BLASTP_BUTTON, {
    req(input$BLASTP_INPUT)
    if (!CONFIGS$UseBlastp) {
      showModal(modalDialog(title = "Error",
        "Blastp has been disabled. Edit the config file to get it working again."
      ))
      return()
    }
    res <- run_blast(input$BLASTP_INPUT)
    if (is.null(res)) {
      showModal(modalDialog(title = "BLASTP",
        "No hits were detected."
      ))
      return()
    } else if (is.data.frame(res[[1]])) {
      output$BLASTP_DOWNLOAD <- downloadHandler(
        filename = "blastp_results.tsv",
        content = function(con) readr::write_tsv(res[[1]][, -ncol(res[[1]])], con)
      )
      output$BLASTP_DOWNLOAD_RAW <- downloadHandler(
        filename = "blastp_results_alignment.txt",
        content = function(con) readr::write_lines(res[[2]], con)
      )
      output$BLASTP_RES_TABLE <- DT::renderDataTable({
        DT::datatable(
          res[[1]],
          extensions = "Buttons",
          escape = FALSE,
          selection = "none",
          options = list(
            lengthMenu = list(c(10, 50, 100, 200), c("10", "50", "100", "200")),
            pageLength = 10,
            dom = "ltip",
            columnDefs = list(
              list(targets = ncol(res), bSortable = FALSE, className = "dt-center")
            )
          )
        )
      })
      output$BLASTP_SHOW_QUERY <- renderText(input$BLASTP_INPUT)
      CURRENT_PAGE$WHICH <- "BLASTP_RES"
    }
  })

  observe({
    updateSelectInput(
      session, "SUBTYPE_SELECTOR",
      label = "",
      choices = names(SEQ_NAMES_LIST[[SELECTED_TYPE$WHICH]]),
      selected = SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
    )
  })

  TREE_STATE <- reactiveValues(SHOW = FALSE)

  observeEvent(input$TREE_BUTTON, {
    if (TREE_STATE$SHOW)
      TREE_STATE$SHOW <- FALSE
    else
      TREE_STATE$SHOW <- TRUE
  })

  observeEvent(CURRENT_PAGE$WHICH, {
    TREE_STATE$SHOW <- FALSE
  })
  output$SHOW_PLOT <- eventReactive(TREE_STATE$SHOW, {
    TREE_STATE$SHOW
  })
  outputOptions(output, "SHOW_PLOT", suspendWhenHidden = FALSE)

  output$PANEL_PLOT <- renderPlot({
    if (TREE_STATE$SHOW)
      switch(CURRENT_PAGE$WHICH,
        WELCOME = plot_tree(),
        INFO = plot_tree(SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]),
        BLASTP_RES = plot_tree()
      )
    else
      NULL
  })

  observeEvent(input$COMMUNITY_BUTTON, {
    req(input$COMMUNITY_INPUT)
    email <- clean_email(input$COMMUNITY_EMAIL)
    comm_text <- input$COMMUNITY_INPUT
    modal_text <- "Your message has been sent."
    if (nchar(email) < 5 || !check_email(email)) {
      modal_text <- paste(
        modal_text,
        "[WARNING: invalid email. Try again with a valid email if you would like",
        "to receive a response.]"
      )
      email <- "anonymous"
    }
    msg("Saving community message locally")
    cat(email, comm_text, sep = "\n",
        file = paste0("community/", as.integer(Sys.time()), ".txt"))
    if (CONFIGS$UseGmail) {
      msg("Sending community message by email")
      gmailr::gm_send_message(
        gmailr::gm_mime() %>%
          gmailr::gm_from(CONFIGS$ServerEmail) %>%
          gmailr::gm_to(CONFIGS$ServerEmail) %>%
          gmailr::gm_subject(paste("Community message from", email)) %>%
          gmailr::gm_text_body(comm_text)
      )
    }
    showModal(modalDialog(title = "Thank you", modal_text))
  })

  output$DOWNLOAD_METADATA <- downloadHandler(
    filename = paste0(
      "ToxinB_", SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]],
      "_Associated_Metadata.tsv"
    ),
    content = function(con) readr::write_tsv(show_metadata(
      SEQ_NAMES_ALL[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]]
    ), con)
  )

  output$DOWNLOAD_ALL <- downloadHandler(
    filename = "ToxinA_sequences.fa",
    content = function(con) {
      readr::write_lines(
        readr::read_lines("downloads/ALL-sequences.fa"), con
      )
    }
  )

  output$DOWNLOAD_TYPE <- downloadHandler(
    filename = function() paste0(
      "ToxinASubtype", SELECTED_TYPE$WHICH, "_sequences.fa"
    ),
    content = function(con) {
      readr::write_lines(
        readr::read_lines(
          paste0("downloads/", SELECTED_TYPE$WHICH, "-sequences.fa")
        ), con
      )
    }
  )

  observeEvent(input$BLASTP_GOTO, {
    CURRENT_PAGE$WHICH <- "INFO"
    SELECTED_TYPE$WHICH <- strsplit(input$BLASTP_GOTO, "_", fixed = TRUE)[[1]][3]
    SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]] <- gsub("BLASTP.GOTO.", "",
      gsub("_", ".", input$BLASTP_GOTO, fixed = TRUE),
      fixed = TRUE
    )
  })

}

msg("Launching app")
server
