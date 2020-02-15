#-------------------------------------------------------------------------------
# server.R
# Last modified: 2020-02-15 10:29:48 (CET)
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
    WHICH = "A"
  )

  SELECTED_SUBTYPE <- reactiveValues(
    A = "A.1",
    B = "B.1",
    C = "C.1",
    D = "D.1",
    E = "E.1",
    F = "F.1",
    G = "G.1",
    H = "H.1",
    I = "I.1"
  )

  output$PANEL_LEFT_CURRENT_TYPE <- renderText({
    paste("<b>Toxinotype:</b>", SELECTED_TYPE$WHICH)
  })

  output$PANEL_LEFT_CURRENT_SUBTYPE <- renderText({
    paste(
      "<b>Currently selected subtype:</b>",
      SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
    )
  })

  output$PANEL_TOP_RIGHT_CURRENT_SUBTYPE <- renderText({
    paste(
      "<b>Currently selected subtype sequence:</b>",
      SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]
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
    SEQS_ALL[[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]]]
  })

  output$CURRENT_ACCESSION <- renderText({
    SEQ_NAMES_ALL[SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]]
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
    if (is.data.frame(res)) {
      output$BLASTP_DOWNLOAD <- downloadHandler(
        filename = "blastp_results.tsv",
        content = function(con) readr::write_tsv(res[, -ncol(res)], con)
      )
      output$BLASTP_RES_TABLE <- DT::renderDataTable({
        DT::datatable(
          res,
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

  output$PANEL_PLOT <- renderPlot({
    switch(CURRENT_PAGE$WHICH,
      WELCOME = plot_tree(),
      INFO = plot_tree(SELECTED_SUBTYPE[[SELECTED_TYPE$WHICH]]),
      BLASTP_RES = plot_tree()
    )
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
