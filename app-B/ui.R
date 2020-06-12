#-------------------------------------------------------------------------------
# ui.R
# Last modified: 2020-06-12 15:59:32 (CEST)
# BJM Tremblay

msg("Loading UI")
ui <- function(request) fluidPage(

  theme = shinytheme("simplex"),

  title = "DiffBase: Toxin B",

  br(), br(), br(),

  navbarPage(
    tags$li(
      class = "dropdown",
      tags$style(".navbar {min-height:75px}")
    ),
    title = div(
      img(src = "DiffBaseLogo.png", height = 53, width = 100),
      HTML(paste0(
        "<a href='", CONFIGS$URL, "'>DiffBase</a>"
      ))
    ),
    id = "NAVBARPAGE",
    position = "fixed-top",
    collapsible = FALSE,
    footer = tags$footer(
      tags$div(HTML(paste0(
        "<p>For questions or concerns, email <a href='mailto:",
        CONFIGS$ServerEmail, "'>", CONFIGS$ServerEmail, "</a>."
      ))),
      align = "center",
      style = "
        position: static;
        bottom: 0;
        width: 100%;
        height: 50px;
        color: grey;
        padding: 10px;
        z-index: 1000;
      "
    ),

    br(), br(),
    tags$div(HTML("<h2>Welcome to the <i>C. difficile</i> toxin B database</h2>")),
    br(), br(),

    fluidRow(

      column(3, # PANEL_LEFT
        wellPanel(

          conditionalPanel(
            condition = "output.CURRENT_PAGE == 'WELCOME'",
            tagList(
              h4("Toxin Subtype (TcdB)"),
              actionLink("BUTTON_B1", "B1:"), GET_MEMBER_COUNT("B1"), br(),
              actionLink("BUTTON_B2", "B2:"), GET_MEMBER_COUNT("B2"), br(),
              actionLink("BUTTON_B3", "B3:"), GET_MEMBER_COUNT("B3"), br(),
              actionLink("BUTTON_B4", "B4:"), GET_MEMBER_COUNT("B4"), br(),
              actionLink("BUTTON_B5", "B5:"), GET_MEMBER_COUNT("B5"), br(),
              actionLink("BUTTON_B6", "B6:"), GET_MEMBER_COUNT("B6"), br(),
              actionLink("BUTTON_B7", "B7:"), GET_MEMBER_COUNT("B7"), br(),
              actionLink("BUTTON_B8", "B8:"), GET_MEMBER_COUNT("B8"), br(),
              actionLink("BUTTON_B9", "B9:"), GET_MEMBER_COUNT("B9"), br(),
              actionLink("BUTTON_B10", "B10:"), GET_MEMBER_COUNT("B10"), br(),
              actionLink("BUTTON_B11", "B11:"), GET_MEMBER_COUNT("B11"), br(),
              actionLink("BUTTON_B12", "B12:"), GET_MEMBER_COUNT("B12"), br(),
              actionLink("BUTTON_sordellii_group", "sordellii_group:"), GET_MEMBER_COUNT("sordellii_group"), br(),
              actionLink("BUTTON_sordellii_TcsL", "sordellii_TcsL:"), GET_MEMBER_COUNT("sordellii_TcsL"),
              br(), br(),
              textInput("APP_SEARCH_BOX", "Query the database"),
              actionLink("SEARCH_BUTTON", "Search"),
              br(), br(),
              downloadLink(
                "DOWNLOAD_ALL", "Download all toxin B sequences"
              ), br(),
              downloadLink(
                "DOWNLOAD_REP", "Download representative sequences"
              )
            )
          ),

          conditionalPanel(
            condition = "output.CURRENT_PAGE == 'INFO'",
            make_type_info()
          ),

          conditionalPanel(
            condition = "output.CURRENT_PAGE == 'BLASTP_RES'",
            tagList(
              actionLink("BUTTON_GO_BACK_TO_WELCOME_BLAST", "Go back"), br(),
              br(),
              tags$b("Query:"), br(),
              br(),
              verbatimTextOutput("BLASTP_SHOW_QUERY")
            )
          ),

          conditionalPanel(
            condition = "output.CURRENT_PAGE == 'SEARCH_RES'",
            actionLink("BUTTON_GO_BACK_TO_WELCOME_SEARCH", "Go back")
          )

        ),
        wellPanel(
          tagList(
            h4("Database information"),
            HTML(paste("<b>Last updated data:</b>", LAST_UPDATE_DATE())),
            br(),
            HTML(paste("<b>Number of sequences:</b>", NUMBER_OF_SEQUENCES())),
            br(),
            HTML(paste("<b>Number of bacterial strains:</b>", NUMBER_OF_STRAINS()))
          )
        ),
        wellPanel(
          tagList(
            h4("Community involvement"),
            paste(
              "Do you have information to provide regarding new papers,",
              "sequences, or data concerning TcdA or B toxins? Feel free",
              "to share it with us:"
            ),
            textAreaInput("COMMUNITY_INPUT",
              "", value = "", placeholder = ""
            ),
            textInput("COMMUNITY_EMAIL",
              "Your email", ""
            ),
            actionLink("COMMUNITY_BUTTON", "Submit")
          )
        )
      ),

      column(9,

        column(12, # PANEL_TOP_RIGHT
          wellPanel(

            conditionalPanel(
              condition = "output.CURRENT_PAGE == 'WELCOME'",
              tagList(
                h4("BLASTP input"),
                textAreaInput("BLASTP_INPUT",
                  "Input a single amino acid sequence string:",
                  value = "", placeholder = ""
                ),
                actionLink("BLASTP_BUTTON", "Submit")
              )
            ),

            conditionalPanel(
              condition = "output.CURRENT_PAGE == 'INFO'",
              make_type_info_more()
            ),

            conditionalPanel(
              condition = "output.CURRENT_PAGE == 'BLASTP_RES'",
              downloadLink("BLASTP_DOWNLOAD", "Download all results"), br(),
              downloadLink("BLASTP_DOWNLOAD_RAW", "Download BLASTP alignment"), br(), br(),
              DT::dataTableOutput("BLASTP_RES_TABLE")
            ),

            conditionalPanel(
              condition = "output.CURRENT_PAGE == 'SEARCH_RES'",
              DT::dataTableOutput("SEARCH_RES_TABLE")
            )

          )
        ),

        column(12, # PANEL_BOTTOM_RIGHT
          wellPanel(
            tagList(
              h4("Phylogenetic tree"),
              actionLink("TREE_BUTTON", "Show/hide tree"),
              br(), br(),
              conditionalPanel(
                condition = "output.SHOW_PLOT",
                plotOutput("PANEL_PLOT")
              )
            )
          )
        )

      )

    )

  )

)

ui
