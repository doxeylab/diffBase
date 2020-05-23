#-------------------------------------------------------------------------------
# ui.R
# Last modified: 2020-05-23 16:42:22 (CEST)
# BJM Tremblay

msg("Loading UI")
ui <- function(request) fluidPage(

  theme = shinytheme("simplex"),

  title = "Diff-Base: Toxin A",

  br(), br(), br(),

  navbarPage(
    tags$li(
      class = "dropdown",
      tags$style(".navbar {min-height:75px}")
    ),
    title = div(
      img(src = "DiffBaseLogo.png", height = 53, width = 100),
      HTML(paste0(
        "<a href='", CONFIGS$URL, "'>Diff-Base</a>"
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
    tags$div(HTML("<h2>Welcome to the <i>C. difficile</i> toxin A database</h2>")),
    br(), br(),

    fluidRow(

      column(3, # PANEL_LEFT
        wellPanel(

          conditionalPanel(
            condition = "output.CURRENT_PAGE == 'WELCOME'",
            tagList(
              h4("Toxin Groups (TcdA)"),
              actionLink("BUTTON_A1", "A1:"), GET_MEMBER_COUNT("A1"), br(),
              actionLink("BUTTON_A2", "A2:"), GET_MEMBER_COUNT("A2"), br(),
              actionLink("BUTTON_A3", "A3:"), GET_MEMBER_COUNT("A3"), br(),
              actionLink("BUTTON_A4", "A4:"), GET_MEMBER_COUNT("A4"), br(),
              actionLink("BUTTON_A5", "A5:"), GET_MEMBER_COUNT("A5"), br(),
              actionLink("BUTTON_A6", "A6:"), GET_MEMBER_COUNT("A6"), br(),
              actionLink("BUTTON_A7", "A7:"), GET_MEMBER_COUNT("A7"), br(),
              actionLink("BUTTON_sordellii_TcsH", "sordellii_TcsH:"), GET_MEMBER_COUNT("sordellii_TcsH"), br(),
              br(),
              downloadLink(
                "DOWNLOAD_ALL", "Download all toxin A sequences"
              ),
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
            # br(),
            # HTML(paste("<b>Number of NCBI IDs:</b>", NUMBER_OF_NCBI_IDs()))
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
              downloadLink("BLASTP_DOWNLOAD", "Download all results"), br(), br(),
              DT::dataTableOutput("BLASTP_RES_TABLE")
            )

          )
        ),

        column(12, # PANEL_BOTTOM_RIGHT
          wellPanel(
            tagList(
              h4("Phylogenetic tree"),
              conditionalPanel(
                condition = "output.CURRENT_PAGE != 'WELCOME'",
                tagList(
                  actionLink("TREE_BUTTON", "Show/hide tree"),
                  br(), br()
                )
              ),
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
