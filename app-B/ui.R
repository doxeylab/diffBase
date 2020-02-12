#-------------------------------------------------------------------------------
# ui.R
# Last modified: 2020-02-12 21:42:26 (CET)
# BJM Tremblay

msg("Loading UI")
ui <- function(request) fluidPage(

  theme = shinytheme("simplex"),

  title = "Diff-Base: Toxin B",

  br(), br(), br(),

  navbarPage(
    tags$li(
      class = "dropdown",
      tags$style(".navbar {min-height:75px}")
    ),
    title = div(
      img(src = "DiffBaseLogo.png", height = 53, width = 100),
      HTML("<a href='http://35.193.113.52:3838/'>Diff-Base</a>")
    ),
    id = "NAVBARPAGE",
    position = "fixed-top",
    collapsible = TRUE,
    footer = tags$footer(
      tags$div(HTML(
        "<p>For questions or concerns, email <a href='mailto:diffbaseserver@gmail.com'>diffbaseserver@gmail.com</a>."
      )),
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
              h4("Toxinotypes (TcdB)"),
              actionLink("BUTTON_A1", "A1:"), "57 members", br(),
              actionLink("BUTTON_A2", "A2:"), "12 members", br(),
              actionLink("BUTTON_B", "B:"), "2 members", br(),
              actionLink("BUTTON_C", "C:"), "11 members", br(),
              actionLink("BUTTON_D", "D:"), "7 members", br(),
              actionLink("BUTTON_E", "E:"), "14 members", br(),
              actionLink("BUTTON_F", "F:"), "2 members", br(),
              actionLink("BUTTON_G", "G:"), "1 member", br(),
              actionLink("BUTTON_H", "H:"), "11 members", br(),
              actionLink("BUTTON_I", "I:"), "1 member", br(),
              actionLink("BUTTON_J", "J:"), "1 member", br(),
              actionLink("BUTTON_K", "K:"), "1 member", br(),
              actionLink("BUTTON_L", "L:"), "5 members", br(),
              br(),
              downloadLink(
                "DOWNLOAD_ALL", "Download all toxin B sequences"
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
          )

        ),
        wellPanel(
          tagList(
            h4("Database information"),
            tags$b("Last updated data:"), br(),
            tags$b("Number of sequences:"), br(),
            tags$b("Number of bacterial strains:"), br(),
            tags$b("Number of NCBI IDs:"), br(),
            tags$b("Number of publications:"), br()
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
              plotOutput("PANEL_PLOT")
            )
          )
        )

      )

    )

  )

)

ui
