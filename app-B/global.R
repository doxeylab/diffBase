#-------------------------------------------------------------------------------
# global.R
# Last modified: 2020-03-05 22:59:43 (CET)
# BJM Tremblay

msg <- function(...) {
  time <- format(as.POSIXlt(Sys.time(), tz = "America/Toronto"))
  message(paste0(c(time, ...), collapse = " "))
}
msg("Loading packages")

msg("  shiny")
library(shiny)
msg("  shinythemes")
library(shinythemes)
msg("  ape")
library(ape)
msg("  ggtree")
suppressMessages(suppressPackageStartupMessages(library(ggtree)))
msg("  tidytree")
suppressPackageStartupMessages(library(tidytree))
msg("  magrittr")
library(magrittr)

CONFIGS <- readr::read_lines("diffBaseConfig.txt") %>%
  Filter(function(x) x != "", .) %>%
  Filter(function(x) !grepl("^\\s+$", x), .) %>%
  paste0(collapse = ",") %>%
  paste0("list(", ., ")") %>%
  parse(text = .) %>%
  eval()

if (!is.null(CONFIGS$UseGmail) && CONFIGS$UseGmail) {
  Sys.setenv(GMAILR_APP = CONFIGS$GmailCredentials)
  gmailr::gm_auth_configure(CONFIGS$GmailCredentials)
  gmailr::gm_auth(email = TRUE, cache = CONFIGS$GmailCache)
}

dir.create("community", showWarnings = FALSE)
dir.create("queries", showWarnings = FALSE)

#-------------------------------------------------------------------------------
msg("Loading functions")

noNA <- function(x) x[!is.na(x)]

plot_tree <- function(tip = NULL) {
  out <- TREE_PLOT
  if (!is.null(tip)) {
    tmp_node <- noNA(TREE_TBL$node[TREE_TBL$label == tip])
    out %<+% TREE_TBL +
      geom_point2(aes(subset = (node == tmp_node)),
                  shape = 21, size = 5, fill = "red")
  } else {
    out
  }
}

make_type_info <- function() {
  tagList(
    actionLink("BUTTON_GO_BACK_TO_WELCOME", "Go back"), br(),
    br(),
    htmlOutput("PANEL_LEFT_CURRENT_TYPE"),
    br(),
    htmlOutput("PANEL_LEFT_CURRENT_SUBTYPE"),
    br(),
    textOutput("CURRENT_ACCESSION"),
    selectInput(
      "SUBTYPE_SELECTOR",
      label = "",
      choices = names(SEQ_NAMES_LIST$A1),
      selected = "A1.1",
      width = "85px"
    ),
    br(),
    downloadLink("DOWNLOAD_TYPE", "Download all subtype sequences")
  )
}

make_type_info_more <- function() {
  tagList(
    tags$b("Information:"), br(),
    br(),
    tags$b("References:"), br(),
    br(),
    htmlOutput("PANEL_TOP_RIGHT_CURRENT_SUBTYPE"),
    br(),
    verbatimTextOutput("PANEL_TOP_RIGHT_CURRENT_SUBTYPE_SEQUENCE")
  )
}

clean_query <- function(query) {
  gsub("\\s+", "", toupper(query))
}

check_query <- function(x) {
  if (grepl("^>", x)) {
    return(list(
      status = FALSE,
      reason = "Please make sure query is not fasta-formatted."
    ))
  }
  if (grepl("[^-ABCDEFGHIKLMNPQRSTUVWYZX*]", x)) {
    return(list(
      status = FALSE,
      reason = paste(
        "An unknown character was found. Only the following characters",
        "(as well as whitespace) are allowed: ABCDEFGHIKLMNPQRSTUVWYZX*-"
      )
    ))
  }
  if (nchar(x) < 10 || nchar(x) > 10000) {
    return(list(
      status = FALSE,
      reason = "Queries must be in between 10-10000 characters long."
    ))
  }
  list(status = TRUE)
}

make_blast_buttons <- function(subject) {
  make_single <- function(subj) {
    subj <- gsub(".", "_", subj, fixed = TRUE)
    as.character(actionLink(
      paste0("BLASTP_GOTO_", subj), "Go",
      onclick = 'Shiny.onInputChange(\"BLASTP_GOTO\", this.id)'
    ))
  }
  vapply(subject, make_single, character(1))
}

run_blast <- function(query, evalue = 1) {
  query <- clean_query(query)
  check <- check_query(query)
  if (!check$status) {
    showModal(modalDialog(title = "Error", check$reason))
  } else {
    msg("Running blastp job")
    d <- as.integer(Sys.time())
    f <- paste0("queries/", d, ".fa")
    o <- paste0("queries/", d, ".res")
    cat(c(">QUERY", query), file = f, sep = "\n")
    cmd <- paste(
      "blastp -query", f, "-db blastdb/ALL-sequences.fa", "-out", o,
      "-outfmt '6 qseqid sseqid evalue mismatch pident length'",
      CONFIGS$BlastpParameters
    )
    if (system(cmd)) {
      msg("Blastp failed")
      showModal(modalDialog(
        title = "BLAST Error", "BLAST failed. Please contact the admin."
      ))
    } else {
      msg("Blastp successful")
      res <- suppressMessages(readr::read_tsv(o, col_names = FALSE))
      msg("Number of hits:", nrow(res))
      colnames(res) <- c("qseqid", "Match", "E-Value", "# of Mismatches", "% Identity", "Coverage")
      res$`Match Coverage %` <- round(
        100 * (res$Coverage / nchar(SEQS_ALL[res$Match])), 1
      )
      res$`Go To Toxin Page` <- make_blast_buttons(res$Match)
      as.data.frame(res)[, -1]
    }
  }
}

clean_email <- function(x) {
  gsub("\\s+", "", x)
}

check_email <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x),
        ignore.case = TRUE)
}

#-------------------------------------------------------------------------------
msg("Loading data")

TREE <- readRDS("data/tree.RDS")

TREE_TBL <- as_tibble(TREE)

SEQ_NAMES_LIST <- list(
  A = readRDS("data/A-names.RDS"),
  B = readRDS("data/B-names.RDS"),
  C = readRDS("data/C-names.RDS"),
  D = readRDS("data/D-names.RDS"),
  E = readRDS("data/E-names.RDS"),
  F = readRDS("data/F-names.RDS"),
  G = readRDS("data/G-names.RDS"),
  H = readRDS("data/H-names.RDS"),
  I = readRDS("data/I-names.RDS"),
  J = readRDS("data/J-names.RDS"),
  K = readRDS("data/K-names.RDS"),
  L = readRDS("data/L-names.RDS")
)
SEQ_NAMES_ALL <- do.call(c, SEQ_NAMES_LIST)
names(SEQ_NAMES_ALL) <- gsub("^[A-Z][.]", "", names(SEQ_NAMES_ALL))
names(SEQ_NAMES_ALL) <- gsub("^[A-Z][1-2][.]", "", names(SEQ_NAMES_ALL))

SEQS_ALL <- gsub("-", "", readRDS("data/ALL-sequences.RDS"), fixed = TRUE)
# SEQS_ALL <- readRDS("data/ALL-sequences.RDS")
# SEQS_ALL_AA <- readRDS("data/ALL-sequences-AAStringSet.RDS)

# SEQS_LIST <- list(  # Note: these are AAStringSet type
#   A1 = readRDS("data/A1-sequences.RDS"),
#   A2 = readRDS("data/A2-sequences.RDS"),
#   B = readRDS("data/B-sequences.RDS"),
#   C = readRDS("data/C-sequences.RDS"),
#   D = readRDS("data/D-sequences.RDS"),
#   E = readRDS("data/E-sequences.RDS"),
#   F = readRDS("data/F-sequences.RDS"),
#   G = readRDS("data/G-sequences.RDS"),
#   H = readRDS("data/H-sequences.RDS"),
#   I = readRDS("data/I-sequences.RDS"),
#   J = readRDS("data/J-sequences.RDS"),
#   K = readRDS("data/K-sequences.RDS"),
#   L = readRDS("data/L-sequences.RDS")
# )

ALL_TYPES <- c(
  "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"
)

for (a_t in ALL_TYPES) {
  if (!file.exists(paste0("downloads/", a_t, "-sequences.fa"))) {
    Biostrings::writeXStringSet(
      Biostrings::AAStringSet(
        SEQS_ALL[names(SEQ_NAMES_LIST[[a_t]])]
      ),
      paste0("downloads/", a_t, "-sequences.fa")
    )
  }
}

ALL_TYPES_SUBTYPES <- lapply(SEQ_NAMES_LIST, function(x) names(x))

TREE2 <- groupOTU(TREE, ALL_TYPES_SUBTYPES)

clades <- data.frame(Label = ALL_TYPES, Node = NA)
for (i in seq_len(nrow(clades))) {
  clades$Node[i] <- MRCA(TREE2, ALL_TYPES_SUBTYPES[[i]])
}

TREE_PLOT <- ggtree(TREE2) +
  layout_dendrogram()
for (i in seq_len(nrow(clades))) {
  TREE_PLOT <- TREE_PLOT +
    geom_cladelabel(
      node = clades$Node[i], label = as.character(clades$Label)[i],
      align = TRUE, offset = -0.148, offset.text = -0.005,
      hjust = 0.4
    )
}
