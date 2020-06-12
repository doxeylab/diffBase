#-------------------------------------------------------------------------------
# global.R
# Last modified: 2020-06-12 10:43:56 (CEST)
# BJM Tremblay

LAST_UPDATE_DATE <- function() "2020-06-12"

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
    # htmlOutput("PANEL_LEFT_CURRENT_SUBTYPE"),
    htmlOutput("CURRENT_ACCESSION"),
    selectInput(
      "SUBTYPE_SELECTOR",
      label = "",
      choices = unname(SEQ_NAMES_LIST$B1),
      selected = "B1.1",
      width = "145px"
    ),
    # br(),
    downloadLink("DOWNLOAD_TYPE", "Download all subtype sequences")
  )
}

show_metadata <- function(subtype) {
  ACCs <- META2ACC$Acc[META2ACC$Subtype == subtype]
  ACCs <- vapply(strsplit(ACCs, ".", fixed = TRUE), function(x) x[1], character(1))
  ACCs <- unique(ACCs)
  if (any(!ACCs %in% names(METADATA))) {
    ACCs2 <- ACCs[!ACCs %in% names(METADATA)]
    ACCs <- ACCs[ACCs %in% names(METADATA)]
    ACCs2 <- tibble(
      Id = NA, Source = NA,
      `Nucleotide Accession` = NA,
      Start = NA, Stop = NA, Strand = NA, Protein = ACCs2,
      `Protein Name` = NA, Organism = NA, Strain = NA, Assembly = NA
    )
    if (!length(ACCs))
      out <- ACCs2
    else
      out <- rbind(do.call(rbind, METADATA[ACCs]), ACCs2)
  } else {
    x <- METADATA[ACCs]
    if (length(x) && sum(vapply(x, nrow, integer(1)))) {
      out <- do.call(rbind, x)
    } else {
      out <- tibble(
        Id = NA, Source = NA,
        `Nucleotide Accession` = NA,
        Start = NA, Stop = NA, Strand = NA, Protein = names(x),
        `Protein Name` = NA, Organism = NA, Strain = NA, Assembly = NA
      )
    }
  }
  out[!is.na(out$Protein), ]
}

make_type_info_more <- function() {
  tagList(
    htmlOutput("PANEL_TOP_RIGHT_CURRENT_SUBTYPE"),
    br(),
    verbatimTextOutput("PANEL_TOP_RIGHT_CURRENT_SUBTYPE_SEQUENCE"),
    br(),
    tags$b("Sequence accessions:"),
    DT::dataTableOutput("PANEL_TOP_RIGHT_METADATA"),
    br(),
    downloadLink("DOWNLOAD_METADATA", "Download sequence accessions table")
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
    o2 <- paste0("queries/", d, "-RAW.res")
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
      if (nrow(res) == 0) return(NULL)
      system(paste(
        "blastp -query", f, "-db blastdb/ALL-sequences.fa", "-out",
        o2, CONFIGS$BlastpParameters
      ))
      res2 <- suppressMessages(readr::read_lines(o2))
      colnames(res) <- c("qseqid", "Match", "E-Value", "# of Mismatches", "% Identity", "Coverage")
      res$`Match Coverage %` <- round(
        100 * (res$Coverage / nchar(as.character(SEQS_ALL)[res$Match])), 1
      )
      res$`Go To Toxin Page` <- make_blast_buttons(res$Match)
      list(as.data.frame(res)[order(res[[5]], decreasing = TRUE), -1], res2)
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

SEQ_NAMES_LIST <- readRDS("data/ALL-names.RDS")
SEQ_NAMES_ALL <- SEQ_NAMES_LIST[[1]]
for (i in seq_along(SEQ_NAMES_LIST)[-1]) {
  SEQ_NAMES_ALL <- c(SEQ_NAMES_ALL, SEQ_NAMES_LIST[[i]])
}

SEQS_ALL <- readRDS("data/ALL-sequences.RDS")

ALL_TYPES <- names(SEQ_NAMES_LIST)

for (a_t in ALL_TYPES) {
  if (!file.exists(paste0("downloads/", a_t, "-sequences.fa"))) {
    Biostrings::writeXStringSet(
      Biostrings::AAStringSet(
        SEQS_ALL[[a_t]]
      ),
      paste0("downloads/", a_t, "-sequences.fa")
    )
  }
}

SEQS_ALL_AA <- SEQS_ALL
SEQS_ALL <- SEQS_ALL[[1]]
for (i in seq_along(SEQS_ALL_AA)[-1]) SEQS_ALL <- c(SEQS_ALL, SEQS_ALL_AA[[i]])

if (!file.exists("downloads/ALL-sequences.fa"))
  Biostrings::writeXStringSet(SEQS_ALL, "downloads/ALL-sequences.fa")

ALL_TYPES_SUBTYPES <- lapply(SEQ_NAMES_LIST, function(x) unname(x))

TREE2 <- groupOTU(TREE, ALL_TYPES_SUBTYPES)

clades <- data.frame(Label = ALL_TYPES, Node = NA, Size = NA_integer_, stringsAsFactors = FALSE)
for (i in seq_len(nrow(clades))) {
  clades$Node[i] <- MRCA(TREE2, ALL_TYPES_SUBTYPES[[i]])
  clades$Size[i] <- length(SEQ_NAMES_LIST[[i]])
}

TREE_PLOT <- ggtree(TREE2, branch.length = "none") +
  # geom_tiplab(size = 1)
  layout_dendrogram()
for (i in seq_len(nrow(clades))) {
  if (clades$Size[i] <= 1) next
  if (nchar(clades$Label[i]) > 3) next
  TREE_PLOT <- TREE_PLOT +
    geom_cladelabel(
      # node = clades$Node[i], label = as.character(clades$Label)[i],
      # node = clades$Node[i], label = gsub("^[A-Z]", "", as.character(clades$Label)[i]),
      node = clades$Node[i],
      label = as.character(clades$Label)[i],
      align = TRUE,
      # offset = 15 + i
      offset = -73.5,
      offset.text = -2.8,
      hjust = 0.5
      # offset = -0.148,
      # offset.text = -0.005,
      # hjust = 0.4
    )
}

METADATA <- readRDS("data/metadata.RDS")
METADATA <- lapply(METADATA, function(x) x[x$Source != "PAT", ])
METADATA <- METADATA[unname(sapply(METADATA, nrow)) > 0]
METADATA_ALL <- do.call(rbind, METADATA)

META2ACC <- readRDS("data/metadata2acc.RDS")
META2ACC$Subtype <- SEQ_NAMES_ALL[META2ACC$Subtype]

REP_SEQS <- suppressMessages(readr::read_tsv("data/repseqs.tsv"))
REP_SEQS$Sequence <- SEQ_NAMES_ALL[REP_SEQS$Rep_identifier]

if (!file.exists("downloads/REPRESENTATIVE-sequences.fa")) {
  Biostrings::writeXStringSet(
    SEQS_ALL[REP_SEQS$Sequence], "downloads/REPRESENTATIVE-sequences.fa"
  )
}

NUMBER_OF_SEQUENCES <- function() length(SEQ_NAMES_ALL)
NUMBER_OF_STRAINS <- function() length(unique(METADATA_ALL$Strain))
NUMBER_OF_NCBI_IDs <- function() length(METADATA)

GET_MEMBER_COUNT <- function(x) {
  x <- length(SEQ_NAMES_LIST[[x]])
  paste(x, if (x > 1) "members" else "member")
}
