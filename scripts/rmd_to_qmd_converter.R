args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript rmd_to_qmd_combined.R input.Rmd output.qmd")
input_file <- args[1]
output_file <- args[2]

# --- First pass: Remove opts_chunk block and collect options ---
rmd_lines <- readLines(input_file)
clean_lines <- character(0)
knitr_opts <- character(0)
in_opts_chunk <- FALSE

for (line in rmd_lines) {
  if (grepl('^opts_chunk\\[\\["set"\\]\\]\\(', line)) {
    in_opts_chunk <- TRUE
    line <- sub('^opts_chunk\\[\\["set"\\]\\]\\(', '', line)
    if (grepl('\\)$', line)) {
      line <- sub('\\)$', '', line)
      in_opts_chunk <- FALSE
    }
    knitr_opts <- c(knitr_opts, line)
    next
  }
  if (in_opts_chunk) {
    if (grepl('\\)$', line)) {
      line <- sub('\\)$', '', line)
      in_opts_chunk <- FALSE
    }
    knitr_opts <- c(knitr_opts, line)
    next
  }
  clean_lines <- c(clean_lines, line)
}

# --- Second pass: Convert chunk options to Quarto #| style ---
qmd_lines <- character(0)
in_chunk <- FALSE
chunk_header <- NULL
chunk_body <- character(0)

parse_chunk_opts <- function(header) {
  header <- trimws(header)
  if (header == "") return(character(0))
  parts <- unlist(strsplit(header, ", *"))
  if (length(parts) > 0 && !grepl("=", parts[1])) {
    parts <- parts[-1]
  }
  opts_lines <- character(0)
  for (opt in parts) {
    kv <- unlist(strsplit(opt, '='))
    if (length(kv) == 2) {
      key <- trimws(kv[1])
      value <- trimws(kv[2])
      value <- gsub('TRUE', 'true', value)
      value <- gsub('FALSE', 'false', value)
      value <- gsub('T', 'true', value)
      value <- gsub('F', 'false', value)
      key <- gsub('fig.width', 'fig-width', key)
      key <- gsub('fig.height', 'fig-height', key)
      key <- gsub('fig.align', 'fig-align', key)
      opts_lines <- c(opts_lines, paste0('#| ', key, ': ', value))
    }
  }
  opts_lines
}

for (i in seq_along(clean_lines)) {
  line <- clean_lines[i]
  if (!in_chunk) {
    if (grepl('^```\\{r', line)) {
      in_chunk <- TRUE
      chunk_body <- character(0)
      chunk_header <- sub('^```\\{r', '', line)
      chunk_header <- sub('\\}$', '', chunk_header)
      chunk_header <- trimws(chunk_header)
      if (chunk_header != "" && !grepl('=', unlist(strsplit(chunk_header, ', *'))[1])) {
        label <- unlist(strsplit(chunk_header, ', *'))[1]
        qmd_lines <- c(qmd_lines, paste0('```{r ', label, '}'))
      } else {
        qmd_lines <- c(qmd_lines, '```{r}')
      }
      opts_lines <- parse_chunk_opts(chunk_header)
      if (length(opts_lines) > 0) {
        qmd_lines <- c(qmd_lines, opts_lines)
      }
    } else {
      qmd_lines <- c(qmd_lines, line)
    }
  } else {
    if (grepl('^```$', line)) {
      qmd_lines <- c(qmd_lines, chunk_body, '```')
      in_chunk <- FALSE
      chunk_header <- NULL
      chunk_body <- character(0)
    } else {
      chunk_body <- c(chunk_body, line)
    }
  }
}

# --- Insert knitr options into YAML ---
if (length(knitr_opts) > 0) {
  knitr_opts <- gsub(',$', '', trimws(knitr_opts))
  knitr_opts <- knitr_opts[knitr_opts != ""]
  yaml_opts <- character(0)
  for (opt in knitr_opts) {
    kv <- unlist(strsplit(opt, '='))
    if (length(kv) == 2) {
      key <- trimws(kv[1])
      value <- trimws(kv[2])
      value <- gsub('TRUE', 'true', value)
      value <- gsub('FALSE', 'false', value)
      value <- gsub('T', 'true', value)
      value <- gsub('F', 'false', value)
      key <- gsub('fig.width', 'fig-width', key)
      key <- gsub('fig.height', 'fig-height', key)
      key <- gsub('fig.align', 'fig-align', key)
      yaml_opts <- c(yaml_opts, paste0('    ', key, ': ', value))
    }
  }
  yaml_start <- which(grepl("^---$", qmd_lines))[1]
  yaml_end <- which(grepl("^---$", qmd_lines))
  yaml_end <- yaml_end[yaml_end > yaml_start][1]
  if (!is.na(yaml_start) && !is.na(yaml_end)) {
    qmd_lines <- append(qmd_lines, c("knitr:", "  opts_chunk:", yaml_opts), after = yaml_end - 1)
  }
}

writeLines(qmd_lines, output_file)
cat("Combined conversion complete!\n") 