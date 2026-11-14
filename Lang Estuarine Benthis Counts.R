install.packages("readr")
library(readr)

raw_path <- file.choose()
raw_hash <- unname(tools::md5sum(raw_path))
out_dir <- file.path(dirname(raw_path), "week3_cleaning_output")
dir.create(out_dir, showWarnings = FALSE)
stamp <- format(Sys.Date(), "%Y%m%d")
topic <- "BENTHIC_COUNT_ESTUARINE"
output_path <- function(label, ext = "csv") {
  file.path(out_dir, paste0(topic, "_", label, "_", stamp, ".", ext))
}

# Read identifiers as text; preserve original values for comparison.
raw <- read_csv(raw_path, col_types = cols(.default = col_character()),
                na = character(), trim_ws = FALSE, show_col_types = FALSE)
expected <- c("PUBLICATION_DATE", "TAXA_ID", "UID", "SITE_ID", "VISIT_NO",
              "YEAR", "DATE_COL", "STUDY", "PSTL_CODE", "NCCA_REG", "INDEX_AP",
              "INDEX_NCCA15", "IS_DISTINCT", "TOTAL", "TARGET_TAXON", "APHIA_ID",
              "ITISTSN", "HABITAT", "PHYLUM", "CLASS", "ORDER", "FAMILY")
stopifnot(setequal(names(raw), expected), nrow(problems(raw)) == 0)
clean <- as.data.frame(raw)
log <- data.frame(Column = character(), Issue = character(),
                  `Action Taken` = character(), Justification = character(),
                  Impact = character(), check.names = FALSE)
add_log <- function(column, issue, action, justification, impact) {
  log[nrow(log) + 1, ] <<- list(column, issue, action, justification, impact)
}
add_log("All", "Raw-data preservation", "Read raw CSV without writing to it; record MD5 checksum.",
        "Keep original observations reproducible and recoverable.",
        paste(nrow(raw), "rows and", ncol(raw), "original columns imported."))

# Trim surrounding spaces and represent blank fields consistently.
for (v in names(clean)) {
  z <- trimws(clean[[v]])
  trimmed <- sum(z != clean[[v]])
  z[z == ""] <- NA_character_
  clean[[v]] <- z
  if (trimmed > 0) add_log(v, "Surrounding whitespace", paste("Trimmed", trimmed, "values."),
                           "Remove formatting differences, not biological information.", "No rows removed.")
  if (anyNA(z)) add_log(v, paste(sum(is.na(z)), "blank fields"),
                        "Represent as NA in R and export as NA; retain rows; no imputation.",
                        "Unknown classifications/IDs cannot be reliably invented. Blank index flags are not assumed to mean N.",
                        "Preserves biological records; missingness remains explicit.")
}

# Normalize label case without altering taxonomic spelling or study membership.
upper_cols <- c("PSTL_CODE", "INDEX_AP", "INDEX_NCCA15", "IS_DISTINCT",
                "HABITAT", "TARGET_TAXON", "PHYLUM", "CLASS", "ORDER", "FAMILY")
case_changes <- 0L
for (v in upper_cols) {
  z <- toupper(clean[[v]])
  case_changes <- case_changes + sum(z != clean[[v]], na.rm = TRUE)
  clean[[v]] <- z
}
add_log(paste(upper_cols, collapse = ", "), "Text case inspection",
        paste("Applied uppercase to code/taxonomic-label fields;", case_changes, "values changed."),
        "Keep codes consistent; preserve supplied spelling and established study/region labels.", "No taxonomic reconciliation attempted.")

# Parse dates explicitly. %b uses English month abbreviations in this CSV.
old_time <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")
clean$PUBLICATION_DATE <- as.Date(clean$PUBLICATION_DATE, format = "%m/%d/%Y")
clean$DATE_COL <- as.Date(clean$DATE_COL, format = "%d-%b-%y")
Sys.setlocale("LC_TIME", old_time)
stopifnot(!anyNA(clean$PUBLICATION_DATE), !anyNA(clean$DATE_COL))
add_log("PUBLICATION_DATE; DATE_COL", "Non-ISO date strings",
        "Convert to Date and export as YYYY-MM-DD.", "Unambiguous dates facilitate sorting and later joins.",
        paste(nrow(clean), "values standardized in each date column; no parsing failures."))

# Convert genuine numeric quantities, but keep taxon and site-visit IDs as text.
for (v in c("TOTAL", "VISIT_NO", "YEAR")) {
  z <- suppressWarnings(as.numeric(clean[[v]]))
  if (anyNA(z) || any(!is.finite(z)) || any(z != floor(z))) {
    stop(paste("Review missing, nonnumeric, or noninteger values in", v))
  }
  clean[[v]] <- as.integer(z)
}
stopifnot(all(clean$TOTAL >= 0), all(clean$VISIT_NO %in% c(1L, 2L)),
          all(clean$YEAR == 2015L), all(format(clean$DATE_COL, "%Y") == "2015"),
          all(clean$IS_DISTINCT %in% c("Y", "N")),
          all(is.na(clean$INDEX_AP) | clean$INDEX_AP %in% c("Y", "N")),
          all(is.na(clean$INDEX_NCCA15) | clean$INDEX_NCCA15 %in% c("Y", "N")),
          !anyNA(clean$UID), !anyNA(clean$TAXA_ID), !anyNA(clean$SITE_ID),
          !anyNA(clean$TARGET_TAXON), all(clean$PSTL_CODE %in% state.abb),
          all(clean$NCCA_REG %in% c("Northeast", "Southeast", "Gulf", "West")),
          all(clean$STUDY %in% c("Estuarine", "Albemarle_Pamlico")),
          all(clean$HABITAT == "EST"))
add_log("TOTAL; VISIT_NO; YEAR; identifiers; codes", "Types and impossible-value inspection",
        "Convert count, visit, and year to integers; retain IDs as text; validate expected categories and dates. Stop on failed checks.",
        "Counts are nonnegative whole individuals; identifiers are labels rather than measurements.",
        "No unit conversion needed; TOTAL remains individuals per taxon record, not density.")

# Audit duplicates without guessing how conflicting observations should be resolved.
full_duplicates <- sum(duplicated(clean))
key_duplicates <- sum(duplicated(clean[c("UID", "TAXA_ID")]))
if (full_duplicates > 0 || key_duplicates > 0) {
  stop("Duplicates detected: investigate before deciding whether to remove or combine records.")
}
add_log("All; UID + TAXA_ID", "Duplicate inspection",
        paste(full_duplicates, "complete duplicate rows;", key_duplicates, "repeated site-visit/taxon keys. None removed."),
        "Repeated visits to the same site are not duplicate observations.", "All site-visits preserved.")

# IQR screening applies to organism counts, NOT numeric-looking identifiers.
quartiles <- quantile(clean$TOTAL, c(0.25, 0.75), type = 7)
iqr <- unname(quartiles[2] - quartiles[1])
lower <- unname(quartiles[1] - 1.5 * iqr)
upper <- unname(quartiles[2] + 1.5 * iqr)
clean$TOTAL_IQR_FLAG <- clean$TOTAL < lower | clean$TOTAL > upper
clean$ZERO_COUNT_FLAG <- clean$TOTAL == 0L
clean$NO_MACROINVERTEBRATES_FLAG <-
  clean$TARGET_TAXON == "NO_MACROINVERTEBRATES_IN_SAMPLE"
add_log("TOTAL; TOTAL_IQR_FLAG", "Statistical outlier screening",
        paste("Q1 =", quartiles[1], "; Q3 =", quartiles[2], "; IQR =", iqr,
              "; fences =", lower, "and", upper, "; flagged", sum(clean$TOTAL_IQR_FLAG), "records. Retain all."),
        "Pooled count distributions may be strongly skewed. A statistical flag alone does not establish a measurement error.",
        "Preserves high-abundance observations for review; flags are not evidence that records are validated or erroneous.")
add_log("TOTAL; TARGET_TAXON; zero/sample flags", "Zero counts and no-organism placeholder",
        paste("Flag and retain", sum(clean$ZERO_COUNT_FLAG), "zero-count rows and",
              sum(clean$NO_MACROINVERTEBRATES_FLAG), "no-macroinvertebrate placeholders."),
        "A no-organism sample is informative, but its placeholder must not count as an observed taxon.",
        "Later richness calculations must exclude placeholders and zero counts while retaining these sample visits in the sample frame.")
add_log("IS_DISTINCT; STUDY; index flags; UID", "Analysis eligibility and repeated visits",
        "Retain all studies, index statuses, distinctness flags and visits; do not compute final richness in this cleaning step.",
        "Study inclusion and distinctness rules need documented methodological decisions; repeated visits are not independent sites.",
        "Biological table is cleaned, but final environmental analysis still requires companion files and eligibility decisions.")

# Missingness report: original blanks versus missing values after cleaning.
missing_report <- data.frame(
  Column = names(raw),
  Raw_blank_fields = vapply(raw, function(z) sum(trimws(z) == ""), integer(1)),
  Clean_missing_fields = vapply(clean[names(raw)], function(z) sum(is.na(z)), integer(1))
)
missing_report$Clean_missing_percent <- round(100 * missing_report$Clean_missing_fields / nrow(clean), 2)

# Final protections: no dropped/reordered rows and no changed organism counts.
stopifnot(nrow(clean) == nrow(raw),
          identical(clean$UID, raw$UID), identical(clean$TAXA_ID, raw$TAXA_ID),
          identical(clean$TOTAL, as.integer(raw$TOTAL)),
          identical(unname(tools::md5sum(raw_path)), raw_hash))
add_log("All", "Final validation", "Verify original row count, identifier order, counts and raw-file MD5 checksum.",
        "Guard against accidental deletion, reordering, count changes or raw-file overwrite.",
        paste(nrow(clean), "rows exported;", ncol(clean), "columns including three new logical audit flags."))

write_csv(clean, output_path("CLEANED"), na = "NA")
write_csv(log, output_path("CLEANING_LOG"), na = "NA")
write_csv(missing_report, output_path("MISSINGNESS_AUDIT"), na = "NA")
write_csv(clean[clean$TOTAL_IQR_FLAG | clean$ZERO_COUNT_FLAG, ],
          output_path("FLAGGED_RECORDS"), na = "NA")
writeLines(c(paste("Raw file:", basename(raw_path)), paste("MD5:", raw_hash)),
           output_path("RAW_CHECKSUM", "txt"))
capture.output(sessionInfo(), file = output_path("R_SESSION", "txt"))

# Required PDF cleaning log: five columns, one readable entry per landscape page.
# Uses grid, included with R; no additional document packages needed.
pdf(output_path("CLEANING_LOG", "pdf"), width = 11.7, height = 8.3)
widths <- c(.15, .17, .25, .23, .20)
lefts <- c(0, head(cumsum(widths), -1))
wrap <- function(z, n) paste(strwrap(z, width = n), collapse = "\n")
for (i in seq_len(nrow(log))) {
  grid::grid.newpage()
  grid::grid.text("NCCA 2015 Estuarine Benthic Counts: Cleaning Log", y = .95,
                  gp = grid::gpar(fontsize = 15, fontface = "bold"))
  grid::grid.text(paste("Entry", i, "of", nrow(log), "| Run date:", stamp), y = .90,
                  gp = grid::gpar(fontsize = 10))
  for (j in seq_along(widths)) {
    x <- .025 + .95 * lefts[j]
    w <- .95 * widths[j]
    grid::grid.rect(x = x, y = .84, width = w, height = .07,
                    just = c("left", "top"), gp = grid::gpar(fill = "grey90"))
    grid::grid.text(names(log)[j], x = x + .008, y = .82,
                    just = c("left", "top"), gp = grid::gpar(fontsize = 10, fontface = "bold"))
    grid::grid.text(wrap(as.character(log[i, j]), floor(widths[j] * 135)),
                    x = x + .008, y = .74, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10, lineheight = 1.2))
  }
}
dev.off()
message("Finished. Outputs saved in: ", out_dir)
message("Required submissions: CLEANED CSV and CLEANING_LOG PDF.")
message("Review the log and flagged records before submission.")