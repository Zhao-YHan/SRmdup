#!/usr/bin/env Rscript

# Properties plotting for a sample
# Input: qc_dir/{stats.tsv, read_lengths.before.tsv, read_lengths.after.tsv,
#                read_coverage.before.tsv, read_coverage.after.tsv}
# Output: {out_prefix}_summary.png, {out_prefix}_summary.eps
# Example:
# Rscript Properties_plot.R --sample_label sample_name --qc_dir input_path --out_dir output_path --out_prefix output_name
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[i], "--")) stop("Invalid argument: ", x[i], call. = FALSE)
    k <- sub("^--", "", x[i])
    out[[k]] <- if (i < length(x) && !startsWith(x[i + 1L], "--")) x[i + 1L] else TRUE
    i <- i + if (isTRUE(out[[k]])) 1L else 2L
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
num_pair <- function(x, d = NULL) if (is.null(x)) d else as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
fread2 <- function(...) { x <- fread(..., data.table = TRUE); setDT(x); x }
pick <- function(dt, cand, label) {
  z <- intersect(cand, names(dt))
  if (!length(z)) stop("Missing column in ", label, ": ", paste(cand, collapse = ", "), call. = FALSE)
  z[1]
}
save_plot <- function(p, prefix, w = 12, h = 14) {
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(prefix, ".eps"), p, width = w, height = h, device = grDevices::cairo_ps, fallback_resolution = 600, bg = "white")
}

theme_qc <- function(size = 12) {
  theme_bw(size) + theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey20"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )
}

COL_STRAND <- c(forward = "#F8766D", reverse = "#00BFC4")
COL_CAT <- c(`both equal` = "#F8766D", `start equal` = "#7CAE00", `end equal` = "#00BFC4", `both ±1` = "#C77CFF")
COL_LEN <- c(total = "#E64B35", kept = "#4DBBD5", deleted = "#00BA38")
COL_COV <- c(Before = "#AFC6FF", After = "#C96B6B")

cli <- parse_args(commandArgs(trailingOnly = TRUE))
sample_label   <- cli$sample_label %||% "Sample"
qc_dir         <- cli$qc_dir %||% "."
out_dir        <- cli$out_dir %||% qc_dir
out_prefix     <- cli$out_prefix %||% "sample"
length_xlim    <- num_pair(cli$length_xlim, c(30, 140))
coverage_xlim  <- num_pair(cli$coverage_xlim, NULL)
coverage_y_max <- if (is.null(cli$coverage_y_max)) NULL else as.numeric(cli$coverage_y_max)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_scope <- function(x) {
  x <- toupper(as.character(x))
  x[x %in% c("FORWARD", "FLAG0", "0")] <- "forward"
  x[x %in% c("REVERSE", "FLAG16", "16")] <- "reverse"
  x[x == "ALL"] <- "ALL"
  x
}

read_frag <- function(path) {
  dt <- fread2(path)
  dt[, length := as.integer(real_len)]
  dt[!is.na(length) & length > 0, .(length)]
}

read_depth <- function(path, label) {
  dt <- fread2(path, header = FALSE)
  dt <- if (ncol(dt) >= 3) {
    dt[, .(chr = as.character(V1), pos = as.integer(V2), depth = as.numeric(V3))]
  } else {
    dt[, .(chr = NA_character_, pos = as.integer(V1), depth = as.numeric(V2))]
  }
  dt <- dt[!is.na(pos) & !is.na(depth) & depth > 0]
  setorder(dt, chr, pos)
  dt[, label := label]
}

read_stats <- function() {
  st <- fread2(file.path(qc_dir, "stats.tsv"))
  st[, scope2 := norm_scope(scope)]
  st[scope2 %in% c("forward", "reverse", "ALL")]
}

plot_qc_ab <- function(st) {
  gap <- 0.55
  sr <- st[scope2 %in% c("forward", "reverse")]
  ct <- rbindlist(list(
    sr[, .(stage = "before", stage_x = 1,       strand = scope2, count = total_input)],
    sr[, .(stage = "after",  stage_x = 1 + gap, strand = scope2, count = total_kept)]
  ))
  ct[, strand := factor(strand, levels = c("reverse", "forward"))]

  tot_before <- sum(sr$total_input)
  tot_after  <- sum(sr$total_kept)
  del_fwd    <- sr[scope2 == "forward", total_deleted]
  del_rev    <- sr[scope2 == "reverse", total_deleted]
  del_tot    <- sum(sr$total_deleted)

  ggplot(ct, aes(stage_x, count, fill = strand)) +
    geom_col(width = 0.28) +
    scale_x_continuous(breaks = c(1, 1 + gap), labels = c("before", "after"), limits = c(0.75, 1 + gap + 0.25)) +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = COL_STRAND, breaks = c("forward", "reverse")) +
    labs(
      title = paste0(sample_label, " — Reads retained by strands"),
      subtitle = paste0(
        sprintf("Total retained: %.1f%% (%s / %s)", 100 * tot_after / tot_before, comma(tot_after), comma(tot_before)),
        "\n",
        sprintf("Deleted: %s (forward), %s (reverse); total %s", comma(del_fwd), comma(del_rev), comma(del_tot))
      ),
      x = NULL, y = "Sequence reads", fill = "Strand"
    ) +
    theme_qc(12)
}

plot_qc_cat <- function(st) {
  pm1 <- pick(st, c("both_±1", "both_+1", "both_pm1", "both_plusminus1", "both_plusminus_1"), "stats.tsv")
  kr  <- intersect(c("kept_rate%", "kept_rate", "retained_rate", "retained_frac"), names(st))
  all_row <- st[scope2 == "ALL"]
  subtxt <- if (length(kr)) {
    x <- all_row[[kr[1]]]
    xn <- suppressWarnings(as.numeric(x))
    if (grepl("%$", as.character(x)[1])) sprintf("Total retained (ALL): %s", as.character(x)[1])
    else if (!is.na(xn) && xn <= 1) sprintf("Total retained (ALL): %.1f%%", 100 * xn)
    else sprintf("Total retained (ALL): %.1f%%", xn)
  } else {
    sprintf("Total retained (ALL): %.1f%%", 100 * all_row$total_kept / all_row$total_input)
  }

  m <- melt(
    st,
    id.vars = c("scope2", "total_input", "total_deleted", "total_kept"),
    measure.vars = c("both_equal", "start_equal", "end_equal", pm1),
    variable.name = "category",
    value.name = "n_deleted"
  )
  m[category == "both_equal", category := "both equal"]
  m[category == "start_equal", category := "start equal"]
  m[category == "end_equal",   category := "end equal"]
  m[category == pm1,            category := "both ±1"]
  m[, category := factor(category, levels = c("both equal", "start equal", "end equal", "both ±1"))]

  ggplot(m[n_deleted > 0], aes(scope2, n_deleted, fill = category)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = COL_CAT) +
    scale_y_log10(labels = comma) +
    labs(
      title = paste0(sample_label, " — Duplicate categories"),
      subtitle = subtxt,
      x = NULL, y = "Deleted reads (log scale)", fill = "Category"
    ) +
    theme_qc(12)
}

plot_fraglen <- function() {
  fb <- read_frag(file.path(qc_dir, "read_lengths.before.tsv"))[, .(status = "total", length)]
  fa <- read_frag(file.path(qc_dir, "read_lengths.after.tsv"))[, .(status = "kept",  length)]
  z  <- rbindlist(list(fb, fa))

  cnt <- z[, .N, by = .(status, length)]
  cnt <- cnt[, {
    g <- data.table(length = min(length):max(length))
    g <- merge(g, .SD, by = "length", all.x = TRUE)
    g[is.na(N), N := 0L]
    g
  }, by = status]
  cnt[, y := frollmean(as.numeric(N), n = 5L, align = "center", fill = NA_real_), by = status]
  cnt[is.na(y), y := as.numeric(N)]

  ymax <- max(cnt$y, na.rm = TRUE)
  wd <- dcast(cnt[, .(status, length, N)], length ~ status, value.var = "N", fill = 0)
  wd[, del := fifelse(total > 0, (total - kept) / total, NA_real_)]
  wd[, del_s := frollmean(del, n = 5L, align = "center", fill = NA_real_)]
  wd[is.na(del_s), del_s := del]

  rlim <- min(1, max(0.2, ceiling(max(wd$del_s, na.rm = TRUE) * 10) / 10))
  wd[, del2 := pmin(pmax(del_s, 0), rlim) * (ymax / rlim)]

  total_n <- nrow(fb)
  kept_n  <- nrow(fa)

  p <- ggplot() +
    geom_line(data = cnt, aes(length, y, color = status), linewidth = 0.8) +
    geom_line(data = wd, aes(length, del2, color = "deleted"), linetype = "dotted", linewidth = 0.95) +
    scale_color_manual(values = COL_LEN, breaks = c("total", "kept", "deleted"), labels = c("Total reads", "Kept reads", "Deleted proportion")) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid", "dotted"), linewidth = c(0.8, 0.8, 0.95)))) +
    scale_y_continuous(
      name = "Number of reads",
      labels = comma,
      sec.axis = sec_axis(~ . * (rlim / ymax), name = "Deleted proportion",
                          breaks = seq(0, rlim, by = if (rlim <= 0.2) 0.05 else 0.1),
                          labels = percent_format(accuracy = 1))
    ) +
    labs(title = paste0(sample_label, " — Read lengths"), x = "Read lengths (bp)", color = NULL) +
    theme_qc(12) +
    theme(legend.position = "bottom") +
    geom_text(
      aes(x = Inf, y = Inf,
          label = sprintf("Retained: %.1f%%  (deleted: %.1f%%)\n(total=%s, kept=%s)",
                          100 * kept_n / total_n, 100 * (1 - kept_n / total_n), comma(total_n), comma(kept_n))),
      hjust = 1.02, vjust = 1.2, size = 3.5, inherit.aes = FALSE
    )

  if (!is.null(length_xlim)) p <- p + coord_cartesian(xlim = length_xlim)
  p
}

plot_cov <- function() {
  dt <- rbindlist(list(
    read_depth(file.path(qc_dir, "reads_coverage.before.tsv"), "Before"),
    read_depth(file.path(qc_dir, "reads_coverage.after.tsv"),  "After")
  ))

  dt[, label := factor(label, levels = c("Before", "After"))]

  ymax <- coverage_y_max %||% max(dt$depth, na.rm = TRUE)

  p <- ggplot(dt, aes(x = pos, ymin = 0, ymax = depth, color = label, linewidth = label)) +
    geom_linerange(key_glyph = draw_key_path) +
    scale_color_manual(values = COL_COV) +
    scale_linewidth_manual(values = c(Before = 0.35, After = 0.55)) +
    scale_y_continuous(
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.03)),
      labels = comma
    ) +
    labs(
      title = paste0(sample_label, " — Depth before vs after duplication removal"),
      x = "Position",
      y = "Depth",
      color = NULL
    ) +
    guides(
      linewidth = "none",
      color = guide_legend(
        override.aes = list(
          linewidth = c(0.8, 0.9)
        )
      )
    ) +
    theme_qc(11) +
    theme(legend.position = "bottom")

  if (!is.null(coverage_xlim)) p <- p + coord_cartesian(xlim = coverage_xlim)
  p
}

st <- read_stats()
summary_plot <- wrap_plots(
  list(
    (plot_qc_ab(st) | plot_qc_cat(st)) + plot_layout(widths = c(1.05, 1.0)),
    plot_fraglen(),
    plot_cov()
  ),
  ncol = 1
)

save_plot(summary_plot, file.path(out_dir, paste0(out_prefix, "_summary")))
cat("Done. Output directory:\n", out_dir, "\n", sep = "")
