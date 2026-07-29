#!/usr/bin/env Rscript
# ============================================================================
# analysis_wilcoxon_fdr.R
#
# Analisis estadistico comparativo para AMRFinderPlus y AMRProfiler, sobre
# Sensitivity y F1-score, a partir de all_metrics_FIXED.tsv.
#
# Metodologia:
# - Cada 'sim' usa el MISMO seed de InSilicoSeq entre datasets y entre
#   herramientas -> diseno PAREADO (Wilcoxon signed-rank, no Mann-Whitney).
# - (A) Entre datasets, por herramienta: wilcox.test(paired=TRUE) pairwise
#   entre los 4 datasets (6 comparaciones), para cada herramienta por
#   separado.
# - (B) Entre herramientas, por dataset: wilcox.test(paired=TRUE),
#   AMRFinderPlus vs AMRProfiler dentro de cada dataset.
# - Correccion Benjamini-Hochberg FDR (p.adjust method="BH") aplicada
#   dentro de cada familia de comparaciones (tool x metric para A;
#   dataset x metric para B), no globalmente.
#
# NOTA: 'control' se excluye porque Sensitivity/F1 son siempre 0 por diseno.
# NOTA: Specificity se excluyo del analisis -- con el diseno actual del
#   experimento (cada sim usa su propio genoma/reads) la clase negativa no
#   puede generar senal real de deteccion, y la metrica no discrimina
#   (siempre ~1.000, ver eval_benchmark.py / compute_specificity.py).
#
# Uso:
#   Rscript analysis_wilcoxon_fdr.R all_metrics_FIXED.tsv stats_output_R/
# ============================================================================

suppressMessages({
  library(ggplot2)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)
metrics_file <- if (length(args) >= 1) args[1] else "all_metrics_FIXED.tsv"
outdir       <- if (length(args) >= 2) args[2] else "stats_output_R"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

METRICS <- c("Sensitivity", "F1")
TOOLS   <- c("AMRFinderPlus", "AMRProfiler")

# -- Carga y limpieza --------------------------------------------------------
df <- read.delim(metrics_file, stringsAsFactors = FALSE)
df <- df[df$dataset != "control", ]
df$sim <- as.character(df$sim)
df$dataset <- as.character(df$dataset)
df$tool <- as.character(df$tool)

datasets <- sort(unique(df$dataset))
cat(sprintf("[info] %d filas cargadas, datasets: %s\n", nrow(df), paste(datasets, collapse = ", ")))

# -- Wilcoxon pareado con manejo de casos degenerados ------------------------
paired_wilcoxon <- function(x, y) {
  d <- x - y
  if (all(d == 0)) return(list(stat = NA, p = 1.0))
  res <- tryCatch(
    suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(stat = NA, p = 1.0))
  list(stat = unname(res$statistic), p = res$p.value)
}

# -- (A) Entre datasets, por herramienta -------------------------------------
between_datasets <- function(df) {
  rows <- list()
  for (tool in TOOLS) {
    for (metric in METRICS) {
      sub <- df[df$tool == tool, c("dataset", "sim", metric)]
      pairs <- combn(datasets, 2, simplify = FALSE)
      pvals <- c()
      keys  <- list()
      for (pr in pairs) {
        d1 <- pr[1]; d2 <- pr[2]
        s1 <- sub[sub$dataset == d1, ]
        s2 <- sub[sub$dataset == d2, ]
        m <- merge(s1, s2, by = "sim", suffixes = c("_1", "_2"))
        m <- m[complete.cases(m), ]
        if (nrow(m) < 2) next
        x <- m[[paste0(metric, "_1")]]
        y <- m[[paste0(metric, "_2")]]
        res <- paired_wilcoxon(x, y)
        pvals <- c(pvals, res$p)
        keys[[length(keys) + 1]] <- list(tool = tool, metric = metric,
                                          group1 = d1, group2 = d2,
                                          statistic = res$stat, n_pairs = nrow(m))
      }
      if (length(pvals) == 0) next
      p_adj <- p.adjust(pvals, method = "BH")
      for (i in seq_along(keys)) {
        k <- keys[[i]]
        rows[[length(rows) + 1]] <- data.frame(
          comparison = "between_datasets", tool = k$tool, metric = k$metric,
          group1 = k$group1, group2 = k$group2, statistic = k$statistic,
          pvalue = pvals[i], pvalue_adj = p_adj[i],
          significant = p_adj[i] < 0.05, n_pairs = k$n_pairs
        )
      }
    }
  }
  do.call(rbind, rows)
}

# -- (B) Entre herramientas, por dataset -------------------------------------
between_tools <- function(df) {
  rows <- list()
  for (ds in datasets) {
    sub <- df[df$dataset == ds, ]
    pvals <- c()
    keys  <- list()
    for (metric in METRICS) {
      s1 <- sub[sub$tool == "AMRFinderPlus", c("sim", metric)]
      s2 <- sub[sub$tool == "AMRProfiler", c("sim", metric)]
      m <- merge(s1, s2, by = "sim", suffixes = c("_f", "_p"))
      m <- m[complete.cases(m), ]
      if (nrow(m) < 2) next
      x <- m[[paste0(metric, "_f")]]
      y <- m[[paste0(metric, "_p")]]
      res <- paired_wilcoxon(x, y)
      pvals <- c(pvals, res$p)
      keys[[length(keys) + 1]] <- list(dataset = ds, metric = metric,
                                        statistic = res$stat, n_pairs = nrow(m))
    }
    if (length(pvals) == 0) next
    p_adj <- p.adjust(pvals, method = "BH")
    for (i in seq_along(keys)) {
      k <- keys[[i]]
      rows[[length(rows) + 1]] <- data.frame(
        comparison = "between_tools", dataset = k$dataset, metric = k$metric,
        group1 = "AMRFinderPlus", group2 = "AMRProfiler", statistic = k$statistic,
        pvalue = pvals[i], pvalue_adj = p_adj[i],
        significant = p_adj[i] < 0.05, n_pairs = k$n_pairs
      )
    }
  }
  do.call(rbind, rows)
}

stats_a <- between_datasets(df)
stats_b <- between_tools(df)

write.table(stats_a, file.path(outdir, "wilcoxon_between_datasets.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(stats_b, file.path(outdir, "wilcoxon_between_tools.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_sig <- sum(stats_a$significant) + sum(stats_b$significant)
n_tot <- nrow(stats_a) + nrow(stats_b)
cat(sprintf("[info] comparaciones significativas (p-adj < 0.05): %d / %d\n", n_sig, n_tot))

sig_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

cat("\n=== (A) Entre datasets, comparaciones significativas ===\n")
sig_a <- stats_a[stats_a$significant, ]
if (nrow(sig_a) == 0) {
  cat("  (ninguna)\n")
} else {
  for (i in seq_len(nrow(sig_a))) {
    r <- sig_a[i, ]
    cat(sprintf("  %-15s %-12s %s vs %s: p-adj=%.4f %s (n=%d)\n",
                r$tool, r$metric, r$group1, r$group2, r$pvalue_adj,
                sig_stars(r$pvalue_adj), r$n_pairs))
  }
}

cat("\n=== (B) AMRFinderPlus vs AMRProfiler, comparaciones significativas ===\n")
sig_b <- stats_b[stats_b$significant, ]
if (nrow(sig_b) == 0) {
  cat("  (ninguna)\n")
} else {
  for (i in seq_len(nrow(sig_b))) {
    r <- sig_b[i, ]
    cat(sprintf("  %-8s %-12s p-adj=%.4f %s (n=%d)\n",
                r$dataset, r$metric, r$pvalue_adj, sig_stars(r$pvalue_adj), r$n_pairs))
  }
}

# ============================================================================
# FIGURAS
# ============================================================================

make_barplots <- function(df, stats_a, outdir) {
  for (tool in TOOLS) {
    plots <- list()
    for (metric in METRICS) {
      sub <- df[df$tool == tool, ]
      sub$dataset <- factor(sub$dataset, levels = datasets)

      p <- ggplot(sub, aes(x = dataset, y = .data[[metric]])) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA, fill = "#4C72B0") +
        geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, color = "black") +
        ylim(-0.05, 1.3) +
        labs(title = paste(tool, "-", metric), x = "Dataset", y = metric) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(face = "bold"))

      pair_stats <- stats_a[stats_a$tool == tool & stats_a$metric == metric & stats_a$significant, ]
      if (nrow(pair_stats) > 0) {
        y_max <- max(sub[[metric]], na.rm = TRUE)
        step <- 0.09
        for (i in seq_len(nrow(pair_stats))) {
          x1 <- match(pair_stats$group1[i], datasets)
          x2 <- match(pair_stats$group2[i], datasets)
          y <- min(1.25, y_max + step * i)
          p <- p +
            annotate("segment", x = x1, xend = x2, y = y, yend = y) +
            annotate("segment", x = x1, xend = x1, y = y, yend = y - 0.015) +
            annotate("segment", x = x2, xend = x2, y = y, yend = y - 0.015) +
            annotate("text", x = (x1 + x2) / 2, y = y + 0.02,
                     label = sig_stars(pair_stats$pvalue_adj[i]), size = 5)
        }
      }
      plots[[metric]] <- p
    }

    out_f <- file.path(outdir, paste0("barplot_", gsub(" ", "_", tool), ".png"))
    png(out_f, width = 500 * length(METRICS), height = 450, res = 100)
    if (length(plots) > 1) {
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        do.call(gridExtra::grid.arrange, c(plots, ncol = length(plots)))
      } else {
        print(plots[[1]])
      }
    } else {
      print(plots[[1]])
    }
    dev.off()
    cat(sprintf("[OK] %s\n", out_f))
  }
}

make_heatmap_pvalues <- function(stats_a, outdir) {
  for (tool in TOOLS) {
    for (metric in METRICS) {
      sub <- stats_a[stats_a$tool == tool & stats_a$metric == metric, ]
      if (nrow(sub) == 0) next

      mat <- matrix(NA, nrow = length(datasets), ncol = length(datasets),
                     dimnames = list(datasets, datasets))
      for (i in seq_len(nrow(sub))) {
        mat[sub$group1[i], sub$group2[i]] <- sub$pvalue_adj[i]
        mat[sub$group2[i], sub$group1[i]] <- sub$pvalue_adj[i]
      }
      melted <- melt(mat, varnames = c("d1", "d2"), value.name = "p_adj")
      melted <- melted[!is.na(melted$p_adj), ]

      p <- ggplot(melted, aes(x = d1, y = d2, fill = p_adj)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.3f", p_adj)), size = 4) +
        scale_fill_gradient(low = "#440154", high = "#fde725", limits = c(0, 1),
                             name = "p-adj\n(BH-FDR)") +
        labs(title = paste(tool, "-", metric, "\n(p-adj, pareado por sim)"),
             x = NULL, y = NULL) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      out_f <- file.path(outdir, sprintf("heatmap_pvalues_%s_%s.png",
                                          gsub(" ", "_", tool), metric))
      ggsave(out_f, p, width = 5.5, height = 5, dpi = 150)
      cat(sprintf("[OK] %s\n", out_f))
    }
  }
}

make_heatmap_tools <- function(stats_b, outdir) {
  if (nrow(stats_b) == 0) return(invisible())
  mat <- dcast(stats_b, dataset ~ metric, value.var = "pvalue_adj")
  melted <- melt(mat, id.vars = "dataset", variable.name = "metric", value.name = "p_adj")

  p <- ggplot(melted, aes(x = metric, y = dataset, fill = p_adj)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", p_adj)), size = 4) +
    scale_fill_gradient(low = "#440154", high = "#fde725", limits = c(0, 1),
                         name = "p-adj\n(BH-FDR)") +
    labs(title = "AMRFinderPlus vs AMRProfiler por dataset\n(p-adj, pareado por sim)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  out_f <- file.path(outdir, "heatmap_pvalues_tool_vs_tool.png")
  ggsave(out_f, p, width = 5.5, height = 5, dpi = 150)
  cat(sprintf("[OK] %s\n", out_f))
}

make_heatmap_summary <- function(df, outdir) {
  agg <- aggregate(df[METRICS], by = list(dataset = df$dataset, tool = df$tool), FUN = mean)

  for (tool in TOOLS) {
    sub <- agg[agg$tool == tool, c("dataset", METRICS)]
    melted <- melt(sub, id.vars = "dataset", variable.name = "metric", value.name = "value")

    p <- ggplot(melted, aes(x = metric, y = dataset, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.3f", value)), size = 4) +
      scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen",
                            midpoint = 0.5, limits = c(0, 1), name = "valor\npromedio") +
      labs(title = paste(tool, "- rendimiento promedio por dataset"), x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    out_f <- file.path(outdir, sprintf("heatmap_summary_%s.png", gsub(" ", "_", tool)))
    ggsave(out_f, p, width = 5, height = 5, dpi = 150)
    cat(sprintf("[OK] %s\n", out_f))
  }
}

make_barplots(df, stats_a, outdir)
make_heatmap_pvalues(stats_a, outdir)
make_heatmap_tools(stats_b, outdir)
make_heatmap_summary(df, outdir)

cat("\n[OK] listo.\n")
