#!/usr/bin/env Rscript
# plot_minorfreq.R - Visualize minor allele frequency distributions

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(RColorBrewer)

# Function to read and process minorfreq data
read_minorfreq <- function(file_path, sample_name = NULL) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Extract sample name from filename if not provided
  if (is.null(sample_name)) {
    sample_name <- gsub("\\.minorfreq\\.txt$", "", basename(file_path))
  }
  
  # Read the data
  data <- read_table(file_path, col_names = c("count", "freq", "pair"), 
                     col_types = "idc", show_col_types = FALSE)
  
  # Add sample name and calculate proportions
  data$sample <- sample_name
  data$total_sites <- sum(data$count)
  data$proportion <- data$count / data$total_sites
  
  # Add transition/transversion classification
  data$type <- ifelse(data$pair %in% c("AG", "CT"), "Transition", "Transversion")
  
  return(data)
}

# Function to create frequency distribution plot
plot_freq_distribution <- function(data, title_prefix = "") {
  p1 <- ggplot(data, aes(x = freq, y = count, fill = pair)) +
    geom_col(position = "stack", alpha = 0.8) +
    scale_fill_brewer(type = "qual", palette = "Set3", name = "Allele Pair") +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = paste0(title_prefix, "Minor Allele Frequency Distribution"),
      subtitle = "Count of sites by frequency and allele pair (log scale)",
      x = "Minor Allele Frequency",
      y = "Number of Sites (log scale)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p1)
}

# Function to create transition/transversion plot
plot_ti_tv <- function(data, title_prefix = "") {
  # Calculate Ti/Tv ratio
  ti_tv_summary <- data %>%
    group_by(sample, type) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = total_count, values_fill = 0) %>%
    mutate(
      Ti_Tv_ratio = Transition / Transversion,
      total = Transition + Transversion
    )
  
  p2 <- ggplot(data, aes(x = freq, y = count, fill = type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Transition" = "#E31A1C", "Transversion" = "#1F78B4"),
                      name = "Mutation Type") +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = paste0(title_prefix, "Transitions vs Transversions"),
      subtitle = paste0("Ti/Tv ratio: ", round(ti_tv_summary$Ti_Tv_ratio, 2)),
      x = "Minor Allele Frequency", 
      y = "Number of Sites (log scale)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p2)
}

# Function to create allele pair comparison
plot_pair_comparison <- function(data, title_prefix = "") {
  pair_summary <- data %>%
    group_by(sample, pair) %>%
    summarise(
      total_sites = sum(count),
      mean_freq = weighted.mean(freq, count),
      .groups = "drop"
    )
  
  p3 <- ggplot(pair_summary, aes(x = reorder(pair, total_sites), y = total_sites, fill = pair)) +
    geom_col(alpha = 0.8) +
    scale_fill_brewer(type = "qual", palette = "Set3", guide = "none") +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = paste0(title_prefix, "Heterozygous Sites by Allele Pair"),
      subtitle = "Total count of heterozygous sites for each nucleotide pair (log scale)",
      x = "Allele Pair",
      y = "Total Sites (log scale)"
    ) +
    theme_minimal() +
    coord_flip()
  
  return(p3)
}

# Function to create frequency spectrum comparison
plot_freq_spectrum <- function(data, title_prefix = "") {
  p4 <- ggplot(data, aes(x = freq, y = proportion, color = pair)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_brewer(type = "qual", palette = "Set3", name = "Allele Pair") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0(title_prefix, "Minor Allele Frequency Spectrum"),
      subtitle = "Proportion of total heterozygous sites",
      x = "Minor Allele Frequency",
      y = "Proportion of Sites (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p4)
}

# Main function to create all plots for a sample
create_sample_plots <- function(minorfreq_file, output_prefix, fmt = "pdf") {
  # Read data
  cat("Processing:", minorfreq_file, "\n")
  data <- read_minorfreq(minorfreq_file)
  
  if (is.null(data) || nrow(data) == 0) {
    warning("No data to plot for", minorfreq_file)
    return()
  }
  
  sample_name <- unique(data$sample)
  title_prefix <- paste0(sample_name, " - ")
  
  # Create individual plots
  p1 <- plot_freq_distribution(data, title_prefix)
  p2 <- plot_ti_tv(data, title_prefix)
  p3 <- plot_pair_comparison(data, title_prefix)
  p4 <- plot_freq_spectrum(data, title_prefix)
  
  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4)
  combined_plot <- combined_plot + 
    plot_annotation(
      title = paste("Minor Allele Frequency Analysis:", sample_name),
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
  
  # Save plots in requested format (pdf/svg/png)
  fmt <- tolower(fmt)
  if (fmt == "pdf") {
    output_file <- paste0(output_prefix, "_minorfreq_plots.pdf")
    ggsave(output_file, combined_plot, width = 16, height = 12, device = pdf)
  } else if (fmt == "svg") {
    output_file <- paste0(output_prefix, "_minorfreq_plots.svg")
    ggsave(output_file, combined_plot, width = 16, height = 12, device = "svg")
  } else {
    output_file <- paste0(output_prefix, "_minorfreq_plots.png")
    ggsave(output_file, combined_plot, width = 16, height = 12, dpi = 300)
  }
  cat("Saved plot:", output_file, "\n")
  
  # Save summary statistics
  summary_stats <- data %>%
    group_by(pair, type) %>%
    summarise(
      total_sites = sum(count),
      mean_freq = weighted.mean(freq, count),
      median_freq = median(rep(freq, count)),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_name)
  
  summary_file <- paste0(output_prefix, "_minorfreq_summary.txt")
  write_tsv(summary_stats, summary_file)
  cat("Saved summary:", summary_file, "\n")
  
  return(list(data = data, plots = combined_plot, summary = summary_stats))
}

# Function to compare multiple samples
compare_samples <- function(minorfreq_files, sample_names = NULL, output_file = "multi_sample_comparison.png") {
  if (length(minorfreq_files) < 2) {
    warning("Need at least 2 samples for comparison")
    return()
  }
  
  # Read all data
  all_data <- map2_dfr(minorfreq_files, sample_names %||% basename(minorfreq_files), 
                       ~read_minorfreq(.x, .y))
  
  if (nrow(all_data) == 0) {
    warning("No data found for comparison")
    return()
  }
  
  # Multi-sample plots
  p1 <- ggplot(all_data, aes(x = freq, y = count, fill = sample)) +
    geom_col(position = "dodge", alpha = 0.7) +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Sample") +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~pair, scales = "free_y") +
    labs(
      title = "Minor Allele Frequency Distribution Comparison",
      subtitle = "Count of sites by frequency across samples (log scale)",
      x = "Minor Allele Frequency",
      y = "Number of Sites (log scale)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- all_data %>%
    group_by(sample, type) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = total_count, values_fill = 0) %>%
    mutate(Ti_Tv_ratio = Transition / Transversion) %>%
    ggplot(aes(x = sample, y = Ti_Tv_ratio, fill = sample)) +
    geom_col(alpha = 0.8) +
    scale_fill_brewer(type = "qual", palette = "Set1", guide = "none") +
    labs(
      title = "Transition/Transversion Ratio Comparison",
      x = "Sample",
      y = "Ti/Tv Ratio"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine and save
  comparison_plot <- p1 / p2
  comparison_plot <- comparison_plot + 
    plot_annotation(
      title = "Multi-Sample Minor Allele Frequency Comparison",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
  
  ggsave(output_file, comparison_plot, width = 16, height = 12, dpi = 300)
  cat("Saved comparison plot:", output_file, "\n")
  
  return(comparison_plot)
}

# Command line interface
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    cat("Usage: Rscript plot_minorfreq.R <minorfreq_file> [--output output_prefix] [--format pdf|svg|png]\n")
    cat("  or: Rscript plot_minorfreq.R --compare file1 file2 ... --output comparison_plot.pdf --format pdf\n")
    quit(status = 1)
  }

  # parse simple args
  fmt_idx <- which(args == "--format")
  fmt <- if (length(fmt_idx) > 0) args[fmt_idx + 1] else "pdf"
  output_idx <- which(args == "--output")

  compare_mode <- "--compare" %in% args

  if (compare_mode) {
    compare_idx <- which(args == "--compare")
    files <- args[(compare_idx + 1):(ifelse(length(output_idx)>0, output_idx - 1, length(args)))]
    out_file <- if (length(output_idx) > 0) args[output_idx + 1] else "multi_sample_comparison"
    compare_samples(files, sample_names = NULL, output_file = paste0(out_file, ".", tolower(fmt)))
  } else {
    input_file <- args[1]
    output_prefix <- if (length(output_idx) > 0) args[output_idx + 1] else gsub("\\.minorfreq\\.txt$", "", input_file)
    create_sample_plots(input_file, output_prefix, fmt = fmt)
  }
}

# Run if called as script
if (!interactive()) {
  main()
}
