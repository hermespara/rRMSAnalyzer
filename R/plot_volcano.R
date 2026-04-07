#' Plot a volcano plot of differential sites.
#'
#' Display a volcano plot of the differential sites finding, with the Delta C-score (or log2 fold change)
#' on the x-axis and the -log10 of the p-value on the y-axis.
#'
#' @details
#' To be considered as differential, a site must follow two conditions :
#'
#'   - Have a significant p-value on the selected statistical
#'    test on c-score between conditions.
#'   - Have a absolute Delta C-score between a condition and the reference
#'   above a certain cutoff.
#'
#' Both the p-value cutoff and c-score range cutoff can be changed with
#' `p_cutoff` and `cscore_cutoff` parameters respectively.
#' @md
#' @param ribo A SummarizedExperiment object.
#' @param factor_column Metadata column used to group samples by.
#' @param reference_level The level in `factor_column` to use as the reference for calculating the Delta C-score. If NULL, the first level is used.
#' @param target_condition The target level in `factor_column` to compare against the reference. If NULL, all levels other than the reference are used.
#' @param p_cutoff Cutoff for the adjusted p-value (or raw p-value depending on `y_axis_p`) to draw a horizontal threshold line.
#' @param cscore_cutoff Cutoff for the absolute Delta C-score (mean difference) to draw vertical threshold lines.
#' @param statistical_test Statistical test used to compute p-values.
#' One of "kruskal" (wilcox is automatically used if there are 2 groups) or "t.test".
#' @param adjust_pvalues_method Method used to adjust p-value (one of p.adjust.methods).
#' @param y_axis_p Which p-value to use for the y-axis. Either "p.adj" or "p.val".
#' @param object_only Return the results of the statistical test and C-score mean range in a dataframe directly, without plotting.
#' @return a ggplot object or a dataframe if `object_only` is TRUE.
#' @export
#' @import ggplot2
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' # ribo_toy <- rename_rna(ribo_toy)
#' # ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' # plot_volcano(ribo_toy, "condition", reference_level = "RNA ref", p_cutoff=0.1)
plot_volcano <- function(ribo, factor_column,
                         reference_level = NULL,
                         target_condition = NULL,
                         p_cutoff = 1e-02,
                         cscore_cutoff = 0.05,
                         statistical_test = "kruskal",
                         adjust_pvalues_method = "fdr",
                         y_axis_p = "p.adj",
                         object_only = FALSE) {
    check_is_se(ribo)
    check_type(factor_column, "character", "factor_column", length = 1)
    if (!is.null(reference_level)) check_type(reference_level, "character", "reference_level", length = 1)
    if (!is.null(target_condition)) {
        check_type(target_condition, "character", "target_condition")
    }
    check_type(p_cutoff, "numeric", "p_cutoff", length = 1)
    check_type(cscore_cutoff, "numeric", "cscore_cutoff", length = 1)
    check_type(statistical_test, "character", "statistical_test", length = 1)
    check_in_set(tolower(statistical_test), c("kruskal", "t.test"), "statistical_test")
    check_type(adjust_pvalues_method, "character", "adjust_pvalues_method", length = 1)
    check_type(y_axis_p, "character", "y_axis_p", length = 1)
    check_in_set(y_axis_p, c("p.adj", "p.val"), "y_axis_p")
    check_type(object_only, "logical", "object_only", length = 1)

    site <- p.val <- p.adj <- delta_cscore <- regulation <- label <- y_value <- target_condition <- NULL

    df_stats <- wrapper_kruskal_test(
        ribo = ribo,
        adjust_pvalues_method = adjust_pvalues_method,
        factor_column = factor_column,
        statistical_test = statistical_test
    )

    # Calculate exact condition means and actual Delta C-scores
    metadata <- as.data.frame(SummarizedExperiment::colData(ribo))
    if (!"samplename" %in% colnames(metadata)) metadata$samplename <- rownames(metadata)

    cscore_matrix <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
    cscore_matrix <- t(cscore_matrix)
    cscore_matrix <- as.data.frame(cscore_matrix[match(metadata[, "samplename"], rownames(cscore_matrix)), ])

    df_mean_each_group <- stats::aggregate(cscore_matrix, list(condition = metadata[, factor_column]), mean, na.rm = TRUE)

    site_means <- t(df_mean_each_group[, -1, drop = FALSE])
    colnames(site_means) <- df_mean_each_group$condition
    site_means <- as.data.frame(site_means)
    site_means$site <- rownames(site_means)

    if (is.null(reference_level)) {
        lvl <- df_mean_each_group$condition
        if ("control" %in% tolower(lvl)) {
            reference_level <- lvl[tolower(lvl) == "control"][1]
        } else if ("ctrl" %in% tolower(lvl)) {
            reference_level <- lvl[tolower(lvl) == "ctrl"][1]
        } else if ("wt" %in% tolower(lvl)) {
            reference_level <- lvl[tolower(lvl) == "wt"][1]
        } else if ("ref" %in% tolower(lvl) || "rna ref" %in% tolower(lvl)) {
            reference_level <- lvl[tolower(lvl) == "ref" | tolower(lvl) == "rna ref"][1]
        } else {
            reference_level <- lvl[1]
        }
        message("No reference_level provided. Automatically using '", reference_level, "' as inference baseline.")
    } else if (!reference_level %in% df_mean_each_group$condition) {
        cli::cli_abort("reference_level '{reference_level}' not found in factor_column")
    }

    target_levels <- setdiff(colnames(site_means)[-ncol(site_means)], reference_level)
    if (!is.null(target_condition)) {
        # Check if provided targets exist
        invalid_targets <- setdiff(target_condition, df_mean_each_group$condition)
        if (length(invalid_targets) > 0) {
            cli::cli_abort("target_condition(s) '{invalid_targets}' not found in factor_column")
        }
        target_levels <- intersect(target_levels, target_condition)
    }
    if (length(target_levels) == 0) {
        cli::cli_abort("Found {length(target_levels)} target conditions valid against the reference.")
    }

    delta_list <- list()
    for (t_lvl in target_levels) {
        temp_df <- data.frame(
            site = site_means$site,
            target_condition = t_lvl,
            delta_cscore = site_means[[t_lvl]] - site_means[[reference_level]],
            mean_target = site_means[[t_lvl]],
            mean_reference = site_means[[reference_level]]
        )
        delta_list[[t_lvl]] <- temp_df
    }
    delta_df <- do.call(rbind, delta_list)
    rownames(delta_df) <- NULL

    df_stats <- merge(df_stats, delta_df, by = "site")

    if (object_only) {
        return(df_stats)
    }

    if (nrow(df_stats) == 0) {
        return(ggplot() +
            annotate("text",
                x = 0, y = 0, size = 8,
                label = "No differential site found !"
            ) +
            theme_void())
    }

    df_stats$regulation <- ifelse(
        df_stats[[y_axis_p]] < p_cutoff & df_stats$delta_cscore < -cscore_cutoff,
        "Downregulated",
        ifelse(
            df_stats[[y_axis_p]] < p_cutoff & df_stats$delta_cscore > cscore_cutoff,
            "Upregulated",
            "Not Significant"
        )
    )
    df_stats$regulation <- factor(
        df_stats$regulation,
        levels = c("Downregulated", "Upregulated", "Not Significant")
    )

    df_stats$y_value <- -log10(df_stats[[y_axis_p]])

    most_signi <- df_stats$site[which(df_stats[[y_axis_p]] < p_cutoff & abs(df_stats$delta_cscore) > cscore_cutoff)]

    df_stats$label <- ""
    df_stats$label[df_stats$site %in% most_signi] <- df_stats$site[df_stats$site %in% most_signi]

    p_volcano <- ggplot(df_stats, aes(x = delta_cscore, y = y_value, color = regulation, label = label))

    if (length(target_levels) > 1) {
        p_volcano <- p_volcano + geom_point(aes(shape = target_condition), alpha = 0.8) +
            labs(shape = "Condition (Target vs Ref)")
    } else {
        p_volcano <- p_volcano + geom_point(alpha = 0.8)
    }

    p_volcano <- p_volcano +
        geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        geom_vline(xintercept = c(-cscore_cutoff, cscore_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
        theme_rRMSAnalyzer() +
        scale_color_manual(values = c(
            "Downregulated" = "darkblue",
            "Upregulated" = "darkred",
            "Not Significant" = "grey"
        )) +
        labs(
            x = paste0("Delta C-score (Target - ", reference_level, ")"),
            y = paste0("-Log10(", y_axis_p, ")"),
            color = "Significance",
            caption = paste(length(unique(df_stats$site)), "sites")
        ) +
        theme(legend.position = "top")

    if (requireNamespace("ggrepel", quietly = TRUE) && any(df_stats$label != "")) {
        p_volcano <- p_volcano + ggrepel::geom_text_repel(
            size = 3,
            max.overlaps = 15,
            show.legend = FALSE
        )
    }

    return(p_volcano)
}
