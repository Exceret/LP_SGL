#' @title Perform Differential Expression Analysis using limma
#' @description
#' Conducts differential expression analysis on bulk RNA-seq data using the limma package.
#' This function identifies genes that are differentially expressed between two phenotypic groups.
#'
#' @param bulk_matrix Bulk expression matrix with genes in rows and samples in columns.
#' Should be normalized and log-transformed (e.g., logCPM).
#' @param phenotype Vector of phenotype labels for each sample. Must contain exactly two unique groups.
#' @param logFC_threshold Minimum absolute log fold change for significance (default: 1)
#' @param pval_threshold Maximum p-value for significance (default: 0.05)
#' @param adjust_method Method for p-value adjustment. Options include "BH", "bonferroni",
#' "holm", etc. (default: "BH" for Benjamini-Hochberg)
#' @param ... Additional argumentss, like verbose (default: TRUE)
#'
#' @return A list containing:
#' \item{all_results}{Complete results from limma's topTable with statistics for all genes}
#' \item{significant_DEGs}{Filtered differentially expressed genes meeting significance thresholds}
#' \item{fit}{The fitted limma model object for further analysis}
#' \item{n_significant}{Total number of significant DEGs}
#' \item{n_upregulated}{Number of significantly upregulated genes}
#' \item{n_downregulated}{Number of significantly downregulated genes}
#' \item{comparison}{Description of the comparison performed}
#'
#' @details
#' This function performs the following steps:
#' 1. Creates a design matrix based on the phenotype groups
#' 2. Fits a linear model using limma's empirical Bayes methods
#' 3. Computes contrasts between the two groups
#' 4. Extracts results and filters based on significance thresholds
#'
#' The method is particularly suited for bulk RNA-seq data and provides robust
#' variance estimation through empirical Bayes moderation.
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' bulk_matrix <- matrix(rnorm(1000*20), nrow=1000, ncol=20)
#' rownames(bulk_matrix) <- paste0("Gene", 1:1000)
#' colnames(bulk_matrix) <- paste0("Sample", 1:20)
#' phenotype <- rep(c("Control", "Treatment"), each=10)
#'
#' # Perform DEG analysis
#' deg_results <- perform_DEG_analysis(
#' bulk_matrix = bulk_matrix,
#' phenotype = phenotype,
#' logFC_threshold = 1,
#' pval_threshold = 0.05
#' )
#'
#' # View significant DEGs
#' head(deg_results$significant_DEGs)
#'
#' # Summary statistics
#' cat("Total significant DEGs:", deg_results$n_significant, "\n")
#' cat("Upregulated:", deg_results$n_upregulated, "\n")
#' cat("Downregulated:", deg_results$n_downregulated, "\n")
#' }
#'
#' @seealso
#' [limma::lmFit()], [limma::eBayes()], [limma::topTable()]
#' @export
perform_DEG_analysis <- function(
    bulk_matrix,
    phenotype,
    logFC_threshold = 1,
    pval_threshold = 0.05,
    adjust_method = "BH",
    ...
) {
    dots <- rlang::list2(...)
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")
    # Input validation
    if (!is.vector(phenotype)) {
        cli::cli_abort(c(
            "x" = "phenotype must be a vector, currently only binomial phenotype is supported"
        ))
    }
    if (length(table(phenotype)) != 2) {
        cli::cli_abort(c(
            "x" = "phenotype must contain exactly two unique groups"
        ))
    }
    adjust_method <- SigBridgeRUtils::MatchArg(
        adjust_method,
        c(
            "BH",
            "BY",
            "bonferroni",
            "holm",
            "hochberg",
            "hommel",
            "fdr",
            "none"
        )
    )

    if (verbose) {
        ts_cli$cli_alert_info(cli::col_green(
            "Performing DEG analysis on bulk RNA-seq data"
        ))
    }
    # Load limma functions
    makeContrasts <- getExportedValue("limma", "makeContrasts")
    lmFit <- getExportedValue("limma", "lmFit")
    contrasts.fit <- getExportedValue("limma", "contrasts.fit")
    eBayes <- getExportedValue("limma", "eBayes")
    topTable <- getExportedValue("limma", "topTable")

    # Create groups
    unique_phenotypes <- unique(phenotype)
    phenotype1_idx <- which(phenotype == unique_phenotypes[1])
    phenotype2_idx <- which(phenotype == unique_phenotypes[2])

    group <- rep(
        c(unique_phenotypes[1], unique_phenotypes[2]),
        c(length(phenotype1_idx), length(phenotype2_idx))
    )

    # Design matrix
    design <- stats::model.matrix(~ 0 + factor(group))
    colnames(design) <- levels(factor(group))
    rownames(design) <- colnames(bulk_matrix)

    # Contrast matrix
    contrast.matrix <- makeContrasts(
        paste0(unique_phenotypes[1], "-", unique_phenotypes[2]),
        levels = design
    )

    # Fit linear model
    fit <- lmFit(bulk_matrix, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    # Extract results
    DEG_all <- topTable(fit2, coef = 1, n = Inf, adjust.method = adjust_method)

    # Filter significant DEGs
    DEGs <- subset(
        DEG_all,
        abs(DEG_all$logFC) > logFC_threshold & DEG_all$P.Value < pval_threshold
    )

    # Return results
    return(list(
        all_results = DEG_all,
        significant_DEGs = DEGs,
        fit = fit2,
        n_significant = nrow(DEGs),
        n_upregulated = sum(DEGs$logFC > 0),
        n_downregulated = sum(DEGs$logFC < 0)
    ))
}
