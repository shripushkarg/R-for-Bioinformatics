library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
BiocManager::install("fgsea")
library("fgsea")


#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts <- as.matrix(read.table(counts_csv, header=TRUE, row.names=1))
  coldata <- read.csv(metafile_csv, row.names=1)
  coldata <- coldata[coldata$timepoint %in% selected_times, ]
  
  coldata$timepoint <- relevel(factor(coldata$timepoint), ref = "vP0")
  ref_samplename <- c(coldata$samplename)
  counts <- as.data.frame(counts) %>% dplyr::select(all_of(ref_samplename))
  counts <- as.matrix(counts)
  
  se <- SummarizedExperiment(assays=list(counts=counts), colData=coldata)
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  results <- as.data.frame(results(dds))
  result_list <- list(dds = dds, results = results)
  return(result_list)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
    row_names <- rownames(deseq2_res)
    labeled_results <- as_tibble(deseq2_res)
    labeled_results <- labeled_results %>% mutate(genes = row_names) %>%
      mutate(volc_plot_status = case_when(
        padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
    ))
    # Remove rows with NAs
    labeled_results <- na.omit(labeled_results)
  
  return(labeled_results)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_plot <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.02, fill = "lightblue", color = "black") +
    labs(x = "pvalue", y = "count") +
    theme_minimal()
  
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  significant_genes <- labeled_results %>% filter(padj <= padj_threshold)
  log2fc_plot <- ggplot(significant_genes, aes(x = log2FoldChange)) +
    geom_histogram(binwidth = 0.2, fill = "lightblue", color = "black") +
    labs(x = "log2FoldChange", y = "count") +
    theme_minimal()
  
  return(log2fc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- head(arrange(labeled_results, padj), num_genes)
  top_gene_ids <- rownames(top_genes)
  valid_gene_ids <- intersect(top_gene_ids, rownames(counts(dds_obj)))
  
  normalized_counts <- counts(dds_obj)[valid_gene_ids, ]
  normalized_counts <- t(normalized_counts)

  data_for_plot <- as.data.frame(normalized_counts)
  data_for_plot$Sample <- rownames(normalized_counts)
  
  p <- ggplot(data_for_plot, aes(x = Sample, y = value, color = Sample)) +
    geom_point() +
    labs(title = "Normalized Counts for Top Genes",
         x = "Sample",
         y = "Normalized Counts") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  volcano_plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point() +
    labs(x = "log2FoldChange", y = "-log10(padj)", color = "volc_plot_status") +
    theme_minimal()
  
  return(volcano_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read.table(id2gene_path, header = FALSE, sep = "\t")
  colnames(id2gene) <- c("genes", "ID")
  
  merged_results <- labeled_results %>%
    left_join(id2gene, by = c("genes" = "genes")) %>%
    arrange(desc(log2FoldChange))
  
  ranked_log2fc <- merged_results$log2FoldChange
  names(ranked_log2fc) <- merged_results$ID
  ranked_log2fc <- na.omit(ranked_log2fc)
  
  return(ranked_log2fc)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  gene_sets_of_interest <- fgsea::gmtPathways(gmt_file_path)
  fgsea_results <- fgsea(gene_sets_of_interest, rnk_list, minSize = min_size, maxSize = max_size)
  fgsea_results <- as_tibble(fgsea_results)
  return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top_positive_pathways <- fgsea_results %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(num_paths)
  
  top_negative_pathways <- fgsea_results %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(num_paths)
  
  top_pathways_combined <- rbind(
    top_positive_pathways, top_negative_pathways
  )
  
  plot <- ggplot(top_pathways_combined, aes(x = reorder(leadingEdge, NES), y = NES, fill = pos_neg)) +
    geom_bar(stat = "identity") +
    labs(title = "fgsea results for Hallmark MSigDB gene set", x = "Normalized Enrichment Score (NES)", y = "Pathway Name") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) # Customize fill colors
  
  return(plot)
}

