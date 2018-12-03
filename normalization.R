fix_low <- function(min_count = 2, min_sample_per_gene = 5, highest_exprs_top_n = 20) {
	filter_low()
	imputation()
}

normalize_check_arg <- function(norm_methods, top_n, rm_gene_type, refer_gene_id) {
	if ('CPM_top' %in% norm_methods && is.null(top_n))
		stop('top_n must be specified for CPM_top method')

	if ('CPM_rm' %in% norm_methods && is.null(rm_gene_type))
		stop('rm_gene_type must be specified for CPM_rm method')

	if ('CPM_refer' %in% norm_methods && is.null(refer_gene_id))
		stop('refer_gene_id must be specified for CPM_refer method')

}

#' @title martix normalization
#' @examples
#' \donotrun{
#'     refer_gene_id = c(
#'         "ENST00000408438.1", "ENST00000385271.1", "ENST00000607334.3",
#'         "ENST00000385059.1", "ENST00000362134.1", "ENST00000385245.1",
#'         "ENST00000385045.1", "ENST00000362117.1", "ENST00000384832.1",
#'         "ENST00000579846.3"
#'     )
#'     norm_mat(
#'         '/path/to/matrix', c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'),
#'         top_n = 20L, rm_gene_type = c('miRNA', 'piRNA'), refer_gene_id = refer_gene_id
#'     )
#' }
normalize <- function(
	counts_mat_path,
	norm_methods = c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'),
	top_n = NULL, rm_gene_type = NULL, refer_gene_id = NULL,
	sample_class_path = NULL, PCA_label_by = NULL, PAC_color_by = NULL,
	refer_gene_name = refer_gene_id,
	tmp_path='.', K = 5, N = 3
    output_dir = '.',output_file = 'norm'
) {
	normalize_check_arg(norm_methods, top_n, rm_gene_type, refer_gene_id)
	mat <- read_mat(counts_mat_path) %>% filter_low(min_count, min_sample_per_gene) %>% imputation(tmp_path, impute_path, K, N)

	if ('SCnorm' %in% norm_methods)    mat_scnorm <- norm_SCnorm(mat)
	if ('TMM' %in% norm_methods)       mat_tmm <- norm_tmm(mat)
	if ('RLE' %in% norm_methods)       mat_rle <- norm_rle(mat)
	if ('CPM' %in% norm_methods)       mat_cpm <- norm_cpm(mat)
	if ('CPM_top' %in% norm_methods)   mat_cpm_top <- norm_cpm_top(mat, top_n)
	if ('CPM_rm' %in% norm_methods)    mat_cpm_rm <- norm_cpm_rm(mat, rm_gene_type)
	if ('CPM_refer' %in% norm_methods) mat_cpm_refer <- norm_cpm_refer(mat, refer_gene_id)

	mat_names <- ls(pattern = '^mat_')

	for(mat_name in mat_names) {
		readr::write_tsv(get(mat_name), paste0(output_dir, '/', output_file, stringr::str_extract(mat_name, '_\\w+', '.tsv')))

	}

}












#' @title SCnorm normalization
#'
#' @param mat integer matrix. counts
#' @param ... other arguments passed on to [SCnorm::SCnorm()]
#'
#' @examples
#' norm_SCnorm(sim_mat*10)
#'
#' @family matrix normalization
#'
#' @export
norm_SCnorm <- function(mat, ...) {
	Conditions = rep(1, ncol(mat));
	sce <- suppressMessages(SCnorm::SCnorm(mat, Conditions, ...));
	SCnorm::normcounts(sce)
}

# norm_scater ------------------

#' @title TMM/RLE normalization by scater package
#'
#' @param mat integer matrix. counts
#'
#' @family matrix normalization
#'
#' @name norm_scater
NULL


#' @rdname  norm_scater
#'
#' @details `norm_tmm()` performs TMM normalization
#'
#' @examples
#' norm_tmm(sim_mat)
#'
#' @export
norm_tmm <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
}


#' @rdname  norm_scater
#'
#' @details `norm_rle()` performs RLE normalization
#'
#' @examples
#' norm_rle(sim_mat)
#'
#' @export
norm_rle <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
}


# norm_cpm ------------------

#' @title CPM normalization by some genes
#'
#' @param mat integer matrix. counts
#' @param row integer or logical. Use which rows (genes) as normalization factor
norm_cpm_impl <- function(mat, row) {
	t(t(mat*1e6) / colSums(mat[row, , drop = F], na.rm = T))
}


#' @title CPM normalization
#' 
#' @description CPM normalization using counts sum of _certain_ genes as scaling factor
#'
#' @param mat integer matrix. counts.
#'
#' @details some functions may throw errors
#'
#' @family matrix normalization
#'
#' @name norm_cpm
NULL



#' @rdname norm_cpm
#'
#' @details `norm_cpm_total()` uses total genes
#'
#' @examples
#' norm_cpm_total(sim_mat)
#'
#' @export
norm_cpm_total <- function(mat) {
	row_all <- nrow(mat) %>% seq_len()

	norm_cpm_impl(mat, row_all)
}

#' @rdname norm_cpm
#'
#' @param top_n integer scalar. see `norm_cpm_top()` below
#'
#' @details `norm_cpm_top()` uses top 20 genes sorted by counts (assuming `top_n = 20L`)
#'
#' @examples
#' norm_cpm_top(sim_mat, 20L)
#'
#' @export
norm_cpm_top <- function(mat, top_n) {
	if (nrow(mat) < top_n)
		stop('too few feature for CPM top n normalization')

	row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
		 {.$ix[seq_len(top_n)]}

	norm_cpm_impl(mat, row_top)
}


#' @rdname norm_cpm
#'
#' @param gene_type character. see `norm_cpm_rm()` below
#'
#' @details `norm_cpm_rm()` uses non-piRNA genes (assuming `gene_type = 'piRNA'`)
#'
#' @examples
#' norm_cpm_rm(sim_mat, c('miRNA', 'piRNA'))
#'
#' @export
norm_cpm_rm <- function(mat, gene_type) {
	if (!all(gene_type %in% rexseek::rna_type$gene_type))
		stop('unknown transcript type to remove for CPM normalization')

	id_rm <- tibble::tibble(gene_type) %>%
		dplyr::left_join(rexseek::rna_type, by = 'gene_type') %>% .$transcript_id
	row_rm <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% id_rm}

	norm_cpm_impl(mat, !row_rm)
}



#' @rdname norm_cpm
#'
#' @param refer_gene_id character. Ensembl transcript id, see `norm_cpm_refer()` below
#'
#' @details `norm_cpm_refer()` uses given reference genes
#'
#' @examples
#' norm_cpm_refer(sim_mat, suggest_refer$id)
#'
#' @export

# mat = sim_mat
# refer_gene_id = suggest_refer$id
norm_cpm_refer <- function(mat, refer_gene_id) {
	row_refer <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% refer_gene_id}
	if (!any(row_refer))
		stop('can\'t find any reference genes in the matrix for CPM normalization')

	norm_cpm_impl(mat, row_refer)
}
