#matrix_path = 'data-raw/external/scirep_sequential_qc.txt'
#classinfo_path = 'scirep_classes.txt'

#' @title read counts matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()]
#'
#' @return integer matrix
#'
#' @details In any case, first part (separated by `|`) of row names must be
#'   Ensembl transcript id
#'
#' @export

# path = 'data-raw/external/scirep_sequential_qc.txt'
read_mat <- function(path, ...) {
	path %>% readr::read_tsv(T, readr::cols(.default = 'c'), ...) %>%
		dplyr::mutate_at(-1, readr::parse_integer) %>%
		dplyr::rename('transcript' = 1) %>% as.data.frame() %>%
		tibble::column_to_rownames('transcript') %>% as.matrix()
}

#' @title sample classinfo
#'
#' @param path string.
#'
#' @return string matrix
#'
#' @details column 1 represents sample name, column 2 represents classinfo
#'
#' @export

# path = 'scirep_classes.txt'
read_classinfo <- function(path, ...) {
    read.table(path, sep = ",", header=T)
}

#' @title filter genes with low expression values
#'
#' @param mat integer matrix.
#' @param min_count, min_sample_per_gene integer scalar. For each gene, it must
#'   contain at least `min_count` reads in at least `min_sample_per_gene`
#'   samples. Otherwise, it would be dropped.
#'
#' @return integer matrix.
#'
#' @examples
#' filter_low(sim_mat)
#'
#' @export
filter_low <- function(mat, min_count = 2, min_sample_per_gene = 5) {
	low_per_row <- rowSums(mat > min_count)
	keeped_row <- low_per_row > min_sample_per_gene
	mat[keeped_row, ]
}

#' @imputation
#'
#' @param mat integer matrix.
#' @param tmp_path where tmp files stores, "data/expression_matrix/" for example.
#' @param out_path where outputs stores, "data/matrix_processing/imputation/" for example.
#' @param K imputation Kcluster
#' @param N imputation ncores
#' @return integer matrix named "scimpute_count.txt" stored in out_path.
#'
#' @examples imputation(mat, "data/expression_matrix/", "data/matrix_processing/imputation/",5,3)
#'
#' @export
imputation <- function(mat,tmp_path=".",impute_path="./imputation/", K = 5, N = 3) {
    suppressMessages(library("scImpute"))
    write.csv(mat, paste(tmp_path,"tmpsave.csv",sep=""))
    scimpute(count_path = paste(tmp_path,"tmpsave.csv",sep=""), infile = "csv",
    outfile = "txt", out_dir = impute_path , Kcluster = K, ncores = N)
    read.table(paste(out_path,"scimpute_count.txt",sep=""))
}


#' @export
plot_highest_exprs <- function(sce, top_n = 20) {
	sce %>% {suppressMessages(scater::calculateQCMetrics(.))} %>%
		scater::plotHighestExprs(n = top_n)
}


# plot_group --------------

plot_group_impl <- function(sce, shape = NULL, color = NULL, plot_fun) {
 	plot_fun(
 		sce,
		shape_by = shape, colour_by = color,
    	run_args = list(exprs_values = 'counts')
	)
}

#' @title plot PCA, TSNE
#'
#' @param sce A SingleCellExperiment object.
#' @param shape, color string. specify a column in `col_data` of [as_SingleCellExperiment()] to shape/color by
#'
#' @name plot_group
NULL


#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#' 
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_PCA <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotPCA)
}



#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#' 
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_TSNE <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotTSNE)
}

# plot CV -------------------------

coef_var_fun <- function(x) {
	sd(x, na.rm = T) / mean(x, na.rm = T)
}


#' @title density plot of coefficient of variation
#'
#' @param mat integer matrix. counts
#' @param refer_gene_id character.
#' @param refer_gene_name character.
#'
#' @return [ggplot2::ggplot()] object
#'
#' @examples NULL
#'
#' @export

# mat = sim_mat
# refer_gene_id = suggest_refer$id
# refer_gene_name = suggest_refer$name
plot_cv_density <- function(mat, refer_gene_id = '', refer_gene_name = refer_gene_id) {
	coef_var_df <- mat %>% apply(1, coef_var_fun) %>%
		{tibble::tibble(id = names(.), value = .)} %>%
		dplyr::mutate(id = stringr::str_extract(id, '[^|]+'))
	coef_var_refer_df <- tibble::tibble(id = refer_gene_id, name = refer_gene_name) %>%
		dplyr::inner_join(coef_var_df, by = 'id')

	plot <- ggplot2::ggplot(coef_var_df) +
		ggplot2::geom_density(ggplot2::aes(value), color = 'blue')

	if (nrow(coef_var_refer_df) > 0L) {

		plot = plot +
			ggplot2::geom_vline(xintercept = coef_var_refer_df$value, color = 'green') +
			ggplot2::geom_point(
				ggplot2::aes(x = value, y = seq_along(value)),
				data = coef_var_refer_df, size = 2, shape = 1
			) +
			ggrepel::geom_label_repel(
				ggplot2::aes(x = value, y = seq_along(value), label = name),
				data = coef_var_refer_df, hjust = 0.5
			)
	}

	return(plot)
}




