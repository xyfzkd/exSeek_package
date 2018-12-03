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

# plot_variance --------------


#' get y axis range of ggplot object.
get_y_range <- function(plot) {
    ggplot2::ggplot_build(plot)$layout$panel_params[[1]]$y.range
}

#' @title generate equally spaced y coordinates, not hit bottom nor top.
#'
#' @details Image `plot`'s y axis extends from 0 to 3, `x` contains 3 values,
#' then we give `c(0.5, 1.5, 2.5)`.
#'
#' @param plot ggplot object.
#' @param x numberic vector
#'
#' @return numeric vector. The same length as `x`
#'
#' @keywords internal
seq_y <- function(plot, x) {
    y_range <- get_y_range(plot)
    by <- diff(y_range) / length(x)
    seq(y_range[1] + by, y_range[2], by) - by/2
}


#' @title plot variance of counts across samples
#'
#' @param mat integer matrix. counts
#' @param refer_gene_id character. Ensembl transcript id, add a vertical line
#'   for each gene to mark the corresponding CV (on x axis). Only genes in
#'   counts matrix would be shown. Usually these genes should be the reference
#'   genes you want to use for normalization.
#' @param refer_gene_name character. Transcript name
#'
#' @return [ggplot2::ggplot()] object.
#'
#' @name plot_variance
NULL





#' @details
#' `plot_cv_density()` produces density plot of coefficient of variation
#'
#' @examples
#' # only one gene exist in the matrix
#' plot_cv_density(sim_mat, suggest_refer$id)
#' plot_cv_density(sim_mat, suggest_refer$id, suggest_refer$name)
#'
#' # the name should be the same length as id
#' plot_cv_density(sim_mat, rownames(sim_mat)[1:6], letters[1:3])
#' # if only part of the genes have name, you can pass the id of other genes
#' plot_cv_density(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))
#'
#' @export
#'
#' @rdname plot_variance

# mat = sim_mat
# refer_gene_id = suggest_refer$id
# refer_gene_name = suggest_refer$name
plot_cv_density <- function(mat, refer_gene_id = '', refer_gene_name = refer_gene_id) {
    cv <- mat %>% apply(1, cv_fun) %>%
    {tibble::tibble(id = names(.), value = .)} %>%
    dplyr::mutate(id = stringr::str_extract(id, '[^|]+'))
    plot <- ggplot2::ggplot(cv, ggplot2::aes(value)) +
    ggplot2::geom_density(color = 'blue') +
    ggplot2::labs(x = 'coefficient of variation')
    
    if (length(refer_gene_id) != length(refer_gene_name)) {
        warning("Ignoring refer_gene_name, since it isn't the same length as refer_gene_id")
        refer_gene_name = refer_gene_id
    }
    cv_refer <- tibble::tibble(id = refer_gene_id, name = refer_gene_name) %>%
    dplyr::inner_join(cv, by = 'id')
    if (nrow(cv_refer) == 0L) {
        warning("None refer gene found in the count matrix")
        return(plot)
    }
    
    plot + ggplot2::geom_vline(xintercept = cv_refer$value, color = 'green') +
    ggplot2::geom_point(
    ggplot2::aes(x = value, y = seq_y(plot, value)),
    data = cv_refer, size = 2, shape = 1
    ) +
    ggrepel::geom_label_repel(
    ggplot2::aes(x = value, y = seq_y(plot, value), label = name),
    data = cv_refer, hjust = 0.5
    )
}



#' @details
#' `plot_cv_density()` produces density plot of coefficient of variation
#'
#' @examples
#' # only one gene exist in the matrix
#' plot_refer_violin(sim_mat, suggest_refer$id)
#' plot_refer_violin(sim_mat, suggest_refer$id, suggest_refer$name)
#'
#' # the name should be the same length as id
#' plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], letters[1:3])
#' # if only part of the genes have name, you can pass the id of other genes
#' plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))
#'
#' @export
#'
#' @rdname plot_variance

# mat = sim_mat
# refer_gene_id = rownames(mat)[1:6]
# refer_gene_name = paste0('gene_', letters[1:6])
plot_refer_violin <- function(mat, refer_gene_id, refer_gene_name = refer_gene_id) {
    if (length(refer_gene_id) != length(refer_gene_name)) {
        warning("Ignoring refer_gene_name, since it isn't the same length as refer_gene_id")
        refer_gene_name = refer_gene_id
    }
    
    refer_gene <- tibble::tibble(id = refer_gene_id, name = refer_gene_name)
    refer_count <- mat %>% tibble::as_tibble(rownames = 'id') %>%
    dplyr::mutate(id = stringr::str_extract(id, '[^|]+')) %>%
    dplyr::inner_join(refer_gene, ., by = 'id') %>% dplyr::select(-id)
    if (nrow(refer_count) == 0L) {
        warning('None refer gene found in the count matrix')
        return(ggplot2::ggplot())
    }
    
    refer_count_long <- refer_count %>%    tidyr::gather('sample', 'count', -1) %>%
    dplyr::mutate_at('name', as.factor)
    g_violin <- refer_count_long %>%
    ggplot2::ggplot(ggplot2::aes(name, log2(count + 0.001))) +
    ggplot2::geom_violin() +
    ggplot2::labs(x = 'reference transcripts', y = quote(log[2](count)))
    
    # max y coordinate of each violin
    y_max <- ggplot2::ggplot_build(g_violin)$data[[1]] %>% tibble::as_tibble() %>%
    dplyr::group_by(x) %>% dplyr::arrange(dplyr::desc(y)) %>% dplyr::slice(1) %>%
    dplyr::ungroup() %>% dplyr::arrange(x) %>% dplyr::select(x, y)
    
    cv_df <- refer_count_long %>%
    dplyr::group_by(name) %>% dplyr::summarise(cv = cv_fun(count)) %>%
    dplyr::arrange(name) %>% dplyr::mutate(x = seq_along(name)) %>%
    dplyr::inner_join(y_max, by = 'x') %>%
    dplyr::mutate(y = y + diff(get_y_range(g_violin)) / 20) %>%
    dplyr::mutate(cv = formatC(cv, digits = 3, format = 'f'))
    
    g_violin + ggplot2::geom_text(ggplot2::aes(x, y, label = cv), cv_df, color = 'blue')
}
