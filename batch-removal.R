norm_select <- function(
){
    #print norm_methods, user select
    print('Here are normalization methods, select one')
    print(c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'))
    #read norm_methods ->norm_methods
    norm_methods <- readline()
    if (norm_methods %in% c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'))
        return(norm_methods)
    #default SCnorm
    return('SCnorm')
}


batch <- function(
    mat,
    classinfo_path,
    batchinfo_path,
    norm_methods = 'SCnorm',
    batch_methods = c('ruv','combat'),
    output_path = './'
){
    suppressMessages(library(EDASeq))
    suppressMessages(library(RUVSeq))
    suppressMessages(library(sva))
    suppressMessages(library(scRNA.seq.funcs))
    
    norm_methods <- norm_select()
    #获取上一步输出结果
    mat <- get(paste('mat_',tolower(norm_methods),sep=''))
    mat_ruv <- ruv(mat,classinfo_path,output_path,2,10)
    mat_combat <- combat(mat,batchinfo_path,output_path,2)
}

ruv <- function(
    mat,
    classinfo_path,
    output_path,
    label_column = 2,
    k = 10
){
    cIdx <- rownames(mat,batchinfo_path,)
    
    sample_info <- read_classinfo(classinfo_path)
    ##根据mat排序
    rownames(sample_info) = sample_info$(names(sample_info)[1])
    sample_info=sample_info[names(mat),]
    rownames(sample_info) <- c()
    
    names(sample_info)[label_column]="label"
    scIdx <- matrix(-1, ncol = max(table(sample_info$label)), nrow = dim(table(sample_info$label)))
    labellist <- names(table(sample_info$label))
    for(i in c(1:dim(table(sample_info$label)))) {
        tmp <- which(sample_info$label == labellist[i])
        scIdx[i, 1:length(tmp)] <- tmp
    }
    ruv <- RUVs(as.matrix(mat), cIdx, k = k, scIdx = scIdx, isLog = TRUE)
    write.table(ruv$normalizedCounts, file=paste(output_path,'mx_ruv_batchremoval.txt',sep=''), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
    ruv$normalizedCounts
}

combat <- function(
    mat,
    batchinfo_path,
    output_path,
    batch_column = 2,
){
    batch_info <-read.csv(batchinfo_path,sep=',',row.names=1)
    batch_info=batch_info[names(m),]
    names(batch_info)[batch_column]="batch"
    mod <- model.matrix(~ 1, data = batch_info)
    combat <- ComBat(
        dat = log(mat+0.001),
        batch = factor(batch_info[,(batch_column-1)]),
        mod = mod,
        par.prior = TRUE,
        prior.plots = FLASE
    )
    mat <- exp(combat)
    write.csv(mat, file=paste(output_path,'mx_combat_batchremoval.txt',sep=''))
    mat
}
