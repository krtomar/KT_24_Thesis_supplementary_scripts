




##################################################################################################################################
##################################################################################################################################



library(ShortRead)
library(rtracklayer)


args = commandArgs(trailingOnly=TRUE)



##################################################################################################################################
##################################################################################################################################


my_gtf <- import(args[1])


my_genes <- my_gtf[my_gtf$type == "gene" & my_gtf$source %in% c("SGD", "PomBase","Ecoli")]
mcols(my_genes) <- mcols(my_genes)[,c("gene_id","gene_name")]

my_chromosomes <- seqlevels(my_genes)[grepl("scer", seqlevels(my_genes))]
my_genes <- keepSeqlevels(my_genes, value = my_chromosomes, pruning.mode = "coarse")


##################################################################################################################################
##################################################################################################################################


resizeVector <- function (x, out_length) 
{
        y <- spline(x = 1:length(x), y = x, n = out_length)$y
        return(y)
}



cov2mat <- function(my_cov, my_gr, id = "gene_id"){
        
        my_gr_left  <- resize(resize(my_gr, width = 1L, fix = "start", ignore.strand = TRUE), width = 1000L, fix = "end", ignore.strand = TRUE)
        my_gr_right <- resize(resize(my_gr, width = 1L, fix = "end",   ignore.strand = TRUE),   width = 1000L, fix = "start", ignore.strand = TRUE)
        
        my_length <- unlist(lapply(my_cov, length))
        my_gr$length <- sapply(seqnames(my_gr), function(x){ my_length[names(my_length) == x]})
        
        to_keep <- (start(my_gr_left) > 0) & (end(my_gr_right) < my_gr$length) 
        
        my_gr       <- my_gr[to_keep]
        my_gr_left  <- my_gr_left[to_keep]
        my_gr_right <- my_gr_right[to_keep]
        
        stopifnot(identical(my_gr_left$gene_id, my_gr_right$gene_id))
        stopifnot(identical(my_gr$gene_id, my_gr_right$gene_id))
        
        my_mat1 <- matrix(as.numeric(unlist(my_cov[my_gr_left])),  nrow = length(my_gr), byrow = TRUE)
        
        my_mat2 <- matrix(as.numeric(unlist(lapply(my_cov[my_gr], function(i){
                resizeVector(x = i, out_length = 2000L)}))),
                nrow = length(my_gr), byrow = TRUE)
        
        my_mat3 <- matrix(as.numeric(unlist(my_cov[my_gr_right])), nrow = length(my_gr), byrow = TRUE)
        
        my_mat <- cbind(my_mat1, my_mat2, my_mat3)
        
        rownames(my_mat) <- mcols(my_gr)[,id]
        
        my_mat[as.character(strand(my_gr)) == "-",] <- t(apply(my_mat[as.character(strand(my_gr)) ==   "-", ], 1, rev))
        
        my_mat
}



##################################################################################################################################
##################################################################################################################################



norm_cov <- coverage(import(args[2]), weight = "score")


mat_GB <- cov2mat(norm_cov, my_genes, "gene_id")


saveRDS(mat_GB, file=args[3])



##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################
