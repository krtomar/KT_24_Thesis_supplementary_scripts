




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


my_TSS <- promoters(my_genes, upstream = 3000, downstream = 3000)


my_anti <- my_genes
strand(my_anti)[strand(my_genes) == "+"] <- "-"
strand(my_anti)[strand(my_genes) == "-"] <- "+"

my_TTS <- promoters(my_anti, upstream = 3000, downstream = 3000)
strand(my_TTS)[strand(my_genes) == "-"] <- "-"
strand(my_TTS)[strand(my_genes) == "+"] <- "+"


my_negatives <- start(my_TSS) < 0 | start(my_TTS) < 0

my_TSS <- my_TSS[!my_negatives]
my_TTS <- my_TTS[!my_negatives]

stopifnot(identical(my_TSS$gene_id, my_TTS$gene_id))


##################################################################################################################################
##################################################################################################################################



cov2mat <- function(my_cov, my_gr, id = "gene_id"){
        
        my_length <- unlist(lapply(my_cov, length))
        
        my_gr$length <- sapply(seqnames(my_gr), function(x){ my_length[names(my_length) == x]})
        my_gr <- my_gr[end(my_gr) < my_gr$length]

        my_mat <- matrix(as.numeric(unlist(my_cov[my_gr])), nrow = length(my_gr), byrow = TRUE)
        rownames(my_mat) <- mcols(my_gr)[,id]
        
        my_mat[as.character(strand(my_gr)) == "-",] <- t(apply(my_mat[as.character(strand(my_gr)) ==   "-", ], 1, rev))
        
        my_mat
}



##################################################################################################################################
##################################################################################################################################



norm_cov <- coverage(import(args[2]), weight = "score")


mat_TSS <- cov2mat(norm_cov, my_TSS, "gene_id")
mat_TTS <- cov2mat(norm_cov, my_TTS, "gene_id")


my_common_ids <- intersect(rownames(mat_TSS), rownames(mat_TTS))

mat_TSS <- mat_TSS[rownames(mat_TSS) %in% my_common_ids,]
mat_TTS <- mat_TTS[rownames(mat_TTS) %in% my_common_ids,]


stopifnot(identical(rownames(mat_TSS), rownames(mat_TTS)))


saveRDS(mat_TSS, file=args[3])
saveRDS(mat_TTS, file=args[4])



##################################################################################################################################
##################################################################################################################################


