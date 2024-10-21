






##################################################################################################################################
##################################################################################################################################




library(DESeq2)
library(sva)


library(ShortRead)
library(rtracklayer)
library(GenomicFeatures)

library(HelpersforDESeq2)

library(org.Sc.sgd.db)
library(topGO)

library(RColorBrewer)
library(pheatmap)


library(ggplot2)

##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)


input_anno <- grep("gtf", args, value = TRUE)
input_table <- grep("SampleTable", args, value = TRUE)

norm_type <- grep("^Scer$|^Spom$", args, value = TRUE)

output_tables <- gsub("\\/log2.*","", grep("log2_", args, value = TRUE))
output_plots <-  gsub("\\/PCA.*","", grep("PCA", args, value = TRUE))
output_GO    <-  gsub("plots/PCA.*","GO",  grep("PCA", args, value = TRUE))

output_sessionInfo <- grep("sessionInfo", args, value = TRUE)





##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



my_gtf <- import(input_anno)


my_genes <- my_gtf[my_gtf$type == "gene" & my_gtf$source %in% c("SGD", "PomBase","Ecoli")]



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        SampleTable         #####################################################




SampleTable <- read.delim(input_table, header = T, stringsAsFactors = F)


SampleTable <- SampleTable[!duplicated(SampleTable$LibraryID), , drop = FALSE]
SampleTable <- SampleTable[order(SampleTable$LibraryID), , drop = FALSE]


SampleTable$Dataset <- gsub(" ","_", SampleTable$Dataset)

SampleTable$Conditions <- paste(SampleTable$Genotype, gsub(" .*", "", SampleTable$Fraction), sep = "_")

SampleTable$fitGroup <- paste(SampleTable$Genotype, SampleTable$Dataset, sep = "_")

##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        count table         #####################################################



my_count_files <-  grep("ReadsPerGene.out.tab", args, value = TRUE)
my_count_files <- my_count_files[order(my_count_files)]

if(identical(SampleTable$LibraryID, gsub(".*\\/|.ReadsPerGene.out.tab|\\..*","",my_count_files))){
        
        my_count_table <- makeCountTable(count_names = SampleTable$LibraryID,
                                         count_files = my_count_files, 
                                         stranded = FALSE)
}

my_count_table <- my_count_table[-1:-4,]



########################################################################################################################################################################  
########################################################################################################################################################################  




my_genes_Scer <- my_genes[grep("scer|ecol", as.character(seqnames(my_genes)))]
my_genes_Spom <- my_genes[grep("spom", as.character(seqnames(my_genes)))]



my_counts_Scer <- my_count_table[rownames(my_count_table) %in% c(my_genes_Scer$gene_id, "Ecoli"),]
my_counts_Spom <- my_count_table[rownames(my_count_table) %in% my_genes_Spom$gene_id,]


write.table(my_counts_Scer, file = paste0(output_tables,"/raw_counts_Scer.txt"), quote = F, sep = "\t", row.names = TRUE, col.names = NA)
write.table(my_counts_Spom, file = paste0(output_tables,"/raw_counts_Spom.txt"), quote = F, sep = "\t", row.names = TRUE, col.names = NA)


########################################################################################################################################################################  
########################################################################################################################################################################  





















##################################################################################################################################
##################################################################################################################################

#################################################      Quality  Control      ##################################################### 




######### cutoffs #########

min_genic_reads = 5e4
min_detected_genes = 1e3


################################################# 


n_all_Scer_reads <- colSums(my_counts_Scer)
n_all_Spom_reads <- colSums(my_counts_Spom)


n_genic_Scer_reads <- colSums(my_counts_Scer[!grepl("RDN|SPRRNA", rownames(my_counts_Scer)),])
n_genic_Spom_reads <- colSums(my_counts_Spom[!grepl("RDN|SPRRNA", rownames(my_counts_Spom)),])
n_genic_Ecoli_reads <- my_counts_Scer[rownames(my_counts_Scer) == "Ecoli",]


n_rRNA_Scer_reads <- colSums(my_counts_Scer[grepl("RDN|SPRRNA", rownames(my_counts_Scer)),])
n_rRNA_Spom_reads <- colSums(my_counts_Spom[grepl("RDN|SPRRNA", rownames(my_counts_Spom)),])


n_detected_Scer_genes <- apply(my_counts_Scer, 2, function(x){   sum(x > 0)   })
n_detected_Spom_genes <- apply(my_counts_Spom, 2, function(x){   sum(x > 0)   })



my_outliers <- (n_genic_Scer_reads < min_genic_reads) | 
        (n_detected_Scer_genes < min_detected_genes) 




##################################################################################################################################
##################################################################################################################################




my_color_palette <- c("#E69F00","#999999","#56B4E9","#009E73")

if(all(identical(colnames(my_counts_Scer), SampleTable$LibraryID),
       identical(colnames(my_counts_Spom), SampleTable$LibraryID))){
        
        my_conditions <- factor(SampleTable$Conditions, levels = unique(SampleTable$Conditions))
}



################################################# 


pdf(paste0(output_plots,"/barplot_QC.pdf"), height = 7, width = 14, useDingbats = F)
par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


bp <- barplot(n_genic_Scer_reads / 1e6, 
              #ylim = c(0,1.5), 
              ylab = "Number of Genic Scer Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_genic_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


bp <- barplot(n_genic_Spom_reads / 1e6, 
              #ylim = c(0,1), 
              ylab = "Number of Genic Spom Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_genic_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 


bp <- barplot((n_genic_Scer_reads/(n_genic_Spom_reads+n_genic_Scer_reads)), 
              #ylim = c(0,1), 
              ylab = "Genic Ratio Scer/(Scer+Spom)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_genic_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 




bp <- barplot(n_rRNA_Scer_reads / 1e6, 
              #ylim = c(0,10), 
              ylab = "Number of rRNA Scer Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_rRNA_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


bp <- barplot(n_rRNA_Spom_reads / 1e6, 
              #ylim = c(0,10), 
              ylab = "Number of rRNA Spom Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_rRNA_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 


bp <- barplot((n_rRNA_Scer_reads/(n_rRNA_Spom_reads+n_rRNA_Scer_reads)), 
              #ylim = c(0,1), 
              ylab = "rRNA Ratio Scer/(Scer+Spom)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_rRNA_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)

plot.new()

legend("bottom", legend = levels(my_conditions), horiz = TRUE, 
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[seq_along(levels(my_conditions))])



################################################# 

par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


bp <- barplot(n_all_Scer_reads / 1e6, 
              #ylim = c(0,10), 
              ylab = "Number of all Scer Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_all_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


bp <- barplot(n_all_Spom_reads / 1e6, 
              #ylim = c(0,10), 
              ylab = "Number of all Spom Reads (Million)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_all_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)


################################################# 


bp <- barplot((n_all_Scer_reads/(n_all_Spom_reads+n_all_Scer_reads)), 
              #ylim = c(0,1), 
              ylab = "all Reads Ratio Scer/(Scer+Spom)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


#abline(h = min_all_reads / 1e6, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 


par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0), new = TRUE)

plot.new()

legend("center", legend = levels(my_conditions), horiz = TRUE, 
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[seq_along(levels(my_conditions))])



################################################# 


par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))

bp <- barplot(n_detected_Scer_genes,
              #ylim = c(0,8e3), 
              ylab = "Number of Detected Scer Genes (Read Counts > 0)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)

#abline(h = min_detected_genes, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)



################################################# 

bp <- barplot(n_detected_Spom_genes,
              #ylim = c(0,8e3), 
              ylab = "Number of Detected Spom Genes (Read Counts > 0)", las = 3, xaxt = "n",
              col = my_color_palette[my_conditions])

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)

#abline(h = min_detected_genes, lty=2)

my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
my_outlier_colors[my_outliers] <- "red3"

par(xpd=NA)
points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
par(xpd = FALSE)




################################################# 


par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0), new = TRUE)

plot.new()

legend("center", legend = levels(my_conditions), horiz = TRUE, 
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[seq_along(levels(my_conditions))])


################################################# 


par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


bp_dat <- rbind(n_genic_Ecoli_reads,
                n_all_Spom_reads,
                n_all_Scer_reads)

bp <- barplot( t(t(bp_dat)/colSums(bp_dat)), 
               #ylim = c(0,1), 
               ylab = "Fraction of all Reads (rel. to all sum)", las = 3, xaxt = "n")

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)



bp_dat <- rbind(n_genic_Ecoli_reads,
                n_genic_Spom_reads,
                n_genic_Scer_reads)

bp <- barplot( t(t(bp_dat)/colSums(bp_dat)), 
               #ylim = c(0,1), 
               ylab = "Fraction of genic Reads (rel. to genic sum)", las = 3, xaxt = "n")

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


bp_dat <- rbind(n_genic_Spom_reads,
                n_genic_Ecoli_reads,
                n_genic_Scer_reads)

bp <- barplot( t(t(bp_dat)/n_all_Spom_reads), col = c("#AEAEAE", "#4D4D4D", "#E6E6E6"),
               #ylim = c(0,0.5), 
               ylab = "Fraction of genic Reads (rel. to all Spom)", las = 3, xaxt = "n")

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


bp_dat <- rbind(n_genic_Spom_reads,
                n_genic_Ecoli_reads,
                n_genic_Scer_reads)

bp <- barplot( t(t(bp_dat)/n_genic_Spom_reads), col = c("#AEAEAE", "#4D4D4D", "#E6E6E6"),
               #ylim = c(0,5), 
               ylab = "Fraction of genic Reads (rel. to genic Spom)", las = 3, xaxt = "n")

axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)


################################################# 


par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0), new = TRUE)

plot.new()

legend("center", legend = c("Scer","Spom","Ecoli"),  horiz = TRUE, 
       bty = "n", cex = 0.8, pt.cex = 1, fill = rev(gray.colors(nrow(bp_dat))))



################################################# 


my_fav_genes <- c("Ecoli","ISW1","CHD1","ISW2","ACT1","PGK1","PMA1","RAD51")


for(norm_set in c("genic","all")){
        
        par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
        
        for(my_fav_gene in my_fav_genes){
                
                numerator <- my_counts_Scer[rownames(my_counts_Scer) %in% my_genes$gene_id[grep(my_fav_gene, my_genes$gene_name)],] 
                denominator <- get(paste0("n_",norm_set,"_",norm_type,"_reads"))
                
                stopifnot(identical(names(numerator), names(denominator)))
                stopifnot(length(numerator) == length(denominator))
                stopifnot(identical(colnames(my_counts_Scer), names(denominator)))
                
                CPM_gene <- numerator / denominator * 1e6 
                
                
                bp <- barplot(CPM_gene, 
                              main = my_fav_gene, ylim = c(0, max(CPM_gene)*1.25),
                              ylab = paste("Counts Per Million",norm_set,norm_type, "Reads"), las = 3, xaxt = "n",
                              col = my_color_palette[my_conditions])
                
                axis(side = 1, at = bp, labels = colnames(my_counts_Scer), cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)
                
                
                
                my_outlier_colors <- rep(rgb(1,1,1,0), length(my_conditions))
                my_outlier_colors[my_outliers] <- "red3"
                
                par(xpd=NA)
                points(x = bp, y = rep(0,length(bp)), pch=17, col = my_outlier_colors)
                par(xpd = FALSE)
                
                #abline(h = min_CPM_CHD1, lty=2)
        }
        
        ################################################# 
        
        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
        
        plot.new()
        
        legend("bottom", legend = levels(my_conditions), horiz = TRUE, 
               bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[seq_along(levels(my_conditions))])
        
        
        ################################################# 
        
}




dev.off()





##################################################################################################################################
##################################################################################################################################






















##################################################################################################################################
##################################################################################################################################

#################################################           Setup            #####################################################




if(all(identical(colnames(my_counts_Scer), SampleTable$LibraryID),
       identical(colnames(my_counts_Spom), SampleTable$LibraryID),
       identical(colnames(my_counts_Scer), names(n_all_Spom_reads)))){
        
        SampleTable$Group <-  paste(SampleTable$Conditions, SampleTable$Time, sep = "_")
        
        dds <- setupDDS(CountTableName = "my_counts_Scer",
                        SampleTableName = "SampleTable",
                        SampleIdName = "LibraryID",
                        ConditionName = "Group",
                        BatchName = "Dataset",
                        n_samples_for_filtering = ceiling(ncol(my_counts_Scer)*0.5)-1,
                        min_number_of_reads = 1)
        
        if(norm_type == "Spom"){
                
                #sF_Spom <- estimateSizeFactorsForMatrix(my_counts_Spom)
                sF_Spom <- n_all_Spom_reads/mean(n_all_Spom_reads)
                
                my_norm_counts_tmp <- t(t(counts(dds))/sF_Spom)
                
                n_sums <- colSums(my_norm_counts_tmp)
                
                my_fitted_vals <- sapply(names(n_sums), function(j){
                        
                        my_fit_group <- SampleTable$fitGroup ==  SampleTable$fitGroup[SampleTable$LibraryID == j]
                        
                        my_ids_pull  <- SampleTable$LibraryID[my_fit_group & SampleTable$Fraction == "pull out fraction"]
                        my_ids_total <- SampleTable$LibraryID[my_fit_group & SampleTable$Fraction == "total RNA"]
                        
                        stopifnot(length(my_ids_pull) == 7)
                        stopifnot(length(my_ids_total) == 7)
                        stopifnot(identical(SampleTable$Time[match(my_ids_pull,  SampleTable$LibraryID)],
                                            SampleTable$Time[match(my_ids_total, SampleTable$LibraryID)]))
                        
                        n_totals <- n_sums[match(my_ids_total, names(n_sums))]
                        t_time   <- SampleTable$Time[match(names(n_totals), SampleTable$LibraryID)]
                        
                        fit <- nls(n_totals ~ x0 * exp(k * t_time),
                                   start = c(x0 = min(n_totals), k = 0.001),
                                   control = nls.control(maxiter = 200))
                        
                        my_fitted <- fitted(fit) 
                        names(my_fitted) <- names(n_totals)
                        
                        my_fitted <- my_fitted/mean(my_fitted)
                        
                        if(j %in% my_ids_pull){
                                return(my_fitted[my_ids_pull == j])
                        }else{
                                return(my_fitted[names(my_fitted) == j])
                        }
                })
                
                my_test <- as.numeric(gsub("L|.*\\.","", names(my_fitted_vals))) - as.numeric(gsub("L|\\..*","", names(my_fitted_vals)))  
                stopifnot(sum(my_test == 0) == 42)
                stopifnot(sum(my_test == 14) == 42)
                
                my_fitted_vals <- my_fitted_vals
                names(my_fitted_vals) <- gsub("\\..*","", names(my_fitted_vals))
                
                stopifnot(identical(names(my_fitted_vals), names(sF_Spom)))
                
                sizeFactors(dds) <- (sF_Spom*my_fitted_vals)
                dds <- DESeq(dds)
                
                rm(list = "my_norm_counts_tmp")
                rm(list = "my_fitted_vals")
        }
}





##################################################################################################################################
##################################################################################################################################

#################################################         Contrasts          #####################################################




contrast_list <- c(lapply(1:(nrow(SampleTable)/2), function(i){c(SampleTable$Group[grep("CHD1", SampleTable$Group)][i],
                                                                 SampleTable$Group[grep("WT", SampleTable$Group)][i])}),
                   lapply(1:(nrow(SampleTable)/2), function(i){c(SampleTable$Group[grep("WT", SampleTable$Group)][i],
                                                                 gsub("total_.*","total_0", gsub("pull_.*","pull_0",
                                                                                                   SampleTable$Group[grep("WT", SampleTable$Group)][i])))}),
                   lapply(1:(nrow(SampleTable)/2), function(i){c(SampleTable$Group[grep("CHD1", SampleTable$Group)][i],
                                                                 gsub("total_.*","total_0", gsub("pull_.*","pull_0",
                                                                                                   SampleTable$Group[grep("CHD1", SampleTable$Group)][i])))}))



for(contrast_name in contrast_list){
        
        if(contrast_name[1] == contrast_name[2]){next()}
        
        res <- getResults(dds = dds,
                          contrast = contrast_name,
                          lfc_cutoff = 0,
                          shrink = F,
                          annotation = "my_genes",
                          anno_symbol = "gene_name",
                          anno_id = "gene_id")
        
        res$gene_id <- rownames(res)
        
        res_name <- paste0("res.", contrast_name[1],"-",contrast_name[2])
        assign(res_name, res)
        
        write.table(res, file = paste0(output_tables,"/",res_name, ".txt"),
                    quote = F, sep = "\t", row.names = T, col.names = NA)
        
        rm(list = "res")
        
}

##################################################################################################################################
##################################################################################################################################
















##################################################################################################################################
##################################################################################################################################

#################################################          MA-plots          #####################################################


res_names <- ls(pattern = "^res\\.")

padj_cutoff <- 0.05


#####################################################


dir.create(paste0(output_plots, "/MAplots"))

NoLabel <- ""
RiboProtGenes <- grep("^RPL|^RPS", my_genes$gene_name, value = TRUE)


for(my_label in c("NoLabel","Labeled","RiboProtGenes")){
        
        for(res_name in res_names){
                
                if(grepl("total", res_name)){
                        Labeled <- c("ACT1", "CHD1", "PGK1","PMA1","RAD51")
                } else{
                        Labeled <- c("ACT1", "CHD1", "PGK1","PMA1","RAD51", "Ecoli")
                }  
                
                
                pdf(paste0(output_plots,"/MAplots","/MA_",res_name,"_",my_label,".pdf"), height = 6, width = 6, useDingbats = F)
                par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                
                res <- get(res_name)
                
                plottingMA(res = res,
                           main_title = gsub("res.","",res_name),
                           main_title_size = 1,
                           selection_ids = get(my_label),
                           selection_id_type = "gene_name",
                           selection_name = my_label,
                           selection_point_size = 0.5,
                           selection_text_label = (my_label == "Labeled"),
                           selection_text_adj =  1.5,
                           selection_text_size = 0.8,
                           selection_shadow = FALSE,
                           xlims = c(0, 5),
                           ylims = c(-10,10),
                           x_axis_by = 1,
                           padj_cutoff = padj_cutoff,
                           show_legend = (my_label != "NoLabel"))
                
                if(my_label == "NoLabel"){
                        
                        legend("topright", 
                               legend = c("all significant", 
                                          "all non-significant"), 
                               col = c(rgb(0.9, 0.6, 0, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), 
                               bg = "white", border = NA, bty = "n", 
                               cex = 0.8, pch = 19)
                        
                        legend("bottomright", 
                               legend = c(paste0("sign. up n=",  sum(res$padj < padj_cutoff & res$log2FoldChange > 0)), 
                                          paste0("sign. down n=",sum(res$padj < padj_cutoff & res$log2FoldChange < 0))), 
                               bg = "white", border = NA, bty = "n", 
                               cex = 0.8, pch = NA)
                }
                
                
                dev.off()
        }
        
}


##################################################################################################################################
##################################################################################################################################


































##################################################################################################################################
##################################################################################################################################

#################################################             rld            #####################################################




rld <- rlog(dds, blind = FALSE)

rlog_norm_counts_raw <- assay(rld)

log2_norm_counts_raw <- log2(counts(dds, normalized = TRUE)+1)



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################      Batch correction      #####################################################



batchVar <- colData(dds)$Batch

modcombat <- model.matrix(~Sample, data = colData(dds))



rlog_norm_counts <- ComBat(dat = rlog_norm_counts_raw,
                           batch = batchVar, mod = modcombat,
                           par.prior = TRUE, prior.plots = FALSE)



log2_norm_counts <- ComBat(dat = log2_norm_counts_raw,
                           batch = batchVar, mod = modcombat,
                           par.prior = TRUE, prior.plots = FALSE)


# rlog_norm_counts <- rlog_norm_counts_raw
#
# log2_norm_counts <- log2_norm_counts_raw

write.table(log2_norm_counts, file = paste0(output_tables,"/log2_norm_counts.txt"),
            quote = F, sep = "\t", row.names = T, col.names = NA)



##################################################################################################################################
##################################################################################################################################




stopifnot(identical(SampleTable$Conditions, as.character(my_conditions)))
stopifnot(identical(SampleTable$LibraryID, colnames(log2_norm_counts)))

annotation_col = data.frame(Conditions = my_conditions)
rownames(annotation_col) <- colnames(log2_norm_counts)

ann_colors = list(Conditions = my_color_palette[seq_along(levels(my_conditions))])
names(ann_colors$Conditions) <- levels(my_conditions)




##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################

#################################################        Correlation         #####################################################



callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}


#################################################


pdf(paste0(output_plots,"/correlation.pdf"), width = 20, height = 20, useDingbats = FALSE)
par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)


for(log2_type in c("log2", "rlog")){
        
        
        my_log2_data <- get(paste0(log2_type,"_norm_counts"))
        
        my_log2_scaled <- my_log2_data - rowMeans(my_log2_data)
        
        
        #################################################
        
        for(cor_type in c("Spearman","Pearson")){
                
                my_corr <- cor(my_log2_scaled, method = tolower(cor_type))
                rownames(my_corr) <- SampleTable$Group[match(rownames(my_corr), SampleTable$LibraryID)]
                
                pheatmap(my_corr,
                         main = paste0(cor_type, "´s R at genes n = ", nrow(my_log2_scaled)," (",log2_type," counts)"),
                         color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100),
                         #breaks = seq(0.92,1, length.out = 101),
                         cellwidth = 12, cellheight = 12,
                         clustering_callback = callback,
                         clustering_method = "ward.D2",
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors)
        }
        
        #################################################
        
        
}


dev.off()



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

#################################################            PCA             #####################################################



xycomps <- list(c(1,2),
                c(1,3),
                c(2,3))


my_group <- factor(SampleTable$Group, levels = unique(SampleTable$Group))
my_group_colors <- paste0(rep(my_color_palette,each=7), rep(c("1A","50","66","8C","B3","D9","FF"), 4))


stopifnot(identical(levels(my_conditions), unique(gsub("_[0-9].*","",levels(my_group)))))
stopifnot(identical(my_color_palette[my_conditions], substr(my_group_colors[my_group], 1, 7)))




pdf(paste0(output_plots,"/PCA.pdf"), height = 3.5, width = 9, useDingbats = F)
par(mfrow=c(1,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))


for(log2_type in c("log2", "rlog")){
        
        
        my_log2_data <- get(paste0(log2_type,"_norm_counts"))
        
        
        #################################################
        
        par(mfrow=c(1,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
        
        for(xycomp in xycomps){
                
                plottingPCA(my_log2_data,
                            xcomp = xycomp[1],
                            ycomp = xycomp[2],
                            conditions = my_group,
                            pca_colors = my_group_colors,
                            main_title = "PCA",
                            quantiles = c(0,1),
                            show_labels = FALSE,
                            point_size = 1,
                            my_xlimits = c(-1,1)*ifelse(xycomp[1] == 2, 40, 200),
                            my_ylimits = c(-40, 40)
                )
                
        }
        
        mtext(text = paste(log2_type, "counts"), side = 1, outer = TRUE, adj = 1, cex = 0.75)
        
        #################################################
        
        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
        
        plot.new()
        
        legend("bottom", legend = levels(my_conditions), horiz = TRUE, 
               bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[seq_along(levels(my_conditions))])
        
        
        
        
        #################################################
}

dev.off()


##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

#################################################        Top Hits         #####################################################



dir.create(paste0(output_plots, "/Heatmaps"))


n_hits <- 50




for(log2_type in c("log2", "rlog")){
        
        
        my_log2_data <- get(paste0(log2_type,"_norm_counts"))
        
        
        #################################################
        
        for(res_name in res_names){
                
                
                pdf(paste0(output_plots,"/Heatmaps/Heatmaps.",gsub("res.","",res_name) ,".", log2_type,".pdf"), 
                    width = 18, height = 11, useDingbats = FALSE)
                par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
                
                
                
                res <- get(res_name)
                res <- res[order(res$pvalue),]
                
                my_fav_genes <- c(res$gene_id[res$log2FoldChange > 0][1:(n_hits/2)],
                                  res$gene_id[res$log2FoldChange < 0][1:(n_hits/2)])
                
                
                my_hm_data <- my_log2_data[rownames(my_log2_data) %in% my_fav_genes,]
                my_hm_data <- my_hm_data - rowMeans(my_hm_data)
                
                my_title <- paste0("Top",(n_hits/2), " up, ","Top",(n_hits/2), " down"," in ", gsub("res.","",res_name))
                
                
                my_hm_rownames <- my_genes$gene_name[match(rownames(my_hm_data),my_genes$gene_id)]
                rownames(my_hm_data)[!is.na(my_hm_rownames)] <- my_hm_rownames[!is.na(my_hm_rownames)]
                
                rm(list = "res")
                rm(list = "my_fav_genes")
                
                pheatmap(my_hm_data,
                         main = paste0(my_title," (relative ",log2_type," counts)"),
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100),
                         breaks = seq(-4,4, length.out = 101),
                         clustering_callback = callback,
                         clustering_method = "ward.D2",
                         cellwidth = 12, cellheight = 12)
                
                rm(list = "my_hm_data")
                
                
                dev.off()
        }
        
        #################################################
}




##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################



GSE69400_data <- read.delim("../GSE69400/analysis/ChIP.normalized_counts.txt")

GSE69400_data <- GSE69400_data[GSE69400_data$name %in% my_genes_Scer$gene_id[my_genes_Scer$gene_biotype == "protein_coding"] ,]

GSE69400_data$log2FC.Pol2_vs_Input_WT <- rowMeans(log2(GSE69400_data[,c("SRR2045644_WT_Rpb3_IP",    "SRR2045645_WT_Rpb3_IP")]+1)) - 
        rowMeans(log2(GSE69400_data[,c("SRR2045642_WT_Rpb3_Input", "SRR2045643_WT_Rpb3_Input")]+1))


GSE69400_data$Qs.Pol2_WT <-  c("q1","q2","q3","q4","q5")[cut(x = GSE69400_data$log2FC.Pol2_vs_Input_WT,
                                                             breaks = quantile(GSE69400_data$log2FC.Pol2_vs_Input_WT, probs = seq(0, 1, 0.20)),
                                                             include.lowest = TRUE)]


##################################################################################################################################
##################################################################################################################################



GSE90998_data <- read.delim("../GSE90998//analysis_chip//ChIP.normalized_counts.txt")

GSE90998_data <- GSE90998_data[GSE90998_data$name %in% my_genes_Scer$gene_id[my_genes_Scer$gene_biotype == "protein_coding"] ,]

GSE90998_data$Qs.CHD1_WT <-  c("q1","q2","q3","q4","q5")[cut(x = GSE90998_data$log2FC.CHD1_vs_mock,
                                                             breaks = quantile(GSE90998_data$log2FC.CHD1_vs_mock, probs = seq(0, 1, 0.20)),
                                                             include.lowest = TRUE)]


##################################################################################################################################
##################################################################################################################################



RNA_data <- read.delim("../GSE90998/analysis_rna//RNA.normalized_counts.txt")

RNA_data <- RNA_data[RNA_data$name %in% my_genes_Scer$gene_id[my_genes_Scer$gene_biotype == "protein_coding"] ,]

RNA_data$RNA_WT <- rowMeans(log2(RNA_data[,c("SRR5085167_wildtype", "SRR5085168_wildtype")]+1)) 

RNA_data$Qs.RNA_WT <-  c("q1","q2","q3","q4","q5")[cut(x = RNA_data$RNA_WT,
                                                       breaks = quantile(RNA_data$RNA_WT, probs = seq(0, 1, 0.20)),
                                                       include.lowest = TRUE)]

RNA_data$Qs.RNA_log2FC_dCHD1vsWT <-  c("q1","q2","q3","q4","q5")[cut(x = RNA_data$log2FC.delCHD1_vs_WT,
                                                                     breaks = quantile(RNA_data$log2FC.delCHD1_vs_WT, probs = seq(0, 1, 0.20)),
                                                                     include.lowest = TRUE)]

##################################################################################################################################
##################################################################################################################################



elife_halflife <- read.delim("external/elife-32536-fig1-data2-v4.txt")

elife_halflife <- elife_halflife[complete.cases(elife_halflife),]

elife_halflife$name <- elife_halflife$gene_id
elife_halflife$halflife <- rowMeans(elife_halflife[,c("halflife_160412_r1", "halflife_160412_r2")])

elife_halflife$Qs.halflife <-  c("q1","q2","q3","q4","q5")[cut(x = elife_halflife$halflife,
                                                               breaks = quantile(elife_halflife$halflife, probs = seq(0, 1, 0.20)),
                                                               include.lowest = TRUE)]


##################################################################################################################################
##################################################################################################################################


ext_data <- merge(GSE69400_data, RNA_data, by = "name")


ext_data <- merge(ext_data, GSE90998_data, by = "name")


ext_data <- merge(ext_data, elife_halflife, by = "name")


##################################################################################################################################
##################################################################################################################################




log2_gene_length <- log2(width(my_genes))
names(log2_gene_length) <- my_genes$gene_id


log2_gene_length <-  log2_gene_length[match(rownames(log2_norm_counts), names(log2_gene_length)) ]

scaled_log2_gene_length <- log2_gene_length - mean(log2_gene_length)


ext_data <- merge(ext_data, 
                  data.frame(log2_gene_length = log2_gene_length,
                             name = names(log2_gene_length)), by = "name")

# ext_data$Qs.gene_length <-  c("q1","q2","q3","q4","q5")[cut(x = ext_data$log2_gene_length,
#                                                             breaks = quantile(ext_data$log2_gene_length, probs = seq(0, 1, 0.20)),
#                                                             include.lowest = TRUE)]


ext_data$Qrs.gene_length <- cut(x = (2^ext_data$log2_gene_length)/1e3,
                                breaks = quantile((2^ext_data$log2_gene_length)/1e3, probs = seq(0, 1, 0.20)),
                                include.lowest = TRUE)



ext_data$Qs.gene_length <- c("q1","q2","q3","q4","q5")[ext_data$Qrs.gene_length]



ext_data <- ext_data[,c("name", "commonName", 
                        "log2FC.Pol2_vs_Input_WT","Qs.Pol2_WT", 
                        "log2FC.CHD1_vs_mock",    "Qs.CHD1_WT",
                        "RNA_WT",                 "Qs.RNA_WT",
                        "log2FC.delCHD1_vs_WT",   "Qs.RNA_log2FC_dCHD1vsWT",
                        "halflife",               "Qs.halflife",
                        "log2_gene_length", "Qrs.gene_length", "Qs.gene_length")]




##################################################################################################################################
##################################################################################################################################


dir.create(gsub("plots","external_data", output_plots), recursive = TRUE, showWarnings = FALSE)


pdf(paste0(gsub("plots","external_data", output_plots),"/boxplots.Qs.pdf"), width = 5, height = 4.75, useDingbats = F)

par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0), cex=0.9)


boxplot(RNA_WT ~ Qs.RNA_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2FC.CHD1_vs_mock ~ Qs.CHD1_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2FC.Pol2_vs_Input_WT ~ Qs.Pol2_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2FC.delCHD1_vs_WT ~ Qs.RNA_log2FC_dCHD1vsWT, data = ext_data, pch=20, col = "darkgrey")
boxplot(halflife ~ Qs.halflife, data = ext_data, pch=20, col = "darkgrey")

boxplot(log2_gene_length ~ Qrs.gene_length, data = ext_data, pch=20, col = "darkgrey", las=3, xlab="")
boxplot(log2_gene_length ~ Qs.RNA_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2_gene_length ~ Qs.CHD1_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2_gene_length ~ Qs.Pol2_WT, data = ext_data, pch=20, col = "darkgrey")
boxplot(log2_gene_length ~ Qs.halflife, data = ext_data, pch=20, col = "darkgrey", las=3)


dev.off()



pdf(paste0(gsub("plots","external_data", output_plots),"/corrs.pdf"), width = 6, height = 6, useDingbats = F)

par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0), cex=0.9)

pheatmap(cor(ext_data[, !grepl("Qs|Qrs|name|commonName", colnames(ext_data))], method = "spearman"),
         main = "Spearman´s R", 
         cellwidth = 16, cellheight = 16,
         col = colorRampPalette(rev(brewer.pal(n = 9, "RdBu")))(101),
         breaks = seq(-1,1, length.out=101))

dev.off()



##################################################################################################################################
##################################################################################################################################


if(all(identical(names(scaled_log2_gene_length), rownames(log2_norm_counts)),
       identical(SampleTable$LibraryID,          colnames(log2_norm_counts)),
       identical(SampleTable$Time[    SampleTable$Fraction ==  "total RNA"], SampleTable$Time[    SampleTable$Fraction ==  "pull out fraction"]),
       identical(SampleTable$Genotype[SampleTable$Fraction ==  "total RNA"], SampleTable$Genotype[SampleTable$Fraction ==  "pull out fraction"]),
       identical(SampleTable$Dataset[ SampleTable$Fraction ==  "total RNA"], SampleTable$Dataset[ SampleTable$Fraction ==  "pull out fraction"]))){
        
        
        ################################################################################
        
        log2_norm_pullout <- log2_norm_counts[, SampleTable$LibraryID[SampleTable$Fraction == "pull out fraction"]]
        log2_norm_total   <- log2_norm_counts[, SampleTable$LibraryID[SampleTable$Fraction ==   "total RNA"]]
        
        write.table(log2_norm_pullout, file = paste0(output_tables,"/log2_norm_pullout.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
        write.table(log2_norm_total,   file = paste0(output_tables,"/log2_norm_total.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
        
        ################################################################################
        
        log2_ip_over_total <- log2_norm_pullout - log2_norm_total
        colnames(log2_ip_over_total) <- paste(colnames(log2_norm_pullout), colnames(log2_norm_total), sep="-")
        
        write.table(log2_ip_over_total, file = paste0(output_tables,"/log2_ip_over_total.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
        
        ################################################################################
        
        log2_lennorm_counts <- log2_norm_counts - scaled_log2_gene_length
        
        log2_lennorm_pullout <- log2_lennorm_counts[, SampleTable$LibraryID[SampleTable$Fraction == "pull out fraction"]]
        log2_lennorm_total   <- log2_lennorm_counts[, SampleTable$LibraryID[SampleTable$Fraction ==   "total RNA"]]
        
        write.table(log2_lennorm_pullout, file = paste0(output_tables,"/log2_lennorm_pullout.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
        write.table(log2_lennorm_total,   file = paste0(output_tables,"/log2_lennorm_total.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
        
        ################################################################################
}    




##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################


gene_class <- "protein_coding"


dir.create(paste0(output_plots, "/Boxplots/"), recursive = TRUE, showWarnings = FALSE)



for(my_plot_type in c("ip_over_total","norm_pullout","norm_total","lennorm_pullout","lennorm_total")){
        
        
        print(my_plot_type)
        
        my_plot_data <- get(paste0("log2_", my_plot_type))
        
        my_order <- c(which(gsub("\\-.*","", colnames(my_plot_data)) %in% SampleTable$LibraryID[SampleTable$Genotype == unique(SampleTable$Genotype)[1]]),
                      which(gsub("\\-.*","", colnames(my_plot_data)) %in% SampleTable$LibraryID[SampleTable$Genotype == unique(SampleTable$Genotype)[2]]))
        
        my_subset <- rownames(my_plot_data) %in% my_genes_Scer$gene_id[my_genes_Scer$gene_biotype == gene_class] 
        
        my_plot_data <- my_plot_data[my_subset, my_order]
        
        ################################################################################
        
        
        pdf(paste0(output_plots,"/Boxplots/Boxplot_reps.",my_plot_type,".pdf"), width = 8, height = 4.75, useDingbats = F)
        
        
        par(mfrow=c(1,1), mar = c(5,5,2,2), oma = c(3,3,3,3), mgp = c(2,1,0), cex=0.9)
        
        my_ylab <- ifelse(my_plot_type == "ip_over_total", "log2(IP/Total)", paste0("log2 ", my_plot_type))
        
        
        my_bp_cond <- factor(SampleTable$Genotype[match(gsub("\\-.*","", colnames(my_plot_data)), SampleTable$LibraryID)],
                             levels = unique(SampleTable$Genotype))
        
        bp <- boxplot(my_plot_data, 
                      main = paste0(gene_class, " (n=", sum(my_subset),")"),
                      col = my_color_palette[my_bp_cond], 
                      las=2, xaxt = "n",
                      ylab = my_ylab, outline=F)
        
        axis(side = 1, at = 1:ncol(my_plot_data), labels = colnames(my_plot_data), las=2, cex.axis=0.6)
        
        abline(v =  0.5+seq(0, ncol(my_plot_data), 
                            ncol(my_plot_data)/length(unique(SampleTable$Genotype))/length(unique(SampleTable$Dataset))), 
               col="grey")  
        
        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
        
        plot.new()
        
        legend("bottom", 
               legend = levels(my_bp_cond), horiz = TRUE, 
               bty = "n", cex = 0.8, pt.cex = 1, 
               fill = my_color_palette)
        
        
        dev.off()
        
}




##################################################################################################################################
##################################################################################################################################


dir.create(paste0(output_plots, "/curves"))





minmax_scale <- function(x){ ((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE)-min(x, na.rm=TRUE))) }


my_fav_genes <- c("ACT1","PGK1","PMA1","RAD51","RDN25-1")

xn <- seq(0,45,0.01)


isLenNorm <- "noLenNorm"
my_plot_type <- "norm_pullout"
islog2 <- "cnt"
my_fav_gene <- my_fav_genes[1]
my_formula <- "logistic"



for(my_plot_type in  c("ip_over_total","norm_pullout","norm_total","lennorm_pullout","lennorm_total")){
        
        
        print(my_plot_type)
        
        my_plot_data <- get(paste0("log2_", my_plot_type))
        
        my_ylab <- ifelse(my_plot_type == "ip_over_total", "ratio log2(IP/Total)", paste0("log2 counts ", my_plot_type))
        my_ylab <- ifelse(islog2 == "log2", my_ylab, gsub("log2","",my_ylab))
        
        
        ################################################################################        
        
        pdf(paste0(output_plots, "/curves/fit.",islog2,".",my_plot_type,".",".pdf"), width = 9, height = 5)
        
        
        for(my_fav_gene in my_fav_genes){
                
                
                par(mfrow=c(1,2), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))
                
                match_samples <- match(gsub("\\-.*","", colnames(my_plot_data)),SampleTable$LibraryID)
                
                
                df <- data.frame(xall =      SampleTable$Time[match_samples],
                                 genotype =  SampleTable$Genotype[match_samples], 
                                 dataset =   SampleTable$Dataset[match_samples])
                
                
                df$genotype <- factor(df$genotype, levels = unique(df$genotype))
                df$dataset <- factor(df$dataset, levels = unique(df$dataset))
                
                
                if(grepl("RDN", my_fav_gene)){
                        
                        df$yall <- as.numeric(my_plot_data[rownames(my_plot_data) ==  my_fav_gene,])
                        
                } else {
                        df$yall <- as.numeric(my_plot_data[rownames(my_plot_data) == my_genes$gene_id[my_genes$gene_name %in% my_fav_gene],])
                }
                
                if(islog2 == "cnt"){
                        df$yall <- 2^df$yall
                }
                
                
                
                for(my_formula in c("logistic","gompertz")){
                        
                        
                        plot(df$xall, df$yall, 
                             main = my_formula, xlab = "Time", ylab = my_ylab,
                             col = my_color_palette[df$genotype], 
                             pch = (15:17)[df$dataset], cex = 0.75)
                        
                        mtext(text = my_fav_gene, side = 3, outer = TRUE, font = 2, cex = 1.5)
                        
                        ################################################################################        
                        
                        for(my_fitgroup in levels(df$genotype)){
                                
                                
                                x <- df$xall[df$genotype == my_fitgroup]
                                y <- df$yall[df$genotype == my_fitgroup]
                                
                                
                                if(my_formula == "logistic"){
                                        
                                        ################################################################################        
                                        
                                        ce <- tryCatch(expr = {
                                                
                                                fit <- nls(y ~ asym/(1+exp(-rate*(x-xmid))),
                                                           start = c(xmid = 20, rate = 0.1, asym = max(y)),
                                                           control = nls.control(maxiter = 200))
                                                
                                                asym <- as.numeric(coef((fit))["asym"])
                                                xmid <- as.numeric(coef((fit))["xmid"])
                                                rate <- as.numeric(coef((fit))["rate"])
                                                
                                                aic <- AIC(fit)
                                                
                                                c(asym, xmid, rate, aic)
                                                
                                        }, error = function(c) c(0,0,0,NA))
                                        
                                        
                                        asym <- ce[1]
                                        xmid <- ce[2] 
                                        rate <- ce[3]
                                        
                                        
                                        if(rate != 0){
                                                
                                                yh <- asym/(1+exp(-rate*(xn-xmid)))
                                                
                                                alow <- asym/(1+exp(-rate*(0-xmid)))
                                                
                                                abline(h = alow, lty = 2,
                                                       col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                slope <- asym * rate / 4
                                                inter <- (asym/2) - slope * xmid
                                                
                                                abline(coef = c(inter, slope), lty=2,
                                                       col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                xlag <- (alow - inter) / slope
                                                
                                                # points(xlag, alow, pch=8,
                                                #        col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                
                                                lines(xn, yh, lwd = 1.5,
                                                      col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                text(x = 1, 
                                                     y = ifelse(my_fitgroup == levels(df$genotype)[2], max(df$yall)*0.99, max(df$yall)*0.92),
                                                     col = ifelse(my_fitgroup == levels(df$genotype)[2],"#999999","#E69F00"),
                                                     adj = 0, label = paste("rate =", round(rate, 4)))
                                                
                                                text(x = 1, 
                                                     y = ifelse(my_fitgroup == levels(df$genotype)[2], max(df$yall)*0.79, max(df$yall)*0.72),
                                                     col = ifelse(my_fitgroup == levels(df$genotype)[2],"#999999","#E69F00"),
                                                     adj = 0, label = paste("lag =", round(xlag, 2)))
                                        } else {
                                                
                                                fit <-  loess(y ~ x, span = 3)
                                                
                                                lines(xn, predict(fit, newdata = data.frame(x = xn)), lwd = 1.5,
                                                      col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                text(x = 1, 
                                                     y = ifelse(my_fitgroup == levels(df$genotype)[2], max(df$yall)*0.99, max(df$yall)*0.92),
                                                     col = ifelse(my_fitgroup == levels(df$genotype)[2],"#999999","#E69F00"),
                                                     adj = 0, label = "loess")
                                        }
                                        
                                        mtext(text = expression(frac(asym,(1+e^{-rate*(t-t[mid])}))), adj = 0, side = 1, line = 4, cex=0.75)
                                        
                                        rm(list = "ce")
                                        
                                        ################################################################################        
                                        
                                } else {
                                        
                                        ################################################################################        
                                        
                                        ce <- tryCatch(expr = {
                                                
                                                # fit <- nls(y ~ a*exp(-exp(b-g*x)),
                                                #            start = c(a = max(y), b = 1, g = 0.1),
                                                #            control = nls.control(maxiter = 200))
                                                
                                                fit <- nls(y ~ a*exp(-b*exp(-g*x)),
                                                           start = c(a = max(y), b = 1, g = 0.1),
                                                           control = nls.control(maxiter = 200))
                                                
                                                a <- as.numeric(coef((fit))["a"])
                                                b <- as.numeric(coef((fit))["b"])
                                                g <- as.numeric(coef((fit))["g"])
                                                #d <- as.numeric(coef((fit))["d"])
                                                
                                                aic <- AIC(fit)
                                                
                                                c(a, b, g, aic)
                                                
                                        }, error = function(c) c(0,0,0,NA))
                                        
                                        a <- ce[1]
                                        b <- ce[2]
                                        g <- ce[3]
                                        #d <- ce[4]
                                        
                                        if(a != 0){
                                                
                                                yh <- a*exp(-b*exp(-g*xn))
                                                
                                                lines(xn, yh, lwd = 1.5,
                                                      col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))  
                                                
                                                text(x = 1, 
                                                     y = ifelse(my_fitgroup == levels(df$genotype)[2], max(df$yall)*0.99, max(df$yall)*0.92),
                                                     col = ifelse(my_fitgroup == levels(df$genotype)[2],"#999999","#E69F00"),
                                                     adj = 0, label = paste("g =", round(g, 4)))
                                                
                                        } else {
                                                
                                                fit <-  loess(y ~ x, span = 3)
                                                
                                                lines(xn, predict(fit, newdata = data.frame(x = xn)), lwd = 1.5,
                                                      col = unique(my_color_palette[df$genotype][df$genotype == my_fitgroup]))
                                                
                                                text(x = 1, 
                                                     y = ifelse(my_fitgroup == levels(df$genotype)[2], max(df$yall)*0.99, max(df$yall)*0.92),
                                                     col = ifelse(my_fitgroup == levels(df$genotype)[2],"#999999","#E69F00"),
                                                     adj = 0, label = "loess")
                                        }
                                        
                                        mtext(text = expression( a*e^{-b*e^{-g*t}}), adj = 0, side = 1, line = 4, cex=0.75)
                                        
                                        rm(list = "ce")
                                        
                                        ################################################################################        
                                }
                                
                                ################################################################################        
                        }
                }
                
                rm(list = "df")
        }  
        
        dev.off()
        
        ################################################################################        
}



##################################################################################################################################
##################################################################################################################################







##################################################################################################################################
##################################################################################################################################


################################################################################        


fitLogistic <- function(x,y){
        
        ce <- tryCatch(expr = {
                
                fit <- nls(y ~ asym/(1+exp(-rate*(x-xmid))),
                           start = c(xmid = 20, rate = 0.1, asym = max(y)),
                           control = nls.control(maxiter = 200))
                
                asym <- as.numeric(coef((fit))["asym"])
                xmid <- as.numeric(coef((fit))["xmid"])
                rate <- as.numeric(coef((fit))["rate"])
                
                aic <- AIC(fit)
                
                c(asym, xmid, rate, aic)
                
        }, error = function(c) c(NA,NA,NA,NA))
        
        return(ce)
}


################################################################################        



calcLag <- function(asym, xmid, rate){
        
        alow <- asym/(1+exp(-rate*(0-xmid)))
        
        slope <- asym * rate / 4
        inter <- (asym/2) - slope * xmid
        
        xlag <- (alow - inter) / slope
        
        return(xlag)
}


################################################################################


fitGompertz <- function(x,y){
        
        ce <- tryCatch(expr = {
                
                # fit <- nls(y ~ a*exp(-exp(b-g*x)),
                #            start = c(a = max(y), b = 1, g = 0.1),
                #            control = nls.control(maxiter = 200))
                
                fit <- nls(y ~ a*exp(-b*exp(-g*x)),
                           start = c(a = max(y), b = 1, g = 0.1),
                           control = nls.control(maxiter = 200))
                
                a <- as.numeric(coef((fit))["a"])
                b <- as.numeric(coef((fit))["b"])
                g <- as.numeric(coef((fit))["g"])
                #d <- as.numeric(coef((fit))["d"])
                
                aic <- AIC(fit)
                
                c(a, b, g, aic)
                
        }, error = function(c) c(NA,NA,NA,NA))
        
        return(ce)
        
}




################################################################################        


dir.create(paste0(output_plots, "/fitparams"), recursive = TRUE, showWarnings = FALSE)


my_gene_id <- "YFL039C"


islog2 <- "cnt"

my_formula <- "logistic"

my_gene_ids <-  c(my_genes$gene_id[my_genes$gene_name %in% my_fav_genes], "RDN25-1")
my_gene_ids <- rownames(log2_ip_over_total)




for(my_plot_type in  c("ip_over_total","norm_pullout","lennorm_pullout")){
        
        
        print(my_plot_type)
        
        my_plot_data <- get(paste0("log2_", my_plot_type))
        
        ################################################################################        
        
        for(my_formula in c("logistic","gompertz")){
                
                
                my_fit_data <- data.frame(gene_id = my_gene_ids,
                                          gene_name = my_genes$gene_name[match(my_gene_ids, my_genes$gene_id)],
                                          matrix(nrow = length(my_gene_ids), ncol = 8))
                
                
                for(i in 1:nrow(my_fit_data)){
                        
                        my_gene_id <- my_fit_data$gene_id[i]
                        
                        match_samples <- match(gsub("\\-.*","", colnames(my_plot_data)),SampleTable$LibraryID)
                        
                        
                        df <- data.frame(xall =      SampleTable$Time[match_samples],
                                         genotype =  SampleTable$Genotype[match_samples], 
                                         dataset =   SampleTable$Dataset[match_samples])
                        
                        
                        df$genotype <- factor(df$genotype, levels = unique(df$genotype))
                        df$dataset <- factor(df$dataset, levels = unique(df$dataset))
                        
                        df$yall <- as.numeric(my_plot_data[rownames(my_plot_data) ==  my_gene_id,])
                        
                        df$yall <- 2^df$yall
                        
                        
                        if(my_formula == "logistic"){
                                
                                params1 <- fitLogistic(df$xall[df$genotype == levels(df$genotype)[1]], df$yall[df$genotype == levels(df$genotype)[1]])
                                params1 <- c(params1, calcLag(params1[1], params1[2], params1[3]))
                                
                                params2 <- fitLogistic(df$xall[df$genotype == levels(df$genotype)[2]], df$yall[df$genotype == levels(df$genotype)[2]])
                                params2 <- c(params2, calcLag(params2[1], params2[2], params2[3]))
                                
                                my_fit_data[i,3:12] <- c(params1, params2)
                                
                                colnames(my_fit_data)[3:12] <- c(paste(c("asym", "xmid", "rate", "aic", "lag"), levels(df$genotype)[1], sep = "_"), 
                                                                 paste(c("asym", "xmid", "rate", "aic", "lag"), levels(df$genotype)[2], sep = "_"))
                                
                        } else {
                                
                                my_fits <- c(fitGompertz(df$xall[df$genotype == levels(df$genotype)[1]], df$yall[df$genotype == levels(df$genotype)[1]]),
                                             fitGompertz(df$xall[df$genotype == levels(df$genotype)[2]], df$yall[df$genotype == levels(df$genotype)[2]]))
                                
                                my_fit_data[i,3:10] <- my_fits
                                
                                colnames(my_fit_data)[3:10] <- c(paste(c("a", "b", "g", "aic"), levels(df$genotype)[1], sep = "_"), 
                                                                 paste(c("a", "b", "g", "aic"), levels(df$genotype)[2], sep = "_"))
                        }
                        
                        rm(list = "df")
                        
                }
                
                write.table(my_fit_data, file = paste0(output_tables,"/fit.",islog2,".",my_plot_type,".",my_formula,".txt"), 
                            quote = F, sep = "\t", row.names = FALSE)
                
                
                if(my_formula == "logistic"){
                        
                        pdf(file = paste0(output_plots,"/fitparams/smooth.",islog2,".",my_plot_type,".",my_formula,".pdf"),
                            width = 5, height = 5, useDingbats = FALSE)
                        
                        par(mfrow=c(1, 1), mar=c(4,4,4,4), mgp = c(2,1,0), oma=c(2,2,2,2))
                        
                        
                        smoothScatter(my_fit_data[,c(5,10)], xlim=c(0,0.8), ylim=c(0,0.8),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        abline(coef = c(0,1))
                        
                        smoothScatter(my_fit_data[,c(3,8)], xlim=c(0,250), ylim=c(0,250),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        abline(coef = c(0,1))
                        
                        smoothScatter(my_fit_data[,c(4,9)], xlim=c(0,30), ylim=c(0,30),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        abline(coef = c(0,1))
                        
                        smoothScatter(my_fit_data[,c(7,12)], xlim=c(0,15), ylim=c(0,15),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        abline(coef = c(0,1))
                        
                        dev.off()
                        
                        ################################################################################        
                        
                        pdf(file = paste0(output_plots,"/fitparams/comparison.",islog2,".",my_plot_type,".",my_formula,".pdf"),
                            width = 5, height = 5, useDingbats = FALSE)
                        
                        par(mfrow=c(1, 1), mar=c(4,4,4,4), mgp = c(2,1,0), oma=c(2,2,2,2))
                        
                        
                        my_merged_data <- merge(my_fit_data, ext_data, by.x = "gene_id", by.y = "name")
                        
                        smoothScatter(my_merged_data[,c(22,10)], xlim=c(0,20), ylim=c(0,0.8),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        smoothScatter(my_merged_data[,c(22,5)], xlim=c(0,20), ylim=c(0,0.8),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        smoothScatter(my_merged_data[,c(22,9)], xlim=c(0,20), ylim=c(0,30),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        smoothScatter(my_merged_data[,c(22,4)], xlim=c(0,20), ylim=c(0,30),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        smoothScatter(my_merged_data[,c(22,8)], xlim=c(0,20), ylim=c(0, 250),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                        
                        smoothScatter(my_merged_data[,c(22,3)], xlim=c(0,20), ylim=c(0,250),
                                      colramp = colorRampPalette(c("white", blues9, "black")),
                                      cex = 0, useRaster=TRUE)
                       
                        dev.off()
                        
                        ################################################################################        
                }
                
                if(my_formula == "logistic" & my_plot_type == "norm_pullout" & islog2 == "cnt"){
                        
                        
                        ext_data <- merge(ext_data, 
                                          data.frame(lag_control = my_fit_data[,12],
                                                     name = my_fit_data[,1]), by = "name")
                        
                        ext_data$Qs.lag_control <-  c("q1","q2","q3","q4","q5")[cut(x = ext_data$lag_control,
                                                                                    breaks = quantile(ext_data$lag_control, probs = seq(0, 1, 0.20), na.rm=TRUE),
                                                                                    include.lowest = TRUE)]
                        
                        write.table(ext_data, file = paste0(gsub("plots","external_data", output_plots),"/ext_data.txt"),
                                    quote = F, sep = "\t", row.names = FALSE)
                        
                }
                
        }  
        
        ################################################################################        
}



##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################


islog2 <- "cnt"


for(my_plot_type in c("ip_over_total","norm_pullout","norm_total","lennorm_pullout","lennorm_total")){
        
        
        print(my_plot_type)
        
        my_plot_data <- get(paste0("log2_", my_plot_type))
        
        dir.create(paste0(output_plots, "/Dotplots/",my_plot_type), recursive = TRUE, showWarnings = FALSE)
        
        
        my_ylab <- ifelse(my_plot_type == "ip_over_total", "ratio log2(IP/Total)", paste0("log2 counts ", my_plot_type))
        my_ylab <- ifelse(islog2 == "log2", my_ylab, gsub("log2","",my_ylab))
        
        
        for(ext_assay in c("Pol2_WT","RNA_WT","CHD1_WT","RNA_log2FC_dCHD1vsWT","gene_length","halflife","lag_control")){
                
                print(ext_assay)
                
                my_df <- data.frame()
                
                ################################################################################        
                
                for(qs in c("q1","q2","q3","q4","q5")){
                        
                        print(qs)
                        
                        my_subset <- rownames(my_plot_data) %in% ext_data$name[ext_data[,paste0("Qs.", ext_assay)] == qs]
                        
                        if(islog2 == "log2"){
                                my_Means <- colMeans(my_plot_data[my_subset,])
                        } else {
                                my_Means <- colMeans(2^my_plot_data[my_subset,])
                        }
                        
                        my_lib_ids <- SampleTable$LibraryID %in% gsub("-.*","",colnames(my_plot_data))
                        
                        stopifnot(identical(SampleTable$LibraryID[my_lib_ids], gsub("\\-.*","",colnames(my_plot_data))))
                        
                        my_df <- rbind(my_df,
                                       data.frame(Genotype = factor(SampleTable$Genotype[my_lib_ids], levels = unique(SampleTable$Genotype[my_lib_ids])),
                                                  Time = SampleTable$Time[my_lib_ids],
                                                  Mean = my_Means,
                                                  Dataset =  SampleTable$Dataset[my_lib_ids],
                                                  Type = my_ylab,
                                                  Q = qs))
                }
                
                ################################################################################        
                
                if(my_plot_type %in% c("ip_over_total","norm_pullout","lennorm_pullout")){
                        
                        
                        for(my_formula in c("logistic","gompertz")){
                                
                                
                                print(my_formula)
                                
                                if(my_formula == "logistic"){
                                        
                                        ################################################################################        
                                        
                                        ggp <- ggplot(my_df, 
                                                      aes(x = Time, y = Mean, col = Genotype, group = Type, shape = Dataset)) + 
                                                theme_bw() +
                                                #ggtitle(label = ext_assay) +
                                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                                      #axis.text.x = element_blank(), 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      strip.background = element_blank(),
                                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                                facet_wrap(~ Q, nrow = 1) +
                                                geom_smooth(data = my_df, 
                                                            method = "nls", formula = y ~  asym/(1+exp(-rate*(x-xmid))),
                                                            method.args = list(start = c(asym = max(my_df$Mean), xmid = 20, rate = 0.1), 
                                                                               control = nls.control(maxiter = 200)), 
                                                            se = FALSE,
                                                            aes(group = Genotype, fill = Genotype), alpha=0.2) +
                                                xlab(label = "Time") +
                                                geom_point(size = 1) +
                                                # coord_cartesian(ylim = c(-6,8)) +
                                                # scale_y_continuous(breaks=seq(-6,8,2)) +
                                                scale_color_manual(values = my_color_palette) +
                                                scale_fill_manual(values = my_color_palette) 
                                        
                                        
                                        ggsave(filename = paste0(output_plots,"/Dotplots/",my_plot_type,"/Dotplot.",islog2,".",my_formula,".",ext_assay,".pdf"),
                                               plot =  ggp,
                                               width = 12, height = 6)
                                        
                                        
                                        
                                } else {
                                        
                                        ################################################################################        
                                        
                                        ggp <- ggplot(my_df, 
                                                      aes(x = Time, y = Mean, col = Genotype, group = Type, shape = Dataset)) + 
                                                theme_bw() +
                                                #ggtitle(label = ext_assay) +
                                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                                      #axis.text.x = element_blank(), 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      strip.background = element_blank(),
                                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                                facet_wrap(~ Q, nrow = 1) +
                                                geom_smooth(data = my_df, 
                                                            method = "nls", formula = y ~ a*exp(-b*exp(-g*x)),
                                                            method.args = list(start = c(a =  max(my_df$Mean), b = 1, g = 0.1), 
                                                                               control = nls.control(maxiter = 200)),
                                                            se = FALSE,
                                                            aes(group = Genotype, fill = Genotype), alpha=0.2) +
                                                xlab(label = "Time") +
                                                geom_point(size = 1) +
                                                # coord_cartesian(ylim = c(-6,8)) +
                                                # scale_y_continuous(breaks=seq(-6,8,2)) +
                                                scale_color_manual(values = my_color_palette) +
                                                scale_fill_manual(values = my_color_palette) 
                                        
                                        
                                        ggsave(filename = paste0(output_plots,"/Dotplots/",my_plot_type,"/Dotplot.",islog2,".",my_formula,".",ext_assay,".pdf"),
                                               plot =  ggp,
                                               width = 12, height = 6) 
                                        
                                        ################################################################################        
                                }    
                        }
                        
                } else {
                        
                        ################################################################################        
                        
                        ggp <- ggplot(my_df, 
                                      aes(x = Time, y = Mean, col = Genotype, group = Type, shape = Dataset)) + 
                                theme_bw() +
                                #ggtitle(label = ext_assay) +
                                theme(text = element_text(size = 12), aspect.ratio = 1, 
                                      #axis.text.x = element_blank(), 
                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(),
                                      plot.title = element_text(hjust = 0.5, face = "bold")) +
                                facet_wrap(~ Q, nrow = 1) +
                                geom_smooth(data = my_df, 
                                            method = "loess", span = 3, se = FALSE,
                                            aes(group = Genotype, fill = Genotype), alpha=0.2) +
                                xlab(label = "Time") +
                                geom_point(size = 1) +
                                # coord_cartesian(ylim = c(-6,8)) +
                                # scale_y_continuous(breaks=seq(-6,8,2)) +
                                scale_color_manual(values = my_color_palette) +
                                scale_fill_manual(values = my_color_palette) 
                        
                        
                        ggsave(filename = paste0(output_plots,"/Dotplots/",my_plot_type,"/Dotplot.",islog2,".loess.",ext_assay,".pdf"),
                               plot =  ggp,
                               width = 12, height = 6) 
                        
                        ################################################################################        
                }
        }
}



##################################################################################################################################
##################################################################################################################################










##################################################################################################################################
##################################################################################################################################






writeLines(capture.output(sessionInfo()), output_sessionInfo)






##################################################################################################################################
##################################################################################################################################

