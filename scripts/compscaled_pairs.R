






##################################################################################################################################
##################################################################################################################################




library(ShortRead)
library(rtracklayer)

library(RColorBrewer)



##################################################################################################################################
##################################################################################################################################




args = commandArgs(trailingOnly=TRUE)



input_table <- grep("SampleTable", args, value = TRUE)

input_anno <- grep("gtf", args, value = TRUE)

input_lnc <- grep("log2_norm_counts", args, value = TRUE)


output_file <- grep("composite.*.pdf", args, value = TRUE)

matrix_files <- grep("_GB.*matrix.rds", args, value = TRUE)
matrix_files <- matrix_files[order(matrix_files)]


stopifnot(length(matrix_files) == 2)


##################################################################################################################################
##################################################################################################################################


minmax_scale <- function(x, idx = c(1, length(x))){ 
        
        ((x - min(x[idx[1]:idx[2]], na.rm=TRUE)) / (max(x[idx[1]:idx[2]], na.rm=TRUE)-min(x[idx[1]:idx[2]], na.rm=TRUE))) 
}


##################################################################################################################################
##################################################################################################################################





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

SampleTable$Samples <- paste(SampleTable$Conditions, SampleTable$Time, sep = "_")


##################################################################################################################################
##################################################################################################################################


high_extremes <- c("YGR190C","YGR191W","YGR192C","YGR193C","YGR194C","RDN25-1","YLR154W-A","YLR154W-B","YLR154W-C","RDN58-1",  
                   "YLR154W-E","YLR154W-F","RDN18-1","RDN5-1","RDN25-2","RDN58-2","RDN18-2","RDN5-2","YLR154C-H","YLR155C", 
                   "YLR160C","YLR161W","RDN5-6", "YLR162W","YLR162W-A","YLR163C","YLR163W-A","YPR078C","YPR079W","YPR080W",  
                   "YPR081C")


transposons <- c("YAR009C", "YAR010C", "YBL100W-B", "YBL100W-B", "YBL100W-A", "YBL016W", "YBL005W-B", "YBL005W-B", "YBL005W-A", "YBR012W-B", "YBR012W-B", "YBR012W-A",
                 "YCL019W", "YCL019W", "YCL020W", "YDR034C-D", "YDR034C-D", "YDR034C-C", "YDR098C-B", "YDR098C-B", "YDR098C-A", "YDR170W-A", "YDR210W-B", "YDR210W-B", "YDR210W-A",
                 "YDR210C-D", "YDR210C-D", "YDR210C-C", "YDR261W-B", "YDR261W-B", "YDR261W-A", "YDR261C-D", "YDR261C-D", "YDR261C-C", "YDR289C", "YDR316W-B", "YDR316W-B", "YDR316W-A", 
                 "YDR365W-B", "YDR365W-B", "YDR365W-A", "YER104W", "YER138C", "YER138C", "YER137C-A", "YER160C", "YER160C", "YER159C-A", "YFL002W-A", "YFL002W-A", "YFL002W-B", "YGR027W-B", 
                 "YGR027W-B", "YGR027W-A", "YGR038C-B", "YGR038C-B", "YGR038C-A", "YGR102C", "YGR109W-B", "YGR109W-B", "YGR109W-A", "YGR161W-B", "YGR161W-B", "YGR161W-A", "YGR161C-D", "YGR161C-D",
                 "YGR161C-C", "YHL009W-B", "YHL009W-B", "YHL009W-A", "YHR154W", "YHR214C-B", "YHR214C-B", "YHR214C-C", "YIL082W-A", "YIL082W-A", "YJL114W", "YJL113W", "YJL113W", "YJL047C", "YJR027W", 
                 "YJR027W", "YJR026W", "YJR029W", "YJR029W", "YJR028W", "YLL002W", "YLR035C-A", "YLR157C-B", "YLR157C-B", "YLR157C-A", "YLR227W-B", "YLR227W-B", "YLR227W-A", "YLR256W-A", "YLR410W-B", 
                 "YLR410W-B", "YLR410W-A", "YLR411W", "YML045W", "YML045W", "YML045W-A", "YML039W", "YML039W", "YML040W", "YMR045C", "YMR045C", "YMR046C", "YMR050C", "YMR050C", "YMR051C", "YNL284C-B", 
                 "YNL284C-B", "YNL284C-A", "YNL054W-A", "YOL103W-B", "YOL103W-B", "YOL103W-A", "YOR142W-B", "YOR142W-B", "YOR142W-A", "YOR144C", "YOR192C-B", "YOR192C-B", "YOR192C-A", "YOR343W-B", 
                 "YOR343W-B", "YOR343W-A", "YPL257W-B", "YPL257W-B", "YPL257W-A", "YPR137C-B", "YPR137C-B", "YPR137C-A", "YPR158W-B", "YPR158W-B", "YPR158W-A", "YPR158C-D", "YPR158C-D", "YPR158C-C", "YPR164W")


##################################################################################################################################
##################################################################################################################################








##################################################################################################################################
##################################################################################################################################

######################################################     Read Mats    ########################################################## 


log2_norm_counts <- read.table(input_lnc)



my_sample_id1 <- gsub(".*\\/|_GB.*", "", matrix_files[1])
my_sample_id2 <- gsub(".*\\/|_GB.*", "", matrix_files[2])

mat_1 <- readRDS(matrix_files[1])
mat_2 <- readRDS(matrix_files[2])


mat_1 <- mat_1[(rownames(mat_1) %in% rownames(log2_norm_counts)),]
mat_2 <- mat_2[(rownames(mat_2) %in% rownames(log2_norm_counts)),]



##################################################################################################################################
##################################################################################################################################



log2_gene_length <- log2(width(my_genes))
names(log2_gene_length) <- my_genes$gene_id




##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

######################################################      Plotting     ######################################################### 


my_condition1 <- which(unique(SampleTable$Genotype) == SampleTable$Genotype[SampleTable$LibraryID == my_sample_id1])
my_condition2 <- which(unique(SampleTable$Genotype) == SampleTable$Genotype[SampleTable$LibraryID == my_sample_id2])

my_color1 <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition1]
my_color2 <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition2]



my_grouping_names <- c("Qs.gene_length")


ylim_table <- read.csv("../metaplot_ranges.csv", header = TRUE)
ylim_table <- ylim_table[ylim_table$experiment == "chd1vswt",]

ni=1


################################################################################


pdf(file = output_file, width = 7.5, height = 4.5)

for(to_remove in c("high_extremes+transposons","high_extremes")){
        
        
        ################################################################################
        
        if(to_remove == "high_extremes"){
                
                matf_1 <- mat_1[!(rownames(mat_1) %in% high_extremes),]
                matf_2 <- mat_2[!(rownames(mat_2) %in% high_extremes),]
                
        } else if(to_remove == "high_extremes+transposons"){
                
                matf_1 <- mat_1[!(rownames(mat_1) %in% c(high_extremes, transposons)),]
                matf_2 <- mat_2[!(rownames(mat_2) %in% c(high_extremes, transposons)),]
        }
        
        ################################################################################
        
        log2_gene_length_filt <- log2_gene_length[match(rownames(matf_1), names(log2_gene_length))]
        
        stopifnot(identical(rownames(matf_1), names(log2_gene_length_filt)))
        
        Qs_gene_length <- c("q1","q2","q3","q4","q5")[cut(x = log2_gene_length_filt,
                                                          breaks = quantile(log2_gene_length_filt, probs = seq(0, 1, 0.20)),
                                                          include.lowest = TRUE)]
        
        
        ext_data <-  data.frame(log2_gene_length = log2_gene_length_filt,
                                Qs.gene_length   = Qs_gene_length,
                                name = names(log2_gene_length_filt), 
                                stringsAsFactors = FALSE)
        
        
        stopifnot(identical(rownames(matf_1), ext_data$name))
        
        ################################################################################
        
        for(ni in seq_along(my_grouping_names)){
                
                dfends <- data.frame(matrix(nrow = 0, ncol = 4))
                colnames(dfends) <- c("LibraryID","Qs","TSS","TTS")
                
                ################################################################################
                
                for(gi in  c("q1","q2","q3","q4","q5")){
                        
                        which_time <- ylim_table$time.of.label == SampleTable$Time[SampleTable$LibraryID == my_sample_id1]
                        which_frac <- ylim_table$Fraction == SampleTable$Fraction[SampleTable$LibraryID == my_sample_id1]
                        which_q <- ylim_table$gene.quintile == gi
                        
                        ylim <- ylim_table$max_mean.read.coverage[which_time & which_frac & which_q]
                        
                        ydata1 <- colMeans(matf_1[rownames(matf_1) %in% ext_data$name[ext_data[,my_grouping_names[ni]] == gi],])
                        ydata2 <- colMeans(matf_2[rownames(matf_2) %in% ext_data$name[ext_data[,my_grouping_names[ni]] == gi],])
                        
                        
                        dfends <- rbind(dfends, 
                                        data.frame(LibraryID = my_sample_id1, 
                                                   Qs = gi, 
                                                   TSS = mean(ydata1[1100:1300], na.rm=TRUE),
                                                   Mid = mean(ydata1[1900:2100], na.rm=TRUE),
                                                   TTS = mean(ydata1[2700:2900], na.rm=TRUE)),
                                        data.frame(LibraryID = my_sample_id2, 
                                                   Qs = gi, 
                                                   TSS = mean(ydata2[1100:1300], na.rm=TRUE),
                                                   Mid = mean(ydata2[1900:2100], na.rm=TRUE),
                                                   TTS = mean(ydata2[2700:2900], na.rm=TRUE)))
                        
                        ydata1 <- zoo::rollmean(c(rep(NA,10), ydata1, rep(NA,10)), k = 21)
                        ydata2 <- zoo::rollmean(c(rep(NA,10), ydata2, rep(NA,10)), k = 21)
                        
                        ################################################################################
                        
                        par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))   
                        
                        plot(ydata1, type = "l", ylim = c(0, ylim),
                             xlab = "", ylab = "", xaxt = "n", lwd = 2,
                             col = my_color1)
                        
                        axis(side = 1, at = c(1, ncol(matf_1)/4,ncol(matf_1)*3/4, ncol(matf_1)), labels = c("-1kb","TSS","TTS","+1kb"))
                        
                        lines(ydata2, lwd = 2,col = my_color2)
                        
                        ################################################################################
                        
                        mtext(text = paste("Mean Coverage"), side = 2, line = 3, cex = 1.2)
                        mtext(text = paste0(my_grouping_names[ni], " - ", gi), side = 3, line = 0.5, cex = 1.2)
                        
                        mtext(text = gsub(".*composite.|.pdf", "", output_file), 
                              side = 3, line = -1, cex = 1.25, font = 2, outer = TRUE)
                        
                        mtext(text = paste0("* ",to_remove, " removed"), 
                              side = 1, line = -2, cex = 0.8, adj = 1, outer = TRUE)
                        
                        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(1, 0, 0, 0), new = TRUE)
                        
                        plot.new()
                        
                        
                        legend("bottom", horiz = TRUE, lwd = 5, cex = 0.8, bty = "n",
                               legend = c(my_sample_id1, my_sample_id2),
                               col    = c(my_color1,     my_color2))
                        
                        ################################################################################
                        
                        par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))   
                        
                        plot(minmax_scale(ydata1, idx = c(750,3250)), type = "l", ylim = c(0, 1.1),
                             xlab = "", ylab = "", xaxt = "n", lwd = 2,
                             col = my_color1)
                        
                        axis(side = 1, at = c(1, ncol(matf_1)/4,ncol(matf_1)*3/4, ncol(matf_1)), labels = c("-1kb","TSS","TTS","+1kb"))
                        
                        lines(minmax_scale(ydata2, idx = c(750,3250)), lwd = 2,col = my_color2)
                        
                        ################################################################################
                        
                        mtext(text = paste("Scaled Mean Coverage"), side = 2, line = 3, cex = 1.2)
                        mtext(text = paste0(my_grouping_names[ni], " - ", gi), side = 3, line = 0.5, cex = 1.2)
                        
                        mtext(text = gsub(".*composite.|.pdf", "", output_file), 
                              side = 3, line = -1, cex = 1.25, font = 2, outer = TRUE)
                        
                        mtext(text = paste0("* ",to_remove, " removed"), 
                              side = 1, line = -2, cex = 0.8, adj = 1, outer = TRUE)
                        
                        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(1, 0, 0, 0), new = TRUE)
                        
                        plot.new()
                        
                        
                        legend("bottom", horiz = TRUE, lwd = 5, cex = 0.8, bty = "n",
                               legend = c(my_sample_id1, my_sample_id2),
                               col    = c(my_color1,     my_color2))
                        ################################################################################
                        
                }
                
                dfends$`log2(TSS/TTS)` <- log2(dfends$TSS / dfends$TTS)
                dfends$`log2(Mid/TTS)` <- log2(dfends$Mid / dfends$TTS)
                
                write.table(dfends, file = gsub("pdf",paste0(my_grouping_names[ni],".rmvd.",to_remove,".txt"),output_file), quote = FALSE, sep = "\t", row.names = FALSE)
        }
}

dev.off()


##################################################################################################################################
##################################################################################################################################


















##################################################################################################################################
##################################################################################################################################
