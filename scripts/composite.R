






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

matrix_files <- grep("_T.*matrix.rds", args, value = TRUE)
matrix_files <- matrix_files[order(matrix_files)]


stopifnot(length(matrix_files) == 2)
stopifnot(identical(gsub("TSS", "",matrix_files[1]), gsub("TTS", "",matrix_files[2])))



my_sample_id <- gsub(".*composite.|.pdf", "", output_file)

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









##################################################################################################################################
##################################################################################################################################

######################################################     Read Mats    ########################################################## 


log2_norm_counts <- read.table(input_lnc)


mat_TSS <- readRDS(grep("TSS", matrix_files, value = TRUE))
mat_TTS <- readRDS(grep("TTS", matrix_files, value = TRUE))


mat_TSS <- mat_TSS[(rownames(mat_TSS) %in% rownames(log2_norm_counts)),]
mat_TTS <- mat_TTS[(rownames(mat_TTS) %in% rownames(log2_norm_counts)),]


high_extremes <- c("YGR190C","YGR191W","YGR192C","YGR193C","YGR194C","RDN25-1","YLR154W-A","YLR154W-B","YLR154W-C","RDN58-1",  
                   "YLR154W-E","YLR154W-F","RDN18-1","RDN5-1","RDN25-2","RDN58-2","RDN18-2","RDN5-2","YLR154C-H","YLR155C", 
                   "YLR160C","YLR161W","RDN5-6", "YLR162W","YLR162W-A","YLR163C","YLR163W-A","YPR078C","YPR079W","YPR080W",  
                   "YPR081C")

mat_TSS <- mat_TSS[!(rownames(mat_TSS) %in% high_extremes),]
mat_TTS <- mat_TTS[!(rownames(mat_TTS) %in% high_extremes),]

stopifnot(identical(rownames(mat_TSS), rownames(mat_TTS)))



##################################################################################################################################
##################################################################################################################################



log2_gene_length <- log2(width(my_genes))
names(log2_gene_length) <- my_genes$gene_id


log2_gene_length <- log2_gene_length[match(rownames(mat_TSS), names(log2_gene_length))]

stopifnot(identical(rownames(mat_TSS), names(log2_gene_length)))


Qs_gene_length <- c("q1","q2","q3","q4","q5")[cut(x = log2_gene_length,
                                                  breaks = quantile(log2_gene_length, probs = seq(0, 1, 0.20)),
                                                  include.lowest = TRUE)]


ext_data <-  data.frame(log2_gene_length = log2_gene_length,
                        Qs.gene_length   = Qs_gene_length,
                        name = names(log2_gene_length), 
                        stringsAsFactors = FALSE)


stopifnot(identical(rownames(mat_TSS), ext_data$name))



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

######################################################      Plotting     ######################################################### 


my_condition <- which(unique(SampleTable$Genotype) == SampleTable$Genotype[SampleTable$LibraryID == my_sample_id])

my_color <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition]

my_color_palette <- colorRampPalette(c("white",my_color,"black"))(7)[c(-1,-7)]


my_site = "TSS"

my_grouping_names <- c("Qs.gene_length")


i=1
gi <- "q1"

################################################################################


for(i in seq_along(my_grouping_names)){
        
        
        pdf(file = output_file, width = 7.5, height = 4.5)
        
        ################################################################################
        
        for(stat in c("Mean")){
                
                par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))
                
                for(my_site in c("TSS", "TTS")){
                        
                        ################################################################################
                        
                        my_mat <- get(paste0("mat_",my_site))
                        
                        plot(colMeans(my_mat), type = "n", ylim = c(0,ifelse(stat == "Mean",5,2.5)),
                             xlab = "", ylab = "", xaxt = "n", lwd = 2,
                             col = my_color_palette[as.numeric(gsub("q","",gi))])
                        
                        axis(side = 1, at = c(1, ncol(my_mat)/2, ncol(my_mat)), labels = c("-3kb",my_site,"+3kb"))
                        
                        
                        if(my_site == "TSS"){
                                mtext(text = paste(stat,"Coverage"), side = 2, line = 3, cex = 1.25)
                        }
                        
                        ################################################################################
                        
                        for(gi in  c("q1","q2","q3","q4","q5")){
                                
                                if(stat == "Mean"){
                                        
                                        ydata <- colMeans(my_mat[rownames(my_mat) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                        
                                } else if(stat == "Median"){
                                        
                                        ydata <- colMedians(my_mat[rownames(my_mat) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                }
                                
                                ydata <- zoo::rollmean(c(rep(NA,250), ydata, rep(NA,250)), k = 501)
                                
                                lines(ydata, lwd = 2, col = my_color_palette[as.numeric(gsub("q","",gi))])
                        }
                }
                
                ################################################################################
                
                mtext(text = paste0(my_sample_id, " - ",SampleTable$Samples[SampleTable$LibraryID == my_sample_id]), 
                      side = 3, line = -1, cex = 1.25, font = 2, outer = TRUE)
                
                par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(1, 0, 0, 0), new = TRUE)
                
                plot.new()
                
                
                legend("bottom", horiz = TRUE, lwd = 5, cex = 0.8, bty = "n",
                       title = my_grouping_names[i],
                       legend =  c("q1","q2","q3","q4","q5"),
                       col = my_color_palette)
                
                ################################################################################
                
        }
        
        dev.off()
        
        
}



##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################
