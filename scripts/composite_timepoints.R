






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


stopifnot(length(matrix_files) == 7)





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



my_sample_id1 <- gsub(".*\\/|_TSS.*", "", matrix_files[1])

mat_1 <- readRDS(matrix_files[1])
mat_1 <- mat_1[(rownames(mat_1) %in% rownames(log2_norm_counts)),]


high_extremes <- c("YGR190C","YGR191W","YGR192C","YGR193C","YGR194C","RDN25-1","YLR154W-A","YLR154W-B","YLR154W-C","RDN58-1",  
                   "YLR154W-E","YLR154W-F","RDN18-1","RDN5-1","RDN25-2","RDN58-2","RDN18-2","RDN5-2","YLR154C-H","YLR155C", 
                   "YLR160C","YLR161W","RDN5-6", "YLR162W","YLR162W-A","YLR163C","YLR163W-A","YPR078C","YPR079W","YPR080W",  
                   "YPR081C")

mat_1 <- mat_1[!(rownames(mat_1) %in% high_extremes),]




##################################################################################################################################
##################################################################################################################################



log2_gene_length <- log2(width(my_genes))
names(log2_gene_length) <- my_genes$gene_id


log2_gene_length <- log2_gene_length[match(rownames(mat_1), names(log2_gene_length))]

stopifnot(identical(rownames(mat_1), names(log2_gene_length)))


Qs_gene_length <- c("q1","q2","q3","q4","q5")[cut(x = log2_gene_length,
                                                  breaks = quantile(log2_gene_length, probs = seq(0, 1, 0.20)),
                                                  include.lowest = TRUE)]


ext_data <-  data.frame(log2_gene_length = log2_gene_length,
                        Qs.gene_length   = Qs_gene_length,
                        name = names(log2_gene_length), 
                        stringsAsFactors = FALSE)


stopifnot(identical(rownames(mat_1), ext_data$name))



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

######################################################      Plotting     ######################################################### 


my_condition1 <- which(unique(SampleTable$Genotype) == SampleTable$Genotype[SampleTable$LibraryID == my_sample_id1])

my_color <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition1]

my_color_palette <- colorRampPalette(c("white",my_color,"black"))((length(matrix_files)+2))[c(-1,-(length(matrix_files)+2))]


my_grouping_names <- c("Qs.gene_length")


i=1

################################################################################


for(i in seq_along(my_grouping_names)){
        
        
        pdf(file = output_file, width = 7.5, height = 4.5)
        
        ################################################################################
        
        for(gi in  c("q1","q2","q3","q4","q5")){
                
                for(stat in c("Mean")){
                        
                        
                        par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))   
                        
                        ################################################################################
                        
                        
                        plot(colMeans(mat_1), type = "n", ylim = c(0,ifelse(stat == "Mean",5,2.5)),
                             xlab = "", ylab = "", xaxt = "n", lwd = 2,
                             col = my_color)
                        
                        axis(side = 1, at = c(1, ncol(mat_1)/2, ncol(mat_1)), labels = c("-3kb","TSS","+3kb"))
                        
                        
                        for(ti in seq_along(matrix_files)){
                                
                                mat_t <- readRDS(matrix_files[ti])
                                
                                if(stat == "Mean"){
                                        
                                        ydata <- colMeans(mat_t[rownames(mat_t) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])

                                } else if(stat == "Median"){
                                        
                                        ydata <- colMedians(mat_t[rownames(mat_t) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                }
                                
                                ydata <- zoo::rollmean(c(rep(NA,250), ydata, rep(NA,250)), k = 501)
                                
                                lines(ydata, lwd = 2, col = my_color_palette[ti])
                        }
                        
                        
                        
                        ################################################################################
                        
                        mtext(text = paste(stat,"Coverage"), side = 2, line = 3, cex = 1.2)
                        mtext(text = paste0(my_grouping_names[i], " - ", gi), side = 3, line = 0.5, cex = 1.2)

                        mtext(text = gsub(".*composite.|.pdf", "", output_file), 
                              side = 3, line = -1, cex = 1.25, font = 2, outer = TRUE)
                        
                        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(1, 0, 0, 0), new = TRUE)
                        
                        plot.new()
                        
                        
                        legend("bottom", horiz = TRUE, lwd = 5, cex = 0.8, bty = "n",
                               legend = gsub(".*\\/|_TSS.*", "", matrix_files),
                               col = my_color_palette)
                        
                        ################################################################################
                        
                }
        }
}

dev.off()




##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################
