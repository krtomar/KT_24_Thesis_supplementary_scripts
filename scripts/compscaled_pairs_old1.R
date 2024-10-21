






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



my_sample_id1 <- gsub(".*\\/|_GB.*", "", matrix_files[1])
my_sample_id2 <- gsub(".*\\/|_GB.*", "", matrix_files[2])

mat_1 <- readRDS(matrix_files[1])
mat_2 <- readRDS(matrix_files[2])


mat_1 <- mat_1[(rownames(mat_1) %in% rownames(log2_norm_counts)),]
mat_2 <- mat_2[(rownames(mat_2) %in% rownames(log2_norm_counts)),]


high_extremes <- c("YGR190C","YGR191W","YGR192C","YGR193C","YGR194C","RDN25-1","YLR154W-A","YLR154W-B","YLR154W-C","RDN58-1",  
                   "YLR154W-E","YLR154W-F","RDN18-1","RDN5-1","RDN25-2","RDN58-2","RDN18-2","RDN5-2","YLR154C-H","YLR155C", 
                   "YLR160C","YLR161W","RDN5-6", "YLR162W","YLR162W-A","YLR163C","YLR163W-A","YPR078C","YPR079W","YPR080W",  
                   "YPR081C")

mat_1 <- mat_1[!(rownames(mat_1) %in% high_extremes),]
mat_2 <- mat_2[!(rownames(mat_2) %in% high_extremes),]




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
my_condition2 <- which(unique(SampleTable$Genotype) == SampleTable$Genotype[SampleTable$LibraryID == my_sample_id2])

my_color1 <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition1]
my_color2 <- c("#E69F00","#999999","#56B4E9","#009E73")[my_condition2]



my_grouping_names <- c("Qs.gene_length")


i=1

################################################################################


for(i in seq_along(my_grouping_names)){
        
        ydata1 <- colMeans(mat_1[rownames(mat_1) %in% ext_data$name[ext_data[,my_grouping_names[i]] == "q1"],])
        ydata2 <- colMeans(mat_2[rownames(mat_2) %in% ext_data$name[ext_data[,my_grouping_names[i]] == "q1"],])
        
        ydata1 <- zoo::rollmean(c(rep(NA,100), ydata1, rep(NA,100)), k = 201)
        ydata2 <- zoo::rollmean(c(rep(NA,100), ydata2, rep(NA,100)), k = 201)
        
        ylim <- ceiling(max(c(ydata1, ydata2)*1.1, na.rm = TRUE)*10)/10
        
        
        pdf(file = output_file, width = 7.5, height = 4.5)
        
        
        ################################################################################
        
        for(gi in  c("q1","q2","q3","q4","q5")){
                
                for(stat in c("Mean")){
                        
                        
                        par(mfrow=c(1,2), oma=c(2,2,2,2), mar=c(4,4,4,2))   
                        
                        ################################################################################
                        
                        if(stat == "Mean"){
                                
                                ydata1 <- colMeans(mat_1[rownames(mat_1) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                ydata2 <- colMeans(mat_2[rownames(mat_2) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                
                        } else if(stat == "Median"){
                                
                                ydata1 <- colMedians(mat_1[rownames(mat_1) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                                ydata2 <- colMedians(mat_2[rownames(mat_2) %in% ext_data$name[ext_data[,my_grouping_names[i]] == gi],])
                        }
                        
                        
                        ydata1 <- zoo::rollmean(c(rep(NA,100), ydata1, rep(NA,100)), k = 201)
                        ydata2 <- zoo::rollmean(c(rep(NA,100), ydata2, rep(NA,100)), k = 201)
                        
                        ################################################################################
                        
                        
                        plot(ydata1, type = "l", ylim = c(0, ylim),
                             xlab = "", ylab = "", xaxt = "n", lwd = 2,
                             col = my_color1)
                        
                        axis(side = 1, at = c(1, ncol(mat_1)/4,ncol(mat_1)*3/4, ncol(mat_1)), labels = c("-1kb","TSS","TTS","+1kb"))
                        
                        lines(ydata2, lwd = 2,col = my_color2)
                        
                        ################################################################################
                        
                        mtext(text = paste(stat,"Coverage"), side = 2, line = 3, cex = 1.2)
                        mtext(text = paste0(my_grouping_names[i], " - ", gi), side = 3, line = 0.5, cex = 1.2)

                        mtext(text = gsub(".*composite.|.pdf", "", output_file), 
                              side = 3, line = -1, cex = 1.25, font = 2, outer = TRUE)
                        
                        par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(1, 0, 0, 0), new = TRUE)
                        
                        plot.new()
                        
                        
                        legend("bottom", horiz = TRUE, lwd = 5, cex = 0.8, bty = "n",
                               legend = c(my_sample_id1, my_sample_id2),
                               col    = c(my_color1,     my_color2))
                        
                        ################################################################################
                        
                }
        }
}

dev.off()




##################################################################################################################################
##################################################################################################################################














##################################################################################################################################
##################################################################################################################################
