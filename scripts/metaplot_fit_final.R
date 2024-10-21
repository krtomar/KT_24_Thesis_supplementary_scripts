
############ to calculate 5'/3' and mid-gene/3' ratio from average read coverage of gene length quintiles ################
##################################################################################################################################
##################################################################################################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggbreak)
library(patchwork)
library(tidyr)
library(nlme)
library(lawstat)
library(lmtest)
rm(list=ls())
##################################################################################################################################
##################################################################################################################################


my_color_palettes <- c('orange2','grey80','olivedrab4','purple4','orange4','grey9','orangered3','grey57')


SampleTable <- read.delim('tables/SampleTable.txt', header=T, stringsAsFactors = F)
SampleTable$Genotype <- gsub('^CHD1$','chd1',SampleTable$Genotype)
genotype <- unique(SampleTable$Genotype)
exp <- c('chd1','WT','SPT6_AA_Rap','SPT6_AA_Control','isw1isw2chd1','isw1isw2','CHD1_AID','Vehicle')
i=1

while(exp[i] != genotype[1]){
        print(i)
        i=i+1
}
cg <- i  
sample_id <- ifelse(genotype[1] == 'CHD1_AID','KripiID','LibraryID')
SampleTable <- SampleTable[!duplicated(SampleTable[,sample_id]), , drop = FALSE]
SampleTable <- SampleTable[order(SampleTable[,sample_id]), , drop = FALSE]


SampleTable$Dataset <- gsub(" ","_", SampleTable$Dataset)
SampleTable$Conditions <- paste(SampleTable$Genotype, gsub(" .*", "", SampleTable$Fraction), sep = "_")
SampleTable$fitGroup <- paste(SampleTable$Genotype, SampleTable$Dataset, sep = "_")


##################################################################################################################################
################################################## major table ################################################################################
N=210
# my_data <- data.frame(LibraryID= integer(N), Qs= character(N), TSS= numeric(N), Mid= numeric(N),
#                       TTS=numeric(N), log2_TSS_by_TTS =numeric(N), log2_Mid_by_TTS=numeric(N))
f_add <- ifelse(genotype[1] == 'chd1', '/analysis/221012_chd1vsWT', 
          ifelse(genotype[1] == 'SPT6_AA_Rap','/analysis/230110_SPT6AA',
          ifelse(genotype[1] == 'isw1isw2chd1','/analysis/230131_TKOvsDKO',
                 '/analysis/230417_Chd1depletion')))
f_add2 <- '/Output_240123/scaledcomp/wmulti' # read coverage of genes without high extremes and transposons

frac <- c('pullout','totalRNA')

for(fr in frac){
files <- list.files(paste0(f_add,f_add2), pattern = paste0("*",fr,".*transposons.txt"), full.names = TRUE)

# the following is a simpler solution:
my_data <- data.frame()

for(i in files){
        print(i)
        f <- read.delim(i, header = TRUE, stringsAsFactors = FALSE)
        my_data <- rbind(my_data, f)
}

colnames(my_data) <- c(sample_id, 'Qs', 'TSS', 'Mid','TTS', 'log2_TSS_by_TTS', 'log2_Mid_by_TTS')

my_data2 <- merge(my_data, SampleTable, by = sample_id) # change to KripiID for chd1 depln dataset
sig_data <- data.frame(Time=integer(),pval=double(),ratio_type=character(),Qs=character(),signf=character())

fol= 'metaplot_ratios'
dir.create(paste0('plots/',fol))
pdf(file =paste0("plots/metaplot_ratios/gene_length_metaplot_ratios_",fr,".pdf"),width = 6, height = 5,useDingbats = FALSE)
par(mfrow=c(1,1), mar = c(4,4,2,2), mgp = c(2.5,1,0), oma=c(2,2,2,2),xpd=F)
##################################################################################################################################
for(c in unique(my_data2$Qs)){
#c='q5' #quintile number

match_samples <- my_data2$Qs == c

### you could try to keep the log2 scale

df <- data.frame(x =  my_data2$Time[match_samples],
                 y1 = as.numeric(2^(my_data2$log2_TSS_by_TTS[match_samples])),
                 y2=  as.numeric(2^(my_data2$log2_Mid_by_TTS[match_samples])),
                 genotype =  my_data2$Genotype[match_samples], 
                 dataset =   my_data2$Dataset[match_samples])

df$genotype <- factor(df$genotype, levels = rev(unique(SampleTable$Genotype)))
df$dataset <- factor(df$dataset, levels = unique(df$dataset))
print(levels(df$genotype))
my_ylab= paste0('TSS/TTS_',c)
my_ylab2= paste0('Mid/TTS_',c)
plot(df$x, df$y1,main= 'TSS/TTS',col = rev(my_color_palettes[cg:(cg+1)])[df$genotype], 
     xlab = "Time", ylab = my_ylab,
     pch = (15:17)[df$dataset])

# please note the colors and the levels which are reversed relative to color scheme used

legend("topright", legend = levels(df$genotype), fill = rev(my_color_palettes[cg:(cg+1)]))

plot(df$x, df$y2, main= 'Mid/TTS' ,col = rev(my_color_palettes[cg:(cg+1)])[df$genotype], 
     xlab = "Time", ylab = my_ylab2,
     pch = (15:17)[df$dataset])

legend("topright", legend = levels(df$genotype), fill = rev(my_color_palettes[cg:(cg+1)]))

print(rev(my_color_palettes[cg:(cg+1)]))
print(levels(df$genotype))

      

###############################

# let´s try simple lm
# use time also as factors

df$time <- factor(df$x)
set.seed(123)
      for (y in df[,c('y1','y2')]){
      
      f_final <- function(df,y){
        f <- lm(y ~ dataset+genotype+time+genotype:time, data = df)
        #summary(f)$coefficients[2,4]
        } 
     
      
      fit <- f_final#ifelse(fn == 1,f1,f2 ) #ifelse(fn == 2,f2,f3))
      
      # in order to get other comparisons you need to relevel the variables
      
      ppt <- data.frame(matrix(nrow = 7,ncol=6))
      colnames(ppt) <- c('Time','pval','estimate','ratio_type','Qs','signf')
      t <- levels(df$time)
      ppt$Qs <- c
      ppt$ratio_type <- ifelse(unique(y == df$y1), 'TSS_by_TTS', 'Mid_by_TSS')
      
      #ppt$model <- paste0('f',fn)
            for( j in 1:7){
            #j=1  
            df$time <- relevel(df$time, ref = t[j])
            
            #reorder the levels of factor time from 0 to jth time point so all values will be compared to this now
            ppt$Time[j] <- t[j] #unique(SampleTable$Time)[j]
            #print(summary(fit(df=df,y=y)))
            ppt$pval[j] <- summary(fit(df=df,y=y))$coefficients[4,4] 
            ppt$signf[j] <- ifelse(ppt$pval[j] < 0.01, 'yes','no')
            ppt$estimate[j] <- summary(fit(df=df,y=y))$coefficients[4,1]
            }
        sig_data <- rbind(sig_data,ppt)
          }
} #for Quintiles
dev.off()
dev.list()

write.table(sig_data, file = paste0('plots/',fol,"/gene_length_metaplot_pval_",genotype[1],'_vs_',genotype[2],"_",fr,".txt"), 
            quote = F, sep = "\t", row.names = FALSE)
rm(files,sig_data,my_data,my_data2,ppt,fol)
} #for fractions pullout or totalRNA

