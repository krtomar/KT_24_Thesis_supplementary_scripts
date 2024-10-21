library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggbreak)
library(patchwork)
rm(list=ls())
my_color_palette <- c("#E69F00","#999999","#56B4E9","#009E73")
my_color_palettes <- c('orange2','grey80','olivedrab4','purple4','orange4','grey9','orangered3','grey57')

# to print all the spom reads norm bt ecoli reads,  plotting counts pombe per million E coli counts.
my_scer <- read.delim('tables/raw_counts_Scer.txt',stringsAsFactors = FALSE,header = TRUE)
my_spom <- read.delim('tables/raw_counts_Spom.txt',stringsAsFactors = FALSE,header = TRUE)
SampleTable <- read.delim('tables/SampleTable.txt', header=T, stringsAsFactors = F)
SampleTable$Genotype <- gsub('^CHD1$','chd1',SampleTable$Genotype)
genotype <- unique(SampleTable$Genotype)
ID <- ifelse(genotype[1] == 'CHD1_AID','KripiID','LibraryID')
exp <- c('chd1','WT','SPT6_AA_Rap','SPT6_AA_Control','isw1isw2chd1','isw1isw2','CHD1_AID','Vehicle')
i=1

while(exp[i] != genotype[1]){
  print(i)
  i=i+1
}
cg <- i  
############### for calcultaing % non specifc scer pulldown by spom spike in at 0 minutes label############### 
scer_0 <- colSums(my_scer[-(my_scer$X == "Ecoli"),
                                 colnames(my_scer) %in% SampleTable[,ID][(SampleTable$Time == 0) & 
                                 (SampleTable$Fraction == 'pull out fraction')]])
zero_pull <- rbind(scer_0,
                   colSums(my_spom[,colnames(my_spom) %in% SampleTable[,ID][(SampleTable$Time == 0) 
                   & (SampleTable$Fraction == 'pull out fraction')]]))
zero_pull <- rbind(zero_pull,zero_pull[1,]/zero_pull[2,])
rownames(zero_pull) <- c('scer','spom','sc_by_sp')
#View(zero_pull)
mean_sc_by_sp_0min <- mean(zero_pull['sc_by_sp',])
####################################################################################
ecoli_reads <- colSums((my_scer[my_scer$X == "Ecoli",][,-1]))/1e6
spom_reads <- colSums(my_spom[,-1])/1e6
frac <- data.frame(row.names = colnames(my_spom)[-1],
                      fraction=(colSums(my_spom[,-1])/1e6)/(colSums((my_scer[my_scer$X == "Ecoli",][,-1]))/1e6),
                      ecoli=colSums((my_scer[my_scer$X == "Ecoli",][,-1]))/1e6,
                      spom=colSums(my_spom[,-1])/1e6,
                      scer = colSums((my_scer[!(my_scer$X == "Ecoli"),][,-1]))/1e6,
                   frac_sp_scer = (colSums(my_spom[,-1])/1e6)/(colSums((my_scer[!(my_scer$X == "Ecoli"),-1])/1e6)))

frac[,1] <-ifelse(is.infinite(frac[,1]), 0, frac[,1]) #change infinite value to finite

if(all(identical(colnames(my_scer[,-1]), SampleTable[,ID]),
       identical(colnames(my_spom[,-1]), SampleTable[,ID]))){
        print('ifworks')
       SampleTable$Conditions <- paste(SampleTable$Genotype, gsub(" .*", "", SampleTable$Fraction), sep = "_")
       my_conditions <-factor(SampleTable$Conditions, levels = unique(SampleTable$Conditions))
        frac$time <- SampleTable$Time
        frac$conditions <- my_conditions
        frac$dataset <- SampleTable$Dataset
        frac$dataset <- factor(frac$dataset)#, levels= unique(frac$dataset))
        }
pull_data <- frac[grepl('_pull',frac$conditions),]
pull_data$dataset <- factor(pull_data$dataset)

total_data <- frac[grepl('_total',frac$conditions),]
total_data$conditions <- droplevels(total_data$conditions) #drop extra levels from pullout data
total_data$dataset <- factor(total_data$dataset)
 
pdf(file="plots/spom_ecoli_fit.pdf", width = 5, height = 5)
par(mfrow=c(1, 1), mar=c(4,4,4,4), mgp = c(2,1,0), oma=c(2,2,2,2),xpd=F)


##########################################################################
##########################################################################
pull_data$dataset <- factor(pull_data$dataset) #, levels = unique(pull_data$dataset))
pull_data$conditions <- factor(pull_data$conditions, levels = unique(pull_data$conditions))
fit1 <- lm(fraction ~ time*conditions + dataset , data = pull_data) #to take into account the different condition(genotypes)
#summary(fit1)

plot(pull_data$time, pull_data$fraction, 
     col = my_color_palettes[cg:(cg+1)][pull_data$conditions], 
     ylab= 'spom/ecoli reads ratio',
     ylim= c(0,500),
     xlab = 'time of label in min',
     pch = (15:17)[pull_data$dataset], cex = 0.75,
     main = "fraction ~ (time * conditions + datasets)")
# Add legend to differentiate between groups
legend("topright", legend = unique(pull_data$conditions), col = my_color_palettes[cg:(cg+1)][unique(pull_data$conditions)], pch = (15:17)[pull_data$dataset])

# Predict mpg for a sequence of horsepower values for both groups
timeol <- seq(min(pull_data$time), max(pull_data$time), length.out = 100)
selected_dataset <- levels(pull_data$dataset)[1]
# For am = 0 (automatic)
pred_gn1 <- predict(fit1, newdata = data.frame(time = timeol, conditions = unique(pull_data$conditions)[1],
                     dataset = selected_dataset))
# For am = 1 (manual)
pred_gn2 <- predict(fit1, newdata = data.frame(time = timeol, conditions = unique(pull_data$conditions)[2],
                   dataset = selected_dataset))

# Add regression lines for both groups
lines(timeol, pred_gn1, col = my_color_palettes[cg][unique(pull_data$conditions)], lwd = 2)  # Line for gn1
lines(timeol, pred_gn2, col = my_color_palettes[(cg+1)][unique(pull_data$conditions)], lwd = 2)  # Line for gn2  

############################for totals spom/scer ratio for same amount of spom added due to cell counting #####################
fit2 <- lm(frac_sp_scer ~ time*conditions + dataset , data = total_data) #to take into account the different condition(genotypes)
#summary(fit2)
plot(total_data$time, total_data$frac_sp_scer, 
     col = my_color_palettes[cg:(cg+1)][total_data$conditions], 
     ylab= 'spom/scer reads ratio in totals',
     ylim= c(0,1.5),
     xlab = 'time of label in min',
     pch = (15:17)[total_data$dataset], cex = 0.75,
     main = "fraction ~ (time * conditions + datasets)")
# Add legend to differentiate between groups
legend("topright", legend = unique(total_data$conditions), col = my_color_palettes[cg:(cg+1)][unique(total_data$conditions)], pch = (15:17)[total_data$dataset])

# Predict fit values for both groups
timeol1 <- seq(min(total_data$time), max(total_data$time), length.out = 100)
selected_dataset1 <- levels(total_data$dataset)[1]
#for genotype 1
pred_gn11 <- predict(fit2, newdata = data.frame(time = timeol1, conditions = unique(total_data$conditions)[1],
                                               dataset = selected_dataset1))
# for genotype 2
pred_gn22 <- predict(fit2, newdata = data.frame(time = timeol1, conditions = unique(total_data$conditions)[2],
                                               dataset = selected_dataset1))

# Add regression lines for both groups
lines(timeol1, pred_gn11, col = my_color_palettes[cg][unique(total_data$conditions)], lwd = 2)  # Line for gn1
lines(timeol1, pred_gn22, col = my_color_palettes[(cg+1)][unique(total_data$conditions)], lwd = 2)  # Line for gn2  

fit_data <- data.frame(fraction= character(),DKOvsTKO_estimate=double(),DKOvsTKO_padj=double(), 
                       DKOvsTKO_estimate_over_time=double(), DKOvsTKO_over_time_padj=double(), 
                       model_R_squared_adj=double(),model=character())
fit_data[1:2,1] <- c('pull out_spom_by_ecoli_reads_perMillion','total RNA_spom_by_scer_reads_perMillion')
fit_data[1:2,2] <- c(summary(fit1)$coefficients[3,1], summary(fit2)$coefficients[3,1]) 
fit_data[1:2,3] <- c(p.adjust(summary(fit1)$coefficients[3,4],method='BH'), p.adjust(summary(fit2)$coefficients[3,4],method = 'BH')) 
fit_data[1:2,4] <- c(summary(fit1)$coefficients[6,1], summary(fit2)$coefficients[6,1]) 
fit_data[1:2,5] <- c(p.adjust(summary(fit1)$coefficients[6,4],method='BH'), p.adjust(summary(fit2)$coefficients[6,4],method = 'BH')) 
fit_data[1:2,6] <- c(summary(fit1)$adj.r.squared, summary(fit2)$adj.r.squared) 
fit_data[1:2,7] <- c(as.character(summary(fit1)$call[2]), as.character(summary(fit2)$call[2])) 

write.table(fit_data, file = paste0('plots/spike-in_levels_',genotype[1],'_vs_',genotype[2],".txt"),
            quote = F, sep = "\t", row.names = FALSE)


################################ for spom raw reads ###################
fit11 <- lm(spom ~ time*conditions + dataset , data = total_data) #to take into account the different condition(genotypes)
#summary(fit1)
plot(total_data$time, total_data$spom, 
     col = my_color_palettes[cg:(cg+1)][total_data$conditions], 
     ylab= 'spom reads in totals',
     #ylim= c(0,500),
     xlab = 'time of label in min',
     pch = (15:17)[total_data$dataset], cex = 0.75,
     main = "fraction ~ (time * conditions + datasets)")
# Add legend to differentiate between groups
legend("topright", legend = unique(total_data$conditions), col = my_color_palettes[cg:(cg+1)][unique(total_data$conditions)], pch = (15:17)[total_data$dataset])

# Predict mpg for a sequence of horsepower values for both groups
timeol1 <- seq(min(total_data$time), max(total_data$time), length.out = 100)
selected_dataset1 <- levels(total_data$dataset)[1]
# For am = 0 (automatic)
pred_gn11 <- predict(fit11, newdata = data.frame(time = timeol1, conditions = unique(total_data$conditions)[1],
                                                 dataset = selected_dataset1))
# For am = 1 (manual)
pred_gn22 <- predict(fit11, newdata = data.frame(time = timeol1, conditions = unique(total_data$conditions)[2],
                                                 dataset = selected_dataset1))

# Add regression lines for both groups
lines(timeol1, pred_gn11, col = my_color_palettes[cg][unique(total_data$conditions)], lwd = 2)  # Line for gn1
lines(timeol1, pred_gn22, col = my_color_palettes[(cg+1)][unique(total_data$conditions)], lwd = 2)  # Line for gn2  


### The estimate for "Intercept" indicates when "time" = 0 for dataset "d1"
### "datasetsd2" and "datasetsd3" indicate the intercept difference relative to dataset "d1"
### "time" estimate is the slope of the line (in this model, there is one slope and 3 intercepts)
### Pr(>|t|) is the p-value to test whether the estimates are different from zero

dev.off()


################################################################################################################
###################################### to plot simple barplots with fraction spombe /ecoli########################
pdf("plots/barplot_QC2.pdf", height = 7, width = 14, useDingbats = F)
par(mfrow=c(1,3), mar = c(6,6,1,1), oma = c(6,3,6,3), mgp = c(2.5,1,0))

cond_pull <- frac$conditions[grepl('pull',frac$conditions)]
cond_tot <- frac$conditions[grepl('total',frac$conditions)]
  ###############
  bp1 <-  barplot(frac$ecoli[grepl('pull',frac$conditions)],
                 ylab = "Number of Ecoli reads in Millions", las = 3, xaxt = "n",
                 col = my_color_palette[cond_pull])
                 #legend=levels(frac$conditions))
  axis(side = 1, at = bp1, labels = rownames(frac)[grepl('pull',frac$condition)], cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)
  legend("topright", legend = unique(cond_pull), horiz = FALSE,
         bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[unique(cond_pull)],inset=c(0,-0.04),xpd=T) #seq along to genotype name in order genotype 
  
  #########################
  bp2 <-  barplot(frac$spom[grepl('pull',frac$conditions)],
                  ylab = "Number of S.pombe reads in Millions", las = 3, xaxt = "n",
                  col= my_color_palette[frac$conditions[grepl('pull',frac$conditions)]]) #color in order of the levels
                  #legend=levels(frac$conditions))
  axis(side = 1, at = bp2, labels = rownames(frac)[grepl('pull',frac$conditions)], cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)
  legend("topright", legend = unique(cond_pull), horiz = FALSE,
          bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[unique(cond_pull)],inset=c(0,-0.04),xpd=T) #seq along to genotype name in order of conditions

  
  #########################
  
   bp <-  barplot(frac$fraction[grepl('pull',frac$conditions)],
         ylab = "Number of S.pombe reads per Million Ecoli ", las = 3, xaxt = "n",
         col= my_color_palette[frac$conditions[grepl('pull',frac$conditions)]])
         #legend=unique(frac$conditions[grepl('pull',frac$conditions)]))
         #legend("topright", legend=levels(frac$conditions[grepl('pull',frac$conditions)]) ,
             #   col=unique(my_color_palette[frac$conditions[grepl('pull',frac$conditions)]])))
axis(side = 1, at = bp, labels = rownames(frac)[grepl('pull',frac$condition)], cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)
legend("topright", legend = unique(cond_pull), horiz = FALSE,
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[unique(cond_pull)],inset=c(0,-0.04),xpd=T) #seq along to genotype name in order of conditions


# Identify samples for genotype1 and genotype2 based on patterns in frac$conditions
genotype1_samples <- frac$dataset[grepl(paste0(genotype[1],'_total'), frac$conditions)]
genotype2_samples <- frac$dataset[grepl(paste0(genotype[2],'_total'), frac$conditions)]

# Ensure that the datasets are ordered by replicates (if applicable)
# Assuming replicates are identifiable within frac$dataset
genotype1_samples_ordered <- sort(genotype1_samples)  # Sorting within genotype1
genotype2_samples_ordered <- sort(genotype2_samples)  # Sorting within genotype2

# Combine the ordered samples: genotype1 first, followed by genotype2
all_samples_ordered <- c(genotype1_samples_ordered, genotype2_samples_ordered)

# Now use this ordered list for x-axis labels in the barplot
bp3 <- barplot(frac$frac_sp_scer[grepl('total', frac$conditions)],
               ylab = "Number of S.p reads per Million Sc", las = 3, xaxt = "n",
               col = my_color_palette[frac$conditions[grepl('total', frac$conditions)]])

# Set axis labels with the ordered samples
axis(side = 1, at = bp3, labels = all_samples_ordered, cex.axis = 0.5, las = 3, lwd.ticks = 0, lwd = 0)

legend("topright", legend = unique(cond_tot), horiz = FALSE,
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[unique(cond_tot)], inset = c(0, -0.04), xpd = T)

################################################################################################################
bp3 <-  barplot(frac$frac_sp_scer[grepl('total',frac$conditions)],
               ylab = "Number of S.pombe reads per Million Scer ", las = 3, xaxt = "n",
               col= my_color_palette[frac$conditions[grepl('total',frac$conditions)]])

axis(side = 1, at = bp3, labels = rownames(frac)[grepl('total',frac$condition)], cex.axis=0.5, las = 3, lwd.ticks = 0, lwd=0)
legend("topright", legend = unique(cond_tot), horiz = FALSE,
       bty = "n", cex = 0.8, pt.cex = 1, fill = my_color_palette[unique(cond_tot)],inset=c(0,-0.04),xpd=T) #seq along to genotype name in order of conditions



dev.off()

