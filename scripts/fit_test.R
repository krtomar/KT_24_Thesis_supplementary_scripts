

##################################################################################################################################
##################################################################################################################################


library(ShortRead)
library(rtracklayer)
library(GenomicFeatures)

library(RColorBrewer)
library(ggplot2)

library(nlme)


##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)


input_anno <- grep("gtf", args, value = TRUE)
input_table <- grep("SampleTable", args, value = TRUE)

norm_type <- grep("^Scer$|^Spom$", args, value = TRUE)

output_file <- grep("paramtests", args, value = TRUE)
output_dir <- gsub("\\/paramtests.*","", output_file)




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



formulaLogi <- function(x, a, b, xp){
        
        y <- a/(1+exp(-b*(x-xp)))       
        
        y
}


formulaExpo <- function(x, a, b, xp){
        
        y <- (a/b)*(1-exp(-b*(x-xp)))    
        
        y
}


formulaLog2Expo <- function(x, a, b, xp){
        
        y <- log2(a) - log2(b) + log2(1 - exp(-b*(x-xp)))
        
        y
}



##################################################################################################################################
##################################################################################################################################


log2type     <- strsplit(gsub(".*paramtests.|.txt","", output_file), "\\.")[[1]][1]

my_data_type <- strsplit(gsub(".*paramtests.|.txt","", output_file), "\\.")[[1]][2]

my_formula   <- strsplit(gsub(".*paramtests.|.txt","", output_file), "\\.")[[1]][3]


print(my_data_type)
print(log2type)
print(my_formula)


##################################################################################################################################
##################################################################################################################################



if(my_formula == "logistic"){
        
        formulaNLS <- formulaLogi
        
} else if(my_formula == "expo"){
        
        formulaNLS <- formulaExpo
}



##################################################################################################################################
##################################################################################################################################



my_data <- read.delim(paste0(output_dir, "/log2_", my_data_type, ".txt"), row.names = 1)

my_gene_ids <- rownames(my_data)


if(log2type == "cnt"){
        
        my_data <- 2^my_data
} 



##################################################################################################################################
##################################################################################################################################



my_fit_data <- data.frame(gene_id = my_gene_ids,
                          gene_name = my_genes$gene_name[match(my_gene_ids, my_genes$gene_id)],
                          matrix(nrow = length(my_gene_ids), ncol = 12))

#my_fit_data <- my_fit_data[1:100,]

my_gene_id <- "YFL039C"


################################################################################        


#for(i in 1:nrow(my_fit_data)){
#for(i in 1:100){


fit_res <- parallel::mclapply(1:nrow(my_fit_data), mc.cores = 16, FUN = function(i){ 
        
        
        my_gene_id <- my_fit_data$gene_id[i]
        
        match_samples <- match(gsub("\\-.*","", colnames(my_data)),SampleTable$LibraryID)
        
        df <- data.frame(x =      SampleTable$Time[match_samples],
                         y = as.numeric(my_data[rownames(my_data) ==  my_gene_id,]),
                         genotype =  SampleTable$Genotype[match_samples], 
                         dataset =   SampleTable$Dataset[match_samples])
        
        df$genotype <- factor(df$genotype, levels = rev(unique(SampleTable$Genotype)))
        df$dataset <- factor(df$dataset, levels = unique(df$dataset))
        
        
        a0  = ifelse(my_formula == "logistic", max(df$y), max(df$y)*0.05)
        b0  = ifelse(my_formula == "logistic", 0.2, 0.05)
        xp0 = ifelse(my_formula == "logistic", 10, 0)
        
        ################################################################################        
        
        my_fits  <- tryCatch(expr = {
                
                
                set.seed(123)
                
                
                fitv1 <- nls(y ~ formulaNLS(x, a, b, xp),
                             data =  df,
                             start = list(a = a0, b = b0, xp = xp0), 
                             control =  nls.control(maxiter = 1000, tol = 0.1),
                             subset = genotype == levels(df$genotype)[1])
                
                
                fitv2 <- nls(y ~ formulaNLS(x, a, b, xp),
                             data =  df,
                             start = list(a = a0, b = b0, xp = xp0), 
                             control =  nls.control(maxiter = 1000, tol = 0.1),
                             subset = genotype == levels(df$genotype)[2])
                
                
                fitfull0 <- nls(y ~ formulaNLS(x, a[genotype], b[genotype], xp[genotype]), 
                                data = df, 
                                control =  nls.control(maxiter = 1000, tol = 0.1),
                                start =  as.data.frame(rbind(coef(fitv1), coef(fitv2))))
                
                
                fitnull <- gnls(y ~ formulaNLS(x, a, b, xp),
                                data =  df,
                                start =  list(a  = mean(coef(fitfull0)[1:2]), 
                                              b  = mean(coef(fitfull0)[3:4]), 
                                              xp = mean(coef(fitfull0)[5:6])),
                                control =  gnlsControl(maxIter = 1000, nlsTol =  0.01),
                                weights = varPower())
                
                
                fitfull <- gnls(y ~ formulaNLS(x, a, b, xp), 
                                data = df, 
                                start =   list(a  = as.numeric(coef(fitfull0)[1:2]), 
                                               b  = as.numeric(coef(fitfull0)[3:4]), 
                                               xp = as.numeric(coef(fitfull0)[5:6])), 
                                control =  gnlsControl(maxIter = 1000, nlsTol =  0.01),
                                weights = varPower(),
                                params = list(a ~ genotype - 1, b ~ genotype -1, xp ~ genotype -1))
                
                
                fita <- gnls(y ~ formulaNLS(x, a, b, xp), 
                                data = df, 
                                start =   list(a  =       mean(coef(fitfull)[1:2]), 
                                               b  = as.numeric(coef(fitfull)[3:4]), 
                                               xp = as.numeric(coef(fitfull)[5:6])), 
                                control =  gnlsControl(maxIter = 1000, nlsTol =  0.01),
                                weights = varPower(),
                                params = list(a ~ 1, b ~ genotype -1, xp ~ genotype -1))
                
                
                fitb <-  gnls(y ~ formulaNLS(x, a, b, xp), 
                                 data = df, 
                                 start =   list(a  = as.numeric(coef(fitfull)[1:2]), 
                                                b  =       mean(coef(fitfull)[3:4]), 
                                                xp = as.numeric(coef(fitfull)[5:6])), 
                                 control =  gnlsControl(maxIter = 1000, nlsTol =  0.01),
                                 weights = varPower(),
                                 params = list(a ~ genotype - 1, b ~ 1, xp ~ genotype -1))
                
                
                fitxp <-  gnls(y ~ formulaNLS(x, a, b, xp), 
                               data = df, 
                               start =   list(a  = as.numeric(coef(fitfull)[1:2]), 
                                              b  = as.numeric(coef(fitfull)[3:4]), 
                                              xp =       mean(coef(fitfull)[5:6])), 
                               control =  gnlsControl(maxIter = 1000, nlsTol =  0.01),
                               weights = varPower(),
                               params = list(a ~ genotype - 1, b ~ genotype -1, xp ~ 1))
                
                
                c(coef(fitfull), 
                  AIC(fitfull),
                  summary(fitfull)$modelStruct$varStruct[1],
                  anova(fitnull, fitfull)[2,"p-value"],
                  anova(fita,    fitfull)[2,"p-value"],
                  anova(fitb,    fitfull)[2,"p-value"],
                  anova(fitxp,   fitfull)[2,"p-value"])
                
        }, error = function(c) rep(NA, 12))
        
        
        ################################################################################        
        
        my_fits
        
        #my_fit_data[i,c(-1,-2)] <- my_fits
        
})


################################################################################        


my_fit_data[,c(-1,-2)] <- Reduce(rbind, fit_res)


colnames(my_fit_data)[c(-1,-2)] <-  c(paste(rep(c("a","b","xp"),each = 2),rev(unique(SampleTable$Genotype)), sep = "_"),
                                      "aic", "power_weight","pval_overall", "pval_a", "pval_b", "pval_xp")

my_fit_data$padj_overall <- p.adjust(my_fit_data$pval_overall, method = "BH")
my_fit_data$padj_a       <- p.adjust(my_fit_data$pval_a,       method = "BH")
my_fit_data$padj_b       <- p.adjust(my_fit_data$pval_b,       method = "BH")
my_fit_data$padj_xp      <- p.adjust(my_fit_data$pval_xp,      method = "BH")


################################################################################        




##################################################################################################################################
##################################################################################################################################




write.table(my_fit_data, file = paste0(output_dir,"/paramtests.",log2type,".",my_data_type,".",my_formula,".txt"), 
            quote = F, sep = "\t", row.names = FALSE)




##################################################################################################################################
##################################################################################################################################















##################################################################################################################################
##################################################################################################################################

