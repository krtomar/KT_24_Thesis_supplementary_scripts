





my_gene_length <- width(my_genes)/1e3
names(my_gene_length) <- my_genes$gene_id

my_gene_length <-  my_gene_length[match(rownames(my_counts_Scer), names(my_gene_length)) ]


if(all(identical(names(my_gene_length),    rownames(my_counts_Scer)),
       identical(names(n_genic_Spom_reads), colnames(my_counts_Scer)),
       identical(names(n_all_Spom_reads), colnames(my_counts_Scer)),
       identical(SampleTable$LibraryID, colnames(my_counts_Scer)))){
    
    #norm_counts <- t(t(my_counts_Scer)/n_genic_Spom_reads)*1e6 #/ my_gene_length
    norm_counts <- t(t(my_counts_Scer)/n_all_Spom_reads)*1e6 #/ my_gene_length
    
} else {
    
    rm(list = "norm_counts")
}


minmax_scale <- function(x){ ((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE)-min(x, na.rm=TRUE))) }



library(nlme)


my_fav_genes <- c("ACT1","PGK1","PMA1","RAD51","RDN25-1")

xn <- seq(0,45,0.01)



pdf("test_fit.pdf", width = 11, height = 6)


par(mfcol=c(2,4), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2.5,1,0))


for(my_fav_gene in my_fav_genes){
    
    
    for(i in seq(1,ncol(norm_counts),7)){
        
        x <- SampleTable$Time[i:(i+6)] 
        
        if(grepl("RDN", my_fav_gene)){
            
            y <- as.numeric(norm_counts[rownames(norm_counts) ==  my_fav_gene,(i:(i+6))])
            
        } else {
            y <- as.numeric(norm_counts[rownames(norm_counts) == my_genes$gene_id[my_genes$gene_name %in% my_fav_gene],(i:(i+6))])
            
        }
        
        #y <- y / max(y)
        #y <- minmax_scale(y)
        
        if(grepl("pull",SampleTable$Conditions[i])){
            
            # fit <- nls(y ~ SSlogis(x, Asymptote, Inflection, Scale),
            #            control = nls.control(maxiter = 200))
            
            asym <- max(y)
            
            fit <- nls(y ~ asym/(1+exp((xmid-x)/scal)),
                       start = c(xmid = 20, scal =5),
                       control = nls.control(maxiter = 200))
            
            # asym <- as.numeric(coef((fit))["Asymptote"])
            # xmid <- as.numeric(coef((fit))["Inflection"])
            # scal <- as.numeric(coef((fit))["Scale"])
            
            xmid <- as.numeric(coef((fit))["xmid"])
            scal <- as.numeric(coef((fit))["scal"])
            
            yh <- asym/(1+exp((xmid-xn)/scal))
            
            
            
            plot(xn, yh, main = unique(SampleTable$Conditions[i:(i+6)]),
                 xlab = "Time", ylab = "Normalized Counts",
                 type="l", lwd=1.5, col = "blue2", ylim = c(0, 1.2*max(y)))
            points(x,y, pch=19, col="red3", cex=0.9)
            
            
            slope <- asym * (1/scal) / 4
            inter <- (asym/2) - slope*xmid
            
            #stopifnot(identical(floor(slope), floor(max(yh[-1] - yh[-length(yh)])*100)))
            
            abline(coef = c(inter, slope), lty=2)
            
            abline(h=asym, lty=2)
            text(x = 0, y = asym, label = paste("Asymptote = ", round(asym,1)), adj=c(0,-1), cex=0.8)
            
            abline(v=xmid, lty=2)
            text(x = ifelse(xmid > 0, xmid, 10), 
                 y = 0, label = paste("Inflection = ", round(xmid,1)), adj=c(-0.25,0), cex=0.8)
            
            abline(v=xmid, lty=2)
            text(x = ifelse(xmid > 0, xmid, 10),
                 y = asym/2, label = paste0(" Scale = ", round(scal,1),
                                                     "\nSlope = ", round(slope,4)), 
                 adj=c(-0.25,0.5), cex=0.8)
            
            
        } else {
            
            fit <- lm(y ~ x)
            
            b0 <- coef((fit))[1]
            b1 <- coef((fit))[2]
            
            
            yh <- b0+b1*xn
            
            
            plot(xn, yh, main = SampleTable$Conditions[i],
                 xlab = "Time", ylab = "Normalized Counts",
                 type="l", col = "blue3", ylim = c(0, max(y)*1.1))
            points(x,y, pch=19, col="red3", cex=0.9)
            
        }
        
        
        # summary(fit)
        # AIC(fit)
        
        
        if(i == 1){mtext(text = my_fav_gene, outer = TRUE, font=2)}
        
        
    }
}

dev.off()













fit <- nls(y ~ a*exp(-exp(b-g*x)+d*x),
           start = c(a = 1, b = 1, g = 0.1, d = 0.01),
           control = nls.control(maxiter = 200))

a <- as.numeric(coef((fit))["a"])
b <- as.numeric(coef((fit))["b"])
g <- as.numeric(coef((fit))["g"])
d <- as.numeric(coef((fit))["d"])


fit <- nls(y ~ a*exp(-b*exp(-c*x)),
           start = c(a = 1, b = 1, c = 0.01),
           control = nls.control(maxiter = 200))


a <- as.numeric(coef((fit))["a"])
b <- as.numeric(coef((fit))["b"])
c <- as.numeric(coef((fit))["c"])


yh <-  a*exp(-b*exp(-c*xn))



plot(xn, yh, 
     main = unique(SampleTable$Conditions[SampleTable$LibraryID %in% gsub("\\-.*","",colnames(my_plot_data)[i:(i+6)])]),
     xlab = "", ylab = "",
     type="l", lwd=1.5, col = "blue2", ylim = c(min(y),max(y) + abs(0.05*max(y))))

points(x,y, pch=19, col="red3", cex=0.9)



