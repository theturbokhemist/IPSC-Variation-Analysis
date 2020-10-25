setwd("/Users/gordid/Desktop/IPSC_Project/Data/GSE79636_RAW/")
files <- list.files()

#import data
IPSC_Generation_Data <- as.data.frame(read.csv(file = "/Users/gordid/Desktop/IPSC_Project/Data/IPSC_Generation_CSV.csv"))
IPSC_Generation_Data$CLONE.ID <- as.character(IPSC_Generation_Data$CLONE.ID)
IPSC_Generation_Data <- IPSC_Generation_Data[gtools::mixedorder(IPSC_Generation_Data$CLONE.ID),]

Subject_Demo_Data <- as.data.frame(read.csv(file = "/Users/gordid/Desktop/IPSC_Project/Data/Subject_Demo_CSV.csv"))
Subject_Demo_Data <- Subject_Demo_Data[,1:5]

#extract group info and frequency of replicates for each group
Pat.ID. <- sub(pattern = "_.*", replacement = "", x = IPSC_Generation_Data[, "CLONE.ID"])
length(Pat.ID.)
freq <- plyr::count(Pat.ID.)[, 2]

#remove all groups with less than 3 replicates
filtered_Pat.IDs <- which(freq <= 2)
removed_groups <- c()

for (i in 1:(length(filtered_Pat.IDs))) {
  
  end <- sum(freq[1:filtered_Pat.IDs[i]])
  start <- end - (freq[filtered_Pat.IDs[i]] - 1)
  
  removed_groups <- c(removed_groups, start:end)
  
}
IPSC_Generation_Data <- IPSC_Generation_Data[-removed_groups,]
dim(IPSC_Generation_Data)

Pat.ID. <- as.factor(Pat.ID.[-removed_groups])
length(Pat.ID.)
IPSC_Generation_Data$Pat.ID. <- Pat.ID.
IPSC_Generation_Data <- IPSC_Generation_Data[, c(14, 1:13)]


Subject_Demo_Data_filtered <- Subject_Demo_Data[Subject_Demo_Data$Pat.ID. %in% unique(IPSC_Generation_Data$Pat.ID.),]
Subject_Demo_Data_filtered <- Subject_Demo_Data_filtered[order(Subject_Demo_Data_filtered$Pat.ID.),]
dim(Subject_Demo_Data_filtered)

freq2 <- plyr::count(IPSC_Generation_Data$Pat.ID.)[, 2]
Subject_Demo_Data_filtered_rep <- Subject_Demo_Data_filtered[rep(1:nrow(Subject_Demo_Data_filtered), times = freq2), ]
Subject_Demo_Data_filtered_rep$Pat.ID. <- Pat.ID

Metadata <- cbind(IPSC_Generation_Data, Subject_Demo_Data_filtered_rep)
Metadata <- Metadata[, -15]
Metadata
####################################################################################################################################
#Generate DGE Object
DGE_object_og <- readDGE(files = files)
DGE_object <- DGE_object_og
DGE_object$samples
DGE_object$counts


#Sort and Add meta data to DGE object
file_info <- DGE_object@.Data[[1]]
CLONE.IDs <- regmatches(file_info$files, regexpr("_.*_", file_info$files)) 
CLONE.IDs <- sub(pattern = "_", replacement = "", x = CLONE.IDs)
CLONE.IDs <- sub(pattern = "_$", replacement = "", x = CLONE.IDs)
file_info <- file_info[gtools::mixedorder(CLONE.IDs), ]
file_info <- file_info[-removed_groups,]

file_info <- cbind(file_info, Metadata)
dim(file_info)
DGE_object@.Data[[1]] <- file_info
DGE_object$samples
# file_info$Pat.ID. <- Metadata$Pat.ID.
# file_info$CLONE.ID <- Metadata$CLONE.ID
DGE_object@.Data[[1]] <- file_info
DGE_object$samples$lib.size

#Sort counts matrix columns so same order as metadata rows
file_names <- colnames(DGE_object@.Data[[2]])
CLONE.IDs <- regmatches(file_names, regexpr("_.*_", file_names)) 
CLONE.IDs <- sub(pattern = "_", replacement = "", x = CLONE.IDs)
CLONE.IDs <- sub(pattern = "_$", replacement = "", x = CLONE.IDs)
DGE_object@.Data[[2]] <- DGE_object@.Data[[2]][,gtools::mixedorder(CLONE.IDs)]
DGE_object@.Data[[2]] <- DGE_object@.Data[[2]][,-removed_groups]
dim(DGE_object@.Data[[2]] )

########
#CPM Calulation
cpm <- cpm(DGE_object)
lcpm_unfiltered <- cpm(DGE_object, log=TRUE)

L <- mean(DGE_object$samples$lib.size) * 1e-6
M <- median(DGE_object$samples$lib.size) * 1e-6
c(L, M)

# cpm <- cpm[,gtools::mixedorder(CLONE.IDs)]

#Convert ensemble IDs to gene names
gene_info <- select(Homo.sapiens, keys = rownames(DGE_object@.Data[[2]]), columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")
gene_info <- gene_info[-which(duplicated(gene_info$ENSEMBL)),]
DGE_object$genes <- gene_info

#Removing Lowly Expressed Genes
ncol(DGE_object$counts)
table(rowSums(DGE_object$counts==0)==length(Pat.ID.))

# Pat.ID. <- DGE_object$samples$Pat.ID.
dim(DGE_object)
keep.exprs <- filterByExpr(DGE_object, group=Pat.ID.)
DGE_object <- DGE_object[keep.exprs,, keep.lib.sizes=FALSE]
dim(DGE_object)

lcpm <- lcpm_unfiltered 
lcpm.cutoff <- log2(10/M + 2/L)
# library(RColorBrewer)
nsamples <- ncol(DGE_object)
samplenames <- DGE_object$samples$CLONE.ID

col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

cpm <- cpm(DGE_object)
lcpm <- cpm(DGE_object, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#Normalising gene expression distributions
# During the sample preparation or sequencing process, external factors that are not of biological interest can affect the expression of individual samples.
# For example, samples processed in the first batch of an experiment can have higher expression overall when compared to samples processed in a second batch.
# It is assumed that all samples should have a similar range and distribution of expression values.
# Normalisation is required to ensure that the expression distributions of each sample are similar across the entire experiment.
par(mfrow=c(1,1))
DGE_object_TMM <- calcNormFactors(DGE_object, method = "TMM")
hist(DGE_object_TMM$samples$norm.factors)
lcpm_TMM <- cpm(DGE_object_TMM, log=TRUE)
cpm_TMM <- cpm(DGE_object_TMM)

# par(mfrow=c(1,2))
# boxplot(lcpm, las=2, col=col, main="")
# title(main="A. Example: Unnormalised data", ylab="Log-cpm")
# DGE_object <- calcNormFactors(DGE_object)
# DGE_object$samples$norm.factors
# ## [1] 0.0577 6.0829 1.2202 1.1648 1.1966 1.0466 1.1505 1.2543 1.1090
# 
# lcpm <- cpm(DGE_object_TMM, log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="B. Example: Normalised data", ylab="Log-cpm")


##############Clustering
#MDS

nb.cols <- nlevels(Pat.ID.)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
levels(mycolors) <-  mycolors
mycolors <- as.character(mycolors)
MDS <- plotMDS(lcpm, labels=Pat.ID., col=mycolors, top = 500)
MDS

title(main="A. Sample groups")

plotMDS(lcpm_TMM, labels=Pat.ID., col=mycolors)
title(main="A. Sample groups")




#PCA
# install.packages("ggfortify")
# library(ggfortify)
# cpm_for_PCA <- as.data.frame(cpm)
# colnames(cpm_for_PCA) <- DGE_object$samples$CLONE.ID
# cpm_for_PCA <- t(cpm_for_PCA)
# 
# pca_cpm <- prcomp(cpm_for_PCA, center = TRUE,scale. = TRUE)
# 
# cpm_for_PCA <- cbind.data.frame(cpm_for_PCA, group = (DGE_object$samples$Pat.ID.))
# 
# autoplot(pca_lcpm, data = cpm_for_PCA, label = TRUE, label.size = 5, colour = "group")

####
lcpm_for_PCA <- as.data.frame(lcpm)
colnames(lcpm_for_PCA) <- DGE_object$samples$CLONE.ID
lcpm_for_PCA <- t(lcpm_for_PCA)

pca_lcpm <- prcomp(lcpm_for_PCA, center = TRUE, scale. = TRUE)

lcpm_for_PCA2 <- cbind.data.frame(lcpm_for_PCA, group = (DGE_object$samples$Pat.ID.))

autoplot(pca_lcpm, data = lcpm_for_PCA2, label = TRUE, label.size = 5, colour = "group")

# col.lane <- lane
# levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
# col.lane <- as.character(col.lane)

#using dist matrix
lcpm_dist <- pdist::pdist(lcpm_for_PCA, lcpm_for_PCA)
as.matrix(lcpm_dist)

fit <- cmdscale(as.matrix(lcpm_dist), eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]

plot(fit$points)

autoplot(isoMDS(as.matrix(lcpm_dist)), colour = 'orange', size = 4, shape = 3)

autoplot(sammon(as.matrix(lcpm_dist)), colour = 'orange', label.size = 3)

##############
# rownames(cpm) <- gene_info$SYMBOL
# 
# counts <- cbind(gene_info, DGE_object@.Data[[2]])
# 
# cpm_filtered <- cpm[-which(duplicated(counts$SYMBOL)),]
# counts <- counts[-which(duplicated(counts$SYMBOL)),]
# 
# cpm_filtered <- cpm_filtered[-which(is.na(counts$SYMBOL)),]
# counts <- counts[-which(is.na(counts$SYMBOL)),]
# View(cpm_filtered)
# 
# colnames(cpm_filtered) <- Metadata$CLONE.ID
# 
# cpm_filtered["NANOG",]

cpm_TMM <- cpm(DGE_object_TMM)
cpm_TMM_mean <- apply(cpm_TMM, 1, FUN = mean)
cpm_TMM_sd <- apply(cpm_TMM, 1, FUN = sd)
cpm_TMM_var <- apply(cpm_TMM, 1, FUN = var)

df_mean_sd <- data.frame(mean = cpm_TMM_mean, sd = cpm_TMM_sd)
df_mean_var <- data.frame(mean = cpm_TMM_mean, var = cpm_TMM_var)


ggplot(data = df_mean_sd, aes(x = mean, y = sd)) + geom_point() + geom_smooth(method=lm)
ggplot(data = df_mean_var, aes(x = mean, y = var)) + geom_point() + xlim(c(0,8000)) + ylim(c(0,8000))

ggplot(data = df_mean_sd, aes(x = mean, y = sd)) + geom_point() + xlim(c(0,1500)) + ylim(c(0,800))

ggplot(data = df_mean_sd, aes(x = mean, y = sd)) + geom_point() + xlim(c(0,400)) + ylim(c(0,180))

View(sort(cpm_mean, decreasing = T))

cpm_sds <- apply(cpm_filtered, 1, FUN = sd)
View(sort(cpm_sds, decreasing = T))
sort(cpm_sds, decreasing = F)[10000]

cpm_sds["COX1"]

cpm_filtered[1:3, 315:317]
############################################################################################################
#generate means matrix 368_2
cpm_means_matrix <- matrix(nrow = nrow(cpm_filtered), ncol = length(freq))

rownames(cpm_means_matrix) <- rownames(cpm_filtered)
colnames(cpm_means_matrix) <- unique(Metadata$Pat.ID.)


start <- 1

for (i in 1:length(freq)) {
  
  end <- sum(freq[1:i])
  
  
  for (j in 1:nrow(cpm_filtered)) {
    
    
    cpm_means_matrix[j,i] <- mean(cpm_filtered[j, c(start:end)])
    
  }
  print(c(start:end))
  
  start <- end + 1
  
}
cpm_means_matrix <- cpm_means_matrix[, -which(freq <= 2)]
dim(cpm_means_matrix)
#generate sds matrix
cpm_sds_matrix <- matrix(nrow = nrow(cpm_filtered), ncol = length(freq))

rownames(cpm_sds_matrix) <- rownames(cpm_filtered)
colnames(cpm_sds_matrix) <- unique(Metadata$Pat.ID.)


start <- 1

for (i in 1:length(freq)) {
  
  end <- sum(freq[1:i])
  
  
  for (j in 1:nrow(cpm_filtered)) {
    
    
    cpm_sds_matrix[j,i] <- sd(c(cpm_filtered[j, c(start:end)]))
    
  }
  print(c(start:end))
  
  start <- end + 1
  
}
cpm_sds_matrix <- cpm_sds_matrix[, -which(freq <= 2)]

############################################################################################################
#generate_bar_plot
generate_bar_plot <- function(gene, metric = "sd" ) {
  
  gene.df <- data.frame(Pat.ID = colnames(cpm_sds_matrix), means = cpm_means_matrix[gene,], sds = cpm_sds_matrix[gene,], CVs = cpm_sds_matrix[gene,]/cpm_means_matrix[gene,])
  
  if (metric == "mean") {
    
    plot <- ggplot(data = gene.df, aes(x = reorder(Pat.ID, means), y = means)) + geom_bar(stat = "identity") + ggtitle(gene) +
      xlab("Pat.ID")
    
  } else if (metric == "sd") {
    
    plot <- ggplot(data = gene.df, aes(x = reorder(Pat.ID, sds), y = sds)) + geom_bar(stat = "identity") + ggtitle(gene) +
      xlab("Pat.ID")
    
  } else if (metric == "cv") {
    
    plot <- ggplot(data = gene.df, aes(x = reorder(Pat.ID, CVs), y = CVs)) + geom_bar(stat = "identity") + ggtitle(gene) +
      xlab("Pat.ID")
    
  }
  
  plot
}

generate_bar_plot(gene = "MEG3", metric = "sd")
