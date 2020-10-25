# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Homo.sapiens")



kilpinen_data <- read.table(file = "/Users/gordid/Desktop/iPSC/Pera/kilpinen/Stemformatics_dataset7254.gct", header = TRUE, skip = 2)
  dim(kilpinen_data)
  
  grep(pattern = "cups_|eofe|wetu|nekd|nusw", colnames(kilpinen_data) )
  kilpinen_data <- kilpinen_data[, -(grep(pattern = "cups_|eofe|wetu|nekd|nusw", colnames(kilpinen_data) ))]
  head(kilpinen_data)
which(duplicated(kilpinen_data[, 1]))

gene_names <- biomaRt::select(Homo.sapiens, keys = as.character(kilpinen_data[, 1]), columns = c("SYMBOL"), keytype = "ENSEMBL")

which(duplicated(gene_names$SYMBOL))
which(is.na(gene_names$SYMBOL))
kilpinen_data

symbols <- gene_names$SYMBOL[-which(duplicated(gene_names$ENSEMBL))]
which(duplicated(symbols))

kilpinen_data_2 <- kilpinen_data[, 3:64]
kilpinen_data_2 <- kilpinen_data_2[-(which(duplicated(symbols))), ]
symbols <- symbols[-which(duplicated(symbols))]

kilpinen_data_2 <- kilpinen_data_2[-(which(is.na(symbols))), ]
symbols <- symbols[-(which(is.na(symbols)))]

row.names(kilpinen_data_2) <- symbols


dim(kilpinen_data_2)
length(symbols)

# pheatmap::pheatmap(kilpinen_data_2)

grep(pattern = paste(c(toupper(Variables_AS), "NODAL", "GDF3", "LEFTY1"), collapse = "|"), x = row.names(kilpinen_data_2) )

row.names(kilpinen_data_2)[grep(pattern = paste(c(toupper(Variables_AS), "NODAL", "GDF3", "LEFTY1"), collapse = "|"), x = row.names(kilpinen_data_2) )]

pheatmap::pheatmap(kilpinen_data_2[grep(pattern = paste(c(toupper(Variables_AS), "NODAL", "GDF3", "LEFTY1"), collapse = "|"), x = row.names(kilpinen_data_2) ), ])

