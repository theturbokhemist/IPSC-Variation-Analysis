setwd("/Users/gordid/Desktop/iPSC/Pera/Variation/GSE79636_RAW")

DGE_variation <- readDGE(list.files(), columns = c(1, 2))

DGE_variation_original <- DGE_variation
DGE_variation <- DGE_variation_original
###
cpm <- cpm(DGE_variation)
lcpm <- cpm(DGE_variation, log=TRUE)
###

Ens_IDs <- rownames(DGE_variation$counts)

Gene_names <- biomaRt::select(Homo.sapiens, keys = Ens_IDs, columns = c("SYMBOL", "TXCHROM"), keytype = "ENSEMBL")

Gene_names1 <- Gene_names[-which(duplicated(Gene_names$ENSEMBL)),]

rownames(cpm) <- Gene_names1$SYMBOL

cpm <- cpm[(!is.na(Gene_names1$SYMBOL)),]
Gene_names2 <- Gene_names1[(!is.na(Gene_names1$SYMBOL)),]

cpm <- cpm[-(which(duplicated(Gene_names2$SYMBOL))),]
Gene_names3 <-Gene_names2[-(which(duplicated(Gene_names2$SYMBOL))),]
###

cpm_sds <- apply(cpm, 1, FUN = sd)
sort(cpm_sds, decreasing = T)

cpm_sds["NANOG"]

