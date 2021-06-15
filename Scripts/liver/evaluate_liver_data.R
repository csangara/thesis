#### EVALUATION OF LIVER DATA ####
setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/Liver/"

source("Scripts/helperFunctions.R")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)

#### 1. INITIALIZING AND PREPROCESSING DATA ####
# Read in spatial data
seurat_obj_visium <- Load10X_Spatial("Data/Liver/JBO01")

# Coarse annotation & transforming the fine to coarse annotation
gen_celltypes <- c('B cells', 'Basophils', 'Cholangiocytes', 'Endothelial cells', 'Fibroblasts', 'Hepatocytes', 'HsPCs', 'ILC1s', 'Kupffer cells',
                   'Mig. cDCs', 'Monocytes & Monocyte-derived cells', 'NK cells', 'Neutrophils', 'T cells', 'cDC1s', 'cDC2s', 'pDCs')
get_gen_annot <- function(celltype){
  conditions <- c(grepl("Endothelial|LSECs", celltype),
                  grepl("Stellate|Mesothelial|Fibroblast", celltype),
                  grepl("Monocyte|Mac", celltype),
                  grepl("CD4|CD8|NKT|Th1|CTLs|TRegs", celltype))
  replacements <- c('Endothelial cells', 'Fibroblasts', 'Monocytes & Monocyte-derived cells', 'T cells')
  
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}

# Colors (exported from python)
colors <- c('#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6',
            '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87', '#5A0007', '#809693', '#6A3A4C', '#1B4400',
            '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900', '#00C2A0', '#FFAA92',
            '#FF90C9', '#B903AA', '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018')

# Get cell types and sort them so they match the order from Scanpy
celltypes <- unique(read.csv(paste0(path, "anndatafineAnnot_celltypes.csv"))$fine_annot)
celltypes <- celltypes[!grepl("Capsular Fibroblasts", celltypes)]
celltypes <- str_replace(celltypes, "NaÃ¯ve", "Naïve")
celltypes <- sort(celltypes, method="radix")
celltypes[29:30] <- rev(celltypes[29:30])
ori_celltypes <- celltypes # Proper names of cell types
ori_celltypes <- str_replace(ori_celltypes, "Portain", "Portal")
celltypes <- str_replace_all(celltypes, " ", ".")
celltypes <- str_replace(celltypes, "\\+", ".")

# Read cell2location results (sc is for model build on using scRNA-seq ref)
# sc = T
# text = ifelse(sc, "_sc", "_both")
# my_c2l <- read.csv(paste0(path, "results/cell2location/liver_deconv_cell2location", text, ".csv"),
#                    row.names=1, quote="")
# colnames(my_c2l) <- str_replace(colnames(my_c2l), "q05_spot_factors", "")
# colnames(my_c2l) <- str_replace(colnames(my_c2l), "NaÃ.ve", "Naïve")
# if (sc){  temp_celltypes <- celltypes[-which(!celltypes %in% colnames(my_c2l))]
# } else {  temp_celltypes <- celltypes }
# sum(is.na(match(temp_celltypes, colnames(my_c2l))))
# my_c2l <- my_c2l[,match(temp_celltypes, colnames(my_c2l))]
# saveRDS(my_c2l, paste0(path, "results/cell2location/liver_deconv_cell2location", text, "_pp.rds"))

# Read RCTD results
# sc = F
# text = ifelse(sc, "_sc", "_both")
# RCTD_deconv <- readRDS(paste0(path, "results/RCTD/liver_deconv_RCTD", text, ".rds"))
# RCTD_deconv <- data.frame(RCTD_deconv)
# colnames(RCTD_deconv) <- str_replace_all(colnames(RCTD_deconv), " ", ".")
# colnames(RCTD_deconv) <- str_replace(colnames(RCTD_deconv), "\\+", ".")
# if (sc){  temp_celltypes <- celltypes[-which(!celltypes %in% colnames(RCTD_deconv))]
# } else {  temp_celltypes <- celltypes }
# sum(is.na(match(temp_celltypes, colnames(RCTD_deconv))))
# RCTD_deconv <- RCTD_deconv[,match(temp_celltypes, colnames(RCTD_deconv))]
# saveRDS(RCTD_deconv, paste0(path, "results/RCTD/liver_deconv_RCTD", text, "_pp.rds"))


#### 2. COMPARING RCTD AND CELL2LOCATION RESULTS ####
c2l_deconv <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_both_pp.rds"))
RCTD_deconv_both <- readRDS(paste0(path, "results/RCTD/liver_deconv_RCTD_both_pp.rds"))
keep_rows <- rownames(c2l_deconv)[rownames(c2l_deconv) %in% rownames(RCTD_deconv_both)]
c2l_deconv <- c2l_deconv[keep_rows, ]
RCTD_deconv_both <- RCTD_deconv_both[keep_rows, ]
c2l_deconv <- c2l_deconv/rowSums(c2l_deconv)

# Correlation
# diag(cor(c2l_deconv, RCTD_deconv_both))
# mean(diag(cor(t(RCTD_deconv_both), t(c2l_deconv), method="spearman")))

## Plot c2l vs RCTD proportions ##
corr_celltypes <- round(diag(cor(c2l_deconv, RCTD_deconv_both)), 2)
df <- data.frame(c2l=stack(c2l_deconv),
                 rctd=stack(RCTD_deconv_both),
                 ind=rep(ori_celltypes, each=nrow(c2l_deconv)),
                 ind_corr=rep(paste0(ori_celltypes, " (", corr_celltypes, ")"),
                              each=nrow(c2l_deconv)))
df <- mutate(df, ind_corr = str_wrap(ind_corr, width = 25))
df$ind_corr <- factor(df$ind_corr, levels=unique(df$ind_corr))

p <- ggplot(df, aes(x=c2l.values, y=rctd.values)) + geom_point() +
  facet_wrap(facets=vars(ind_corr), scales="free") +
  ylab('RCTD') + xlab('cell2location')
png("Data/Liver/plots/liver_celltypecorr_RCTDboth2.png", width=297, height=210, units="mm", res=200)
print(p)
dev.off()

## Pie chart (hepatocytes vs others) ##
deconvs <- list("cell2location"=c2l_deconv, "RCTD"=RCTD_deconv_both)
for (i in 1:2){ #1 is c2l, 2 is RCTD
  mean_props <- colMeans(deconvs[[i]])
  hepa_prop <- mean_props[which(names(mean_props) == "Hepatocytes")]
  piechart <- data.frame(group=c("Hepatocytes", "Others"), value=c(hepa_prop, 1-hepa_prop))
  
  png(paste0("Data/Liver/plots/pie_chart_", names(deconvs)[i], ".png"),
      width=210, height=150, units="mm", res=200)
  ggplot(piechart, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    scale_fill_manual(values=c(colors[which(new_cell_types == "Hepatocytes")], "gray")) + 
    coord_polar("y", start=0) + theme_void() + labs(fill = "Cell type") +
    ggtitle(names(deconvs)[i])
  dev.off()
}

## Bar plot (excluding hepatocytes) ##
# Fine annotation #
deconvs <- list("cell2location"=c2l_deconv, "RCTD"=RCTD_deconv_both)
df <- data.frame(reshape2::melt(lapply(deconvs, colMeans)),
                 df_celltypes = celltypes, df_ori_celltypes = ori_celltypes,
                 df_gen_celltypes = sapply(ori_celltypes, get_gen_annot, USE.NAMES=FALSE))
hepa_loc <- which(celltypes=="Hepatocytes")
df <- df[-which(df$df_celltypes=="Hepatocytes"),]
p <- ggplot(df, aes(fill=factor(df_celltypes, levels=celltypes), y=value, x=L1)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual(values=colors[-hepa_loc], labels=ori_celltypes[-hepa_loc]) +
  ylab("Relative proportions") + xlab("Method") +labs(fill="Cell types") +
  theme_classic()
png("Data/Liver/plots/barplot_fineannot.png", width=175, height=110, units="mm", res=200)
print(p)
dev.off()

# Coarse annotation #
# Sum up all cell types belonging to the same coarse annotation
gen_deconvs <- lapply(deconvs, function(deconv) {
  gen_deconv <- deconv %>% `colnames<-`(sapply(colnames(deconv), get_gen_annot))
  gen_deconv <- t(rowsum(t(gen_deconv), group=colnames(gen_deconv)))
  gen_deconv <- gen_deconv[, sort(colnames(gen_deconv), method="radix")]
  return (gen_deconv)
})

gen_colors <- colors[c(1, 2, 6, 12, 27, 8, 9, 10, 11, 15, 18, 19, 23, 29, 33:35)]
gen_df <- data.frame(reshape2::melt(lapply(gen_deconvs, function(k) apply(k, 2, mean))),
                 gen_df_celltypes = gen_celltypes)
hepa_loc <- which(gen_celltypes=="Hepatocytes")
gen_df <- gen_df[-which(gen_df$gen_df_celltypes=="Hepatocytes"),]
p <- ggplot(gen_df, aes(fill=factor(gen_df_celltypes, levels=gen_celltypes), y=value, x=L1)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual(values=gen_colors[-hepa_loc]) +
  ylab("Relative proportions") + xlab("Method") +labs(fill="Cell types") +
  theme_classic()
png("Data/Liver/plots/barplot_coarseannot.png", width=150, height=100, units="mm", res=200)
print(p)
dev.off()

#### 3. DOWNSTREAM ANALYSIS #####
deconv_obj <- RCTD_deconv_both # c2l_deconv or RCTD_deconv_both
seurat_obj_visium2 <- seurat_obj_visium[,rownames(deconv_obj)]
cluster_assignments <- read.csv("Data/Liver/clusters_list.csv", row.names=1)
cluster_assignments <- cluster_assignments[match(rownames(deconv_obj), rownames(cluster_assignments)),,drop=FALSE]

## Moran's I ##
# very similar results whether you use "row" and "col" or "imagerow" and "imagecol"
coords_dist <- as.matrix(dist(seurat_obj_visium2@images$slice1@coordinates[c("imagerow", "imagecol")]))
coords_dist <- 1/coords_dist
diag(coords_dist) <- 0
all_morans <- apply(deconv_obj, 2, ape::Moran.I, coords_dist)
p_morans <- reshape2::melt(all_morans) %>% tidyr::spread(., L2, value)
p_morans$padj <- p.adjust(p_morans$p.value)
p_morans[order(p_morans$observed, decreasing=TRUE),][1:5,] # Order by most spatially autocorrelated

## Correlation test ##
deconv_obj_cor <- cor(deconv_obj)
# deconv_obj_cor <- cor(gen_deconvs[["cell2location"]]) # if you just want the coarse annot
p.mat <- corrplot::cor.mtest(mat = deconv_obj_cor, conf.level = 0.95) # Test
pmat_adj <- matrix(p.adjust(p.mat$p), ncol=ncol(p.mat$p))

# How many significant?
(sum(p.mat$p < 0.05) - 35)/2
(sum(pmat_adj < 0.05) - 35)/2

# Writing out the significant correlations
indices <- which(pmat_adj< 0.05, arr.ind=TRUE)
psigs <- cbind(ori_celltypes[indices[,1]], ori_celltypes[indices[,2]],
               sapply(1:nrow(indices), function(k) deconv_obj_cor[indices[k,1], indices[k,2]]))
psigs <- data.frame(psigs) %>% filter(X1 != X2)
psigs <- psigs[order(as.numeric(psigs$X3), decreasing=TRUE),]
psigs <- psigs[seq(1, nrow(psigs), 2),] # Only get even rows, because correlation matrix is doubled
write.table(psigs,"Data/Liver/sigcorr.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# In case you want to remove rows and columns that are not sig, use these instead!
new_pmat <- pmat_adj %>% `diag<-`(1)
new_pmat <- new_pmat[rowSums(new_pmat < 0.05) > 0, colSums(new_pmat < 0.05) > 0]
new_cor <- deconv_obj_cor[rowSums(new_pmat < 0.05) > 0, colSums(new_pmat < 0.05) > 0]

col3 <- colorRampPalette(c(rep("red", 10), "white", rep("blue", 10)))
png("Data/Liver/plots/correlationplot_upper_5_adj.png", width=2000, height=1500)
corrplot::corrplot(
  deconv_obj_cor,
  p.mat = pmat_adj,
  type="lower", order="hclust", cl.pos="n", tl.col="gray",
  sig.level=.05, method="number", col=col3(20), insig="blank")
dev.off()

## Coefficient of variation ##
df <- data.frame(deconv=stack(deconv_obj),
                 cluster=cluster_assignments)

# Find the mean within each cluster and celltype, then divide std between clusters by mean
cv <- df %>% group_by(deconv.ind, clusters) %>% summarise(mean_prop=mean(deconv.values)) %>%
  group_by(deconv.ind) %>% summarise(sd_cluster=sd(mean_prop)/mean(mean_prop))
cv <- cv[order(cv$sd_cluster, decreasing=T),]
(top_10_cv <- cv$deconv.ind[1:10])

# If you want to plot it in a dot plot
df <- df[df$deconv.ind %in% top_10_cv,]
ggplot(df, aes(x=factor(clusters), y=factor(deconv.ind, levels=rev(top_10_cv)),
               size=deconv.values)) + geom_point()

## Spatial map of cell types ##
seurat_obj_visium2 <- seurat_obj_visium[,rownames(deconv_obj)]
celltypes2 <- str_replace(celltypes, "Portain", "Portal")
deconv_obj2 <- deconv_obj %>% `colnames<-`(celltypes2)

decon_df <- deconv_obj2 %>% data.frame() %>% tibble::rownames_to_column("barcodes")
seurat_obj_visium2@meta.data <- seurat_obj_visium2@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

# Plot celltypes of interest
coi <- celltypes2[grepl("Portal|Central|LSEC|Lymphatic|KC|Stellate|Cholang|Fibro", celltypes2)]
coi <- coi[c(3,2,6,7,1,5,4,8)]
p <- SpatialFeaturePlot(seurat_obj_visium2, features=coi,
                   ncol=3,  pt.size.factor = 2)
png("Data/Liver/plots/spatialmap_oi_big.png", width=250, height=297, units="mm", res=200)
print(p)
dev.off()

# Remaining cell types
filtered_celltypes <- celltypes2[!celltypes2 %in% coi]
for (i in 0:1){
  range = seq(i*15+1, 15*(i+1))
  p <- SpatialFeaturePlot(seurat_obj_visium2, features=filtered_celltypes[range],
                          ncol=3,  pt.size.factor = 2)
  png(paste0("Data/Liver/plots/spatialmap_big", i+1, ".png"), width=250, height=297, units="mm", res=200)
  print(p)
  dev.off()
}

## Heatmap ##
celltypes2 <- str_replace(celltypes, "Portain", "Portal")
df <- data.frame(deconv=stack(deconv_obj),
                 cluster=cluster_assignments)
df <- df[grepl("Portain|Central|LSEC|Lymphatic|KC|Stellate|Cholang|Fibro", df$deconv.ind),]
df_heatmap <- reshape2::dcast(df, deconv.ind~clusters, value.var="deconv.values", fun.aggregate=mean) %>%
  tibble::column_to_rownames("deconv.ind")
df_heatmap <- df_heatmap[c(3,2,6,7,1,5,4,8),]

png("Data/Liver/plots/heatmap.png", width=210, height=150, units="mm", res=200)
gplots::heatmap.2(as.matrix(df_heatmap), scale="row", Colv=FALSE, Rowv=FALSE, dendrogram="none")
#gplots::heatmap.2(as.matrix(df_heatmap), scale="row")
dev.off()

#### 4. OTHERS (CAN IGNORE) ####
## Plotting the clusters ##
seurat_obj_visium2$clusters <- cluster_assignments
SpatialPlot(seurat_obj_visium2, group.by="clusters", pt.size.factor=1.8)

## Spatial scatterpie (pie chart for each spot) ##
# Process data into right format
seurat_obj_visium2 <- seurat_obj_visium[,rownames(deconv_obj)]
deconv_obj_temp <- deconv_obj
deconv_obj_temp$barcodes <- rownames(deconv_obj_temp)

seurat_obj_visium2@meta.data <- seurat_obj_visium2@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(deconv_obj_temp, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

metadata_ds <- data.frame(seurat_obj_visium2@meta.data)
spatial_coord <- data.frame(seurat_obj_visium2@images[["slice1"]]@coordinates) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::mutate(imagerow_scaled = imagerow * seurat_obj_visium2@images[["slice1"]]@scale.factors$hires,
                imagecol_scaled = imagecol * seurat_obj_visium2@images[["slice1"]]@scale.factors$hires) %>%
  dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("ID"), by = "ID")

# If we want to visualize the image, but kind of a PITA so I did this in photoshop
# img <- png::readPNG(paste0(path, "JBO01/spatial/tissue_hires_image.png"))
# img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
# annotation_custom(grob = img_grob, xmin = 0, xmax = 2000, ymin = 0, ymax = -2000) + # Add this in ggplot

# If you wanna test it out with a few spots first
# set.seed(10)
# pct_sample <- 0.05
# spot_index <- sample.int(ncol(seurat_obj_visium2), ncol(seurat_obj_visium2)*pct_sample)
# spatial_coord2 <- spatial_coord[spot_index,]

## Plot spatial scatterpie plot
png("Data/Liver/plots/scatterpie_fullscale_legend3.png", width=210, height=150, units="mm", res=300)
ggplot() + scatterpie::geom_scatterpie(data = spatial_coord, aes(x = imagecol, y = imagerow),
                                       cols = celltypes, color = NA, alpha = 1, pie_scale = 0.5) +
  scale_y_reverse() +  theme_half_open(11, rel_small = 1) +
  scale_fill_manual(values=colors, labels=ori_celltypes) +
  theme_void() + theme(legend.direction="vertical", legend.position="right")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
# + theme(legend.position = "none") # Remove the legend
dev.off()

## Results of using only single-cell reference ##
c2l_deconv_sc <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_sc_pp.rds"))
RCTD_deconv_sc <- readRDS(paste0(path, "results/RCTD/liver_deconv_RCTD_sc_pp.rds"))

celltype_sc <- celltypes[-to_remove]
ori_celltypes_sc <- ori_celltypes[-to_remove]
sum(is.na(match(celltype_sc, colnames(RCTD_deconv_sc))))

keep_rows <- rownames(c2l_deconv_sc)[rownames(c2l_deconv_sc) %in% rownames(RCTD_deconv_sc)]
c2l_deconv_sc <- c2l_deconv_sc[keep_rows, ]
RCTD_deconv_sc <- RCTD_deconv_sc[keep_rows, ]
c2l_deconv_sc <- c2l_deconv_sc/rowSums(c2l_deconv_sc)
corr_celltypes_sc <- round(diag(cor(c2l_deconv_sc, RCTD_deconv_sc)), 2)

# Correlation
# diag(cor(c2l_deconv_sc, RCTD_deconv_sc))
# mean(diag(cor(t(RCTD_deconv_sc), t(c2l_deconv_sc), method="spearman")))

df <- data.frame(c2l=stack(c2l_deconv_sc),
                 rctd=stack(RCTD_deconv_sc),
                 ind=rep(ori_celltypes_sc, each=nrow(c2l_deconv_sc)),
                 ind_corr=rep(paste0(ori_celltypes_sc, " (", corr_celltypes_sc, ")"),
                              each=nrow(c2l_deconv_sc)))

png("Data/Liver/plots/liver_celltypecorr_RCTDsc_mine.png", width=1500, height=1000)
p <- ggplot(df, aes(x=c2l.values, y=rctd.values)) + geom_point() +
  facet_wrap(facets=vars(ind_corr), scales="free") +
  ylab('RCTD') + xlab('cell2location')
print(p)
dev.off()

## Comparing with another person's (Lotte) results ##
# Preprocessing Lotte's results
# c2l_deconv <- read.csv(paste0(path, "results/cell2location/lotte_liver_deconv_cell2location_both.csv"),
#                        row.names=1, quote="")
# rownames(c2l_deconv) <- str_replace(rownames(c2l_deconv), '"', "")
# c2l_deconv <- c2l_deconv[,grepl("q05_spot_factors", colnames(c2l_deconv))]
# colnames(c2l_deconv) <- str_replace(colnames(c2l_deconv), "q05_spot_factors", "")
# colnames(c2l_deconv) <- str_replace(colnames(c2l_deconv), "NaÃ.ve", "Naïve")
# sum(is.na(match(celltypes, colnames(c2l_deconv))))
# c2l_deconv <- c2l_deconv[,match(celltypes, colnames(c2l_deconv))]
# saveRDS(c2l_deconv, paste0(path, "results/cell2location/lotte_liver_deconv_cell2location_both.rds"))
# 
# to_remove <- which(grepl("HsPCs", celltypes))
# c2l_deconv_sc <- c2l_deconv[,-to_remove]
# saveRDS(c2l_deconv_sc, paste0(path, "results/cell2location/lotte_liver_deconv_cell2location_sc.rds"))

lotte_deconv <- readRDS(paste0(path, "results/cell2location/lotte_liver_deconv_cell2location_both.rds"))
# lotte_deconv <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_allslides_both_pp.rds"))
my_deconv <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_both_pp.rds"))

keep_rows <- rownames(lotte_deconv)[rownames(lotte_deconv) %in% rownames(my_deconv)]
lotte_deconv <- lotte_deconv[keep_rows, ]
my_deconv <- my_deconv[keep_rows, ]
lotte_deconv <- lotte_deconv/rowSums(lotte_deconv)
my_deconv <- my_deconv/rowSums(my_deconv)
corr_celltypes <- round(diag(cor(my_deconv, lotte_deconv)), 2)

# Correlation
# mean(diag(cor(t(my_deconv), t(lotte_deconv), method="spearman")))
# diag(cor(my_deconv, lotte_deconv))

df <- data.frame(my_c2l=stack(my_deconv),
                 lotte_c2l=stack(lotte_deconv),
                 ind=rep(ori_celltypes, each=nrow(lotte_deconv)),
                 ind_corr=rep(paste0(ori_celltypes, " (", corr_celltypes, ")"),
                              each=nrow(lotte_deconv)))

p <- ggplot(df, aes(x=my_c2l.values, y=lotte_c2l.values)) + geom_point() +
  facet_wrap(facets=vars(ind_corr), scales="free") +
  xlab('Mine') + ylab('Mine (allslides)')
png("Data/Liver/plots/liver_celltypecorr_mineallvsone.png", width=1500, height=1000)
print(p)
dev.off()