setwd("D:/Work (Yr 2 Sem 1)/Thesis/")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/Liver/"

source("Scripts/helperFunctions.R")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)

#### READ IN SPATIAL DATA ####
seurat_obj_visium <- Load10X_Spatial("Data/Liver/JBO01")

colors <- c('#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6', '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87', '#5A0007', '#809693', '#6A3A4C', '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA', '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018')

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

#### PREPROCESS RESULTS #####
celltypes <- unique(read.csv(paste0(path, "anndatafineAnnot_celltypes.csv"))$fine_annot)
celltypes <- celltypes[!grepl("Capsular Fibroblasts", celltypes)]
celltypes <- str_replace(celltypes, "NaÃ¯ve", "Naïve")
celltypes <- sort(celltypes, method="radix")
celltypes[29:30] <- rev(celltypes[29:30])
ori_celltypes <- celltypes
ori_celltypes <- str_replace(ori_celltypes, "Portain", "Portal")
celltypes <- str_replace_all(celltypes, " ", ".")
celltypes <- str_replace(celltypes, "\\+", ".")

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

###### Let's evaluate lotte and my results
lotte_deconv <- readRDS(paste0(path, "results/cell2location/lotte_liver_deconv_cell2location_both.rds"))
my_deconv <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_both_pp.rds"))
lotte_deconv <- readRDS(paste0(path, "results/cell2location/liver_deconv_cell2location_allslides_both_pp.rds"))

keep_rows <- rownames(lotte_deconv)[rownames(lotte_deconv) %in% rownames(my_deconv)]
lotte_deconv <- lotte_deconv[keep_rows, ]
my_deconv <- my_deconv[keep_rows, ]
lotte_deconv <- lotte_deconv/rowSums(lotte_deconv)
my_deconv <- my_deconv/rowSums(my_deconv)

mean(diag(cor(t(my_deconv), t(lotte_deconv), method="spearman")))
diag(cor(my_deconv, lotte_deconv))
corr_celltypes <- round(diag(cor(my_deconv, lotte_deconv)), 2)

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

############## RCTD_BOTH AND LOTTE/MINE ################
mine = T
text <- ifelse(mine, "liver_deconv_cell2location_both_pp.rds", "lotte_liver_deconv_cell2location_both.rds")
c2l_deconv <- readRDS(paste0(path, "results/cell2location/", text))
RCTD_deconv_both <- readRDS(paste0(path, "results/RCTD/liver_deconv_RCTD_both_pp.rds"))

keep_rows <- rownames(c2l_deconv)[rownames(c2l_deconv) %in% rownames(RCTD_deconv_both)]
c2l_deconv <- c2l_deconv[keep_rows, ]
RCTD_deconv_both <- RCTD_deconv_both[keep_rows, ]
c2l_deconv <- c2l_deconv/rowSums(c2l_deconv)

# Evaluate
diag(cor(c2l_deconv, RCTD_deconv_both))
mean(diag(cor(t(RCTD_deconv_both), t(c2l_deconv), method="spearman")))
corr_celltypes <- round(diag(cor(c2l_deconv, RCTD_deconv_both)), 2)

#df_text <- data.frame(
#  label = round(diag(cor(c2l_deconv, RCTD_deconv_both)), 2),
#  ind   = ori_celltypes
#)
#print(p + geom_text(data=df_text, aes(x=-Inf, y=-Inf, label=label), hjust=-5, vjust=-12))

df <- data.frame(c2l=stack(c2l_deconv),
                 rctd=stack(RCTD_deconv_both),
                 ind=rep(ori_celltypes, each=nrow(c2l_deconv)),
                 ind_corr=rep(paste0(ori_celltypes, " (", corr_celltypes, ")"),
                              each=nrow(c2l_deconv)))

p <- ggplot(df, aes(x=c2l.values, y=rctd.values)) + geom_point() +
  facet_wrap(facets=vars(ind_corr), scales="free") +
  ylab('RCTD') + xlab('cell2location')
png("plots/liver_celltypecorr_RCTDboth.png", width=1500, height=1000)
print(p)
dev.off()

##### BARPLOT #####
# First make a pie chart
method_no <- 2 #1 is c2l, 2 is RCTD
deconvs <- list("cell2location"=c2l_deconv, "RCTD"=RCTD_deconv_both)
mean_props <- apply(deconvs[[method_no]], 2, mean)
hepa_prop <- mean_props[which(names(mean_props) == "Hepatocytes")]

piechart <- data.frame(group=c("Hepatocytes", "Others"), value=c(hepa_prop, 1-hepa_prop))
png(paste0("Data/Liver/plots/pie_chart_", names(deconvs)[method_no], ".png"),
    width=210, height=150, units="mm", res=200)
ggplot(piechart, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values=c(colors[which(new_cell_types == "Hepatocytes")], "gray")) + 
  coord_polar("y", start=0) + theme_void() + labs(fill = "Cell type") +
  ggtitle(names(deconvs)[method_no])
dev.off()



# Plot proportions excluding hepatocytes
df <- data.frame(reshape2::melt(lapply(deconvs, function(k) apply(k, 2, mean))),
                 df_celltypes = celltypes, df_ori_celltypes = ori_celltypes,
                 df_gen_celltypes = sapply(ori_celltypes, get_gen_annot, USE.NAMES=FALSE))
hepa_loc <- which(celltypes=="Hepatocytes")
df <- df[-which(df$df_celltypes=="Hepatocytes"),]
ggplot(df, aes(fill=factor(df_celltypes, levels=celltypes), y=value, x=L1)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual(values=colors[-hepa_loc], labels=ori_celltypes[-hepa_loc]) +
  ylab("Relative proportions") + xlab("Method") +labs(fill="Cell types")

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
ggplot(gen_df, aes(fill=factor(gen_df_celltypes, levels=gen_celltypes), y=value, x=L1)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual(values=gen_colors[-hepa_loc]) +
  ylab("Relative proportions") + xlab("Method") +labs(fill="Cell types")

##### SPATIAL SCATTERPIE #####
# Process data into right format
seurat_obj_visium2 <- seurat_obj_visium[,rownames(c2l_deconv)]
c2l_temp <- c2l_deconv
c2l_temp$barcodes <- rownames(c2l_temp)

seurat_obj_visium2@meta.data <- seurat_obj_visium2@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(c2l_temp, by = "barcodes") %>%
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

# Downsample spatial coord object for testing
set.seed(10)
pct_sample <- 0.05
spot_index <- sample.int(ncol(seurat_obj_visium2), ncol(seurat_obj_visium2)*pct_sample)
spatial_coord2 <- spatial_coord[spot_index,]

## Plot spatial scatterpie plot
png("Data/Liver/plots/scatterpie_fullscale_legend3.png", width=210, height=150, units="mm", res=300)
ggplot() + scatterpie::geom_scatterpie(data = spatial_coord2, aes(x = imagecol, y = imagerow),
                              cols = celltypes, color = NA, alpha = 1, pie_scale = 0.5) +
  scale_y_reverse() +  theme_half_open(11, rel_small = 1) +
  scale_fill_manual(values=colors, labels=ori_celltypes) +
  theme_void() + theme(legend.direction="vertical", legend.position="right")+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
# + theme(legend.position = "none") # Remove the legend
dev.off()

#### DOT PLOT? #####
cluster_assignments <- read.csv("Data/Liver/clusters_list.csv", row.names=1)
cluster_assignments <- cluster_assignments[match(rownames(c2l_deconv), rownames(cluster_assignments)),,drop=FALSE]

df <- data.frame(c2l=stack(c2l_deconv),
                 cluster=cluster_assignments)

# First find the mean within each cluster and celltype, then find the standard deviation (CV) between clusters
sd_between_clusters <- df %>% group_by(c2l.ind, clusters) %>% summarise(mean_prop=mean(c2l.values)) %>%
  group_by(c2l.ind) %>% summarise(sd_cluster=sd(mean_prop)/mean(mean_prop))
sd_between_clusters <- sd_between_clusters[order(sd_between_clusters$sd_cluster, decreasing=T),]
top_10_variable <- sd_between_clusters$c2l.ind[1:10]

df <- df[df$c2l.ind %in% top_10_variable,]
ggplot(df, aes(x=factor(clusters), y=factor(c2l.ind, levels=rev(top_10_variable)), size=c2l.values)) + geom_point()

##### RCTD SC #####
mine = F
text <- ifelse(mine, "liver_deconv_cell2location_sc_pp.rds", "lotte_liver_deconv_cell2location_sc.rds")
c2l_deconv_sc <- readRDS(paste0(path, "results/cell2location/", text))
RCTD_deconv_sc <- readRDS(paste0(path, "results/RCTD/liver_deconv_RCTD_sc_pp.rds"))

celltype_sc <- celltypes[-to_remove]
ori_celltypes_sc <- ori_celltypes[-to_remove]
sum(is.na(match(celltype_sc, colnames(RCTD_deconv_sc))))

keep_rows <- rownames(c2l_deconv_sc)[rownames(c2l_deconv_sc) %in% rownames(RCTD_deconv_sc)]
c2l_deconv_sc <- c2l_deconv_sc[keep_rows, ]
RCTD_deconv_sc <- RCTD_deconv_sc[keep_rows, ]
c2l_deconv_sc <- c2l_deconv_sc/rowSums(c2l_deconv_sc)

# Evaluate
diag(cor(c2l_deconv_sc, RCTD_deconv_sc))
mean(diag(cor(t(RCTD_deconv_sc), t(c2l_deconv_sc), method="spearman")))

corr_celltypes_sc <- round(diag(cor(c2l_deconv_sc, RCTD_deconv_sc)), 2)

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
