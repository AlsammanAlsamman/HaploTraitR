haplotype_clusters
LDsInfo

# read one cluster
ld_folder<-LDsInfo$ld_folder
i<-1
out_info<-LDsInfo$out_info
chr<-out_info[i,][1]
cls<-out_info[i,][2]
ld_matrix<-read.csv(file.path(ld_folder, paste(chr, cls, "ld_matrix.csv", sep="_")), row.names=1)
getwd()
write.csv(ld_matrix, "ld_matrix.csv")

# plot the LD matrix
library(ggplot2)
library(reshape2)

ld_matrix_mod<-ld_matrix
# remove upper triangle
ld_matrix_mod[upper.tri(ld_matrix_mod)]<-0
# metl to cr
# conver the matrix to a long format
ld_matrix_melt<-as.matrix(ld_matrix)
colnames(ld_matrix_melt)<-rownames(ld_matrix)
# melt the matrix to be row and column, and value
ld_matrix_melt<-melt(ld_matrix_melt)
colnames(ld_matrix_melt)<-c("snp1", "snp2", "value")
# add columns for sig_1 and sig_2
ld_matrix_melt$sig_1<-"black"
ld_matrix_melt$sig_2<-"black"
# if snp1 in gwas rsids, change sig_1 to red
ld_matrix_melt$sig_1[ld_matrix_melt$snp1 %in% gwas$rsid]<-"red"
# if snp2 in gwas rsids, change sig_2 to red
ld_matrix_melt$sig_2[ld_matrix_melt$snp2 %in% gwas$rsid]<-"red"
ld_matrix_melt$value<-as.numeric(ld_matrix_melt$value)
# plot the matrix
ggplot(ld_matrix_melt, aes(x=snp1, y=snp2, fill=value))+
  geom_tile()+
  scale_fill_gradient(low="white", high="blue")+
  geom_text(aes(label=round(value, 2)), size=3)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_color_manual(values=c("black", "red"))+
  theme(legend.position = "none")+
  geom_tile(aes(color=sig_1))+
  geom_tile(aes(color=sig_2))






snps<-rownames(ld_matrix)
snps.pos<-as.numeric(gsub(".*:", "", snps))
# convert to scale from 0 to 1
snps.pos<-(snps.pos-min(snps.pos))/(max(snps.pos)-min(snps.pos))
snps.labels<-snps
# create a data frame for the snps
snps.data<-data.frame(x=snps.pos, y=rep(1, length(snps.pos)), label=snps.labels)
snps.data
# create aggplot line and add the snps as ticks and text
library(ggplot2)
library(ggrepel)


library(ggplot2)
library(ggrepel)

ggplot() +
  # Add points to the plot
  geom_point(data = snps.data, aes(x = x, y = y), color = "black") +
  # Add vertical annotation labels with `geom_label_repel`
  geom_label_repel(
    data = snps.data,
    aes(x = x, y = y, label = label),
    size = 2,

    point.padding = 0.5,
#    force = 100,
    segment.size = 0.2,
    segment.color = "grey50",
    direction = "y",    # Align labels vertically above points
    angle = 90,         # Rotate text vertically
    nudge_y = 0.1,       # Offset labels slightly above points
    # dashed line
segment.inflect = TRUE,
hjust             = 0,
segment.curvature = -0.2,
force             = 0.5

  ) +
  # Clean up background, y-axis labels, and ticks
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),      # Remove y-axis labels
    axis.ticks = element_blank(),       # Remove all axis ticks
    panel.grid = element_blank(),       # Remove grid lines
    panel.background = element_blank(), # Remove background color
    legend.position = "none",
    # remove y axis title
    axis.title.y = element_blank()
  )
