library(RcppCNPy) # for reading python numpy (.npy) files

setwd("~/Documents/GitHub/Ecological_Genomics_23/PopGenomics/results/")

list.files()

s<-npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))

p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

names(pval) = c("p_PC1","p_PC2")

cutoff=1e-3   

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]
outliers_PC2 <- p_filtered[which(pval$p_PC2<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff for each PC axis?
dim(outliers_PC1)[1]
dim(outliers_PC2)[1]

# write them out to a file
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)

clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T)

fviz_eig(clim_PCA)

fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

dimdesc(clim_PCA)[1:2]

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio12"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)

