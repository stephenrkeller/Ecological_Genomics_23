setwd("") # set your path to your results folder in your repo where you saved your diversity stats file

list.files() # list out the files in this folder to make sure you're in the right spot.

# First let's read in the diversity stats
theta <- read.table("_.thetas",sep="\t",header=T)

theta$tWsite = theta$tW/theta$nSites #scales the theta-W by the number of sites
theta$tPsite = theta$tP/theta$nSites #scales the theta-Pi by the number of sites

summary(theta)

# You can order the contig list to show you the contigs with the highest values of Tajima's D, or the lowest

head(theta[order(theta$Tajima, decreasing = TRUE),]) # top 10 Tajima's D values

head(theta[order(theta$Tajima, decreasing = FALSE),]) # bottom 10 Tajima's D values

#You can also look for contigs that have combinations of high Tajima's D and low diversity -- these could represent outliers for selection
#theta[which(theta$Tajima>1.5 & theta$tPsite<0.001),]


sfs<-scan('9999_.sfs')
sfs<-sfs[-c(1,which(sfs==0))]
sfs<-sfs/sum(sfs)

# Be sure to replace "9999" with your pop code in the "main" legend below
barplot(sfs,xlab="Chromosomes",
        names=1:length(sfs),
        ylab="Proportions",
        main="Pop 9999 Site Frequency Spectrum",
        col='blue')

# Put the nucleotide diversities, Tajima's D, and SFS into a 4-panel figure
par(mfrow=c(2,2))
hist(theta$tWsite, xlab="theta-W", main="Watterson's theta")
hist(theta$tPsite, xlab="theta-Pi", main="Pairwise Nucleotide Diversity")
hist(theta$Tajima, xlab="D", main="Tajima's D")
barplot(sfs,names=1:length(sfs),main='Site Frequency Spectrum')


# To reset the panel plotting, execute the line below:
dev.off()

# Get the total number of sites sequenced:

sum(theta$nSites)

# plot the 2D SFS (for ex, with black vs. red spruce

###############################################################################
# below code taken from https://github.com/isinaltinkaya/WPSG2022/blob/main/fst_pbs.md


##run in R                      
yc<-scan("yri.ceu.ml")
yj<-scan("yri.jpt.ml")
jc<-scan("jpt.ceu.ml")
source("plot2dSFS.R")
plot2<-function(s,...){
  dim(s)<-c(21,21)
  s[1]<-NA
  s[21,21]<-NA
  s<-s/sum(s,na.rm=T)
  
  pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")
  pplot(s/sum(s,na.rm=T),pal=pal,...)
}

plot2(yc,ylab="YRI",xlab="CEU")
x11()
plot2(yj,ylab="YRI",xlab="JPT")
x11()
plot2(jc,ylab="JPT",xlab="CEU")

)

# Let's see how the Fst and PBS varies between different regions of the genome my using a sliding windows approach (windows site of 50kb)


# $REALSFS fst index yri.saf.idx jpt.saf.idx ceu.saf.idx -fstout yri.jpt.ceu -sfs yri.jpt.ml -sfs yri.ceu.ml -sfs jpt.ceu.ml
# $REALSFS fst stats2 yri.jpt.ceu.fst.idx -win 50000 -step 10000 >slidingwindowBackground


r<-read.delim("slidingwindowBackground",as.is=T,head=T)
names(r)[-c(1:4)] <- c("wFst_YRI_JPT","wFst_YRI_CEU","wFst_JPT_CEU","PBS_YRI","PBS_JPT","PBS_CEU")

head(r) #print the results to the screen

#plot the distribution of Fst
mmax<-max(c(r$wFst_YRI_JPT,r$wFst_YRI_CEU,r$wFst_JPT_CEU),na.rm=T)
par(mfcol=c(3,2))
hist(r$wFst_YRI_JPT,col="lavender",xlim=c(0,mmax),br=20)
hist(r$wFst_YRI_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$wFst_JPT_CEU,col="hotpink",xlim=c(0,mmax),br=20)

mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)

#plot the distribution of PBS
mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)
hist(r$PBS_YRI,col="lavender",xlim=c(0,mmax),br=20)
hist(r$PBS_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$PBS_JPT,col="hotpink",xlim=c(0,mmax),br=20)


