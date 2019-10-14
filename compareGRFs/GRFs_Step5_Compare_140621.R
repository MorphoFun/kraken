################# GRF DATA COMPARISON #########################################
# Running analyses to determine how the GRF data differs between the forelimb and pectoral fin
# Af  = Ambystoma tigrinum salamander (had to use Af for Ambystoma forelimb since Megan used At)
# Pb = Periophthalmus barbarus mudskipper
# Pw = Pleurodeles waltl newt

# 10-13-11: added code to evaluate interspecific differences of the net GRF magnitude at peak net GRF

# Establishing the location of the folder that stores the data file I'll be looking at and manipulating

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")

# Peak Net GRF data
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Peak Net GRF')
PeakNetPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Peak Net GRF')
PeakNetFiles <- list.files(PeakNetPath, pattern=".csv", full=TRUE)

# Getting data
Pec.Af.PeakNet.Use <- PeakNetFiles[which(substring(PeakNetFiles, 114, 131)=='Pec_Af_PeakNetGRF_')]
Pec.Af.PeakNet <- data.frame(read.csv(Pec.Af.PeakNet.Use, header=TRUE))

Pec.Pb.PeakNet.Use <- PeakNetFiles[which(substring(PeakNetFiles, 114, 131)=='Pec_Pb_PeakNetGRF_')]
Pec.Pb.PeakNet <- data.frame(read.csv(Pec.Pb.PeakNet.Use, header=TRUE))

Pel.Af.PeakNet.Use <- PeakNetFiles[which(substring(PeakNetFiles, 114, 131)=='Pel_Af_PeakNetGRF_')]
Pel.Af.PeakNet <- data.frame(read.csv(Pel.Af.PeakNet.Use, header=TRUE))

Pec.Pw.PeakNet.Use <- PeakNetFiles[which(substring(PeakNetFiles, 114, 131)=='Pec_Pw_PeakNetGRF_')]
Pec.Pw.PeakNet <- data.frame(read.csv(Pec.Pw.PeakNet.Use, header=TRUE))

Pel.Pw.PeakNet.Use <- PeakNetFiles[which(substring(PeakNetFiles, 114, 131)=='Pel_Pw_PeakNetGRF_')]
Pel.Pw.PeakNet <- data.frame(read.csv(Pel.Pw.PeakNet.Use, header=TRUE))


# Subsetting file to only analyze the variables of interest
Pec.Af.PeakNet.Subset <- data.frame(Pec.Af.PeakNet[,c(1,10,12:20)], "Af.Pec"); names(Pec.Af.PeakNet.Subset)[12] <- "group"
Pec.Pb.PeakNet.Subset <- data.frame(Pec.Pb.PeakNet[,c(1,10,12:20)], "Pb.Pec"); names(Pec.Pb.PeakNet.Subset)[12] <- "group"
Pel.Af.PeakNet.Subset <- data.frame(Pel.Af.PeakNet[,c(1,10,12:20)], "Af.Pel"); names(Pel.Af.PeakNet.Subset)[12] <- "group"
Pec.Pw.PeakNet.Subset <- data.frame(Pec.Pw.PeakNet[,c(1,10,12:20)], "Pw.Pec"); names(Pec.Pw.PeakNet.Subset)[12] <- "group"
Pel.Pw.PeakNet.Subset <- data.frame(Pel.Pw.PeakNet[,c(1,10,12:20)], "Pw.Pel"); names(Pel.Pw.PeakNet.Subset)[12] <- "group"

# Combining all the datasets together
Full.PeakNet <- rbind(Pec.Af.PeakNet.Subset, Pec.Pb.PeakNet.Subset, Pel.Af.PeakNet.Subset, Pec.Pw.PeakNet.Subset, Pel.Pw.PeakNet.Subset)

# Datasets for comparison
Pec.PeakNet.Subset <- rbind(Pec.Af.PeakNet.Subset, Pec.Pb.PeakNet.Subset, Pec.Pw.PeakNet.Subset)
FL.PeakNet.Subset <- rbind(Pec.Af.PeakNet.Subset, Pec.Pw.PeakNet.Subset)
Af.PeakNet.Subset <- rbind(Pec.Af.PeakNet.Subset, Pel.Af.PeakNet.Subset)
Pw.PeakNet.Subset <- rbind(Pec.Pw.PeakNet.Subset, Pel.Pw.PeakNet.Subset)
Pel.PeakNet.Subset <- rbind(Pel.Af.PeakNet.Subset, Pel.Pw.PeakNet.Subset)
PwPb.PeakNet.Subset <- rbind(Pec.Pw.PeakNet.Subset, Pec.Pb.PeakNet.Subset)


############## CANONICAL DISCRIMINANT ANALYSIS   ##########
#Centering and scaling data
#Full.PeakNet.Scale <- data.frame(scale(Full.PeakNet[,1:7], center=T, scale=T))

library(candisc)
## Comparison between groups
# Must begin with a linear model formulation
attach(Full.PeakNet)
Full.Y <- cbind(PercentStance, APAngle.Convert, MLAngle.Convert, InterpV.BW, InterpML.BW, InterpHz.BW, NetGRF.BW)
#Full.Y <- data.frame(Full.Y)
#names(Full.Y) <- names(Full.PeakNet)[1:7]

### Comparison of all appendages
# Running analysis 
Full.PeakNet$Group <- paste(Full.PeakNet$Species,".", Full.PeakNet$Appendage, sep="")

Full.LM <- manova(Full.Y~as.factor(Group), data=Full.PeakNet)
Full.CDA <- candisc(Full.LM, type="II")   # Type II error
#  Computing the MANOVA statistics
#  Currently set to use Wilks Lambda, but can do others
#  test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
summary(Full.LM, test=c("Wilks"))
detach(Full.PeakNet)

Full.CDA.Loadings <- Full.CDA$coeffs.std  # you want the standardized coefficients rather than raw

# Checking to see if Type II error produces same results as Type I error for MANOVA
# requires the car package for Manova()
Full.MANOVA.TypeII <- Manova(lm(Full.Y~Group, data=Full.PeakNet), type="II")

quartz()
plot(Full.CDA)

# Testing significance of variables against the axes; 1st column of DFA is the group, so DF1 is in column 2
DF1.Cor <- sapply(Full.PeakNet[,1:7], FUN=function(x) cor.test(x, Full.CDA$scores[,2], method="spearman",exact=F))
DF2.Cor <- sapply(Full.PeakNet[,1:7], FUN=function(x) cor.test(x, Full.CDA$scores[,3], method="spearman",exact=F))
DF3.Cor <- sapply(Full.PeakNet[,1:7], FUN=function(x) cor.test(x, Full.CDA$scores[,4], method="spearman",exact=F))
DF4.Cor <- sapply(Full.PeakNet[,1:7], FUN=function(x) cor.test(x, Full.CDA$scores[,5], method="spearman",exact=F))


# Calculating centroid of each group
pw.Pec.DFA.Centroid <- sapply(Full.CDA$scores[Full.PeakNet$Group=="pw.Pec",2:3], FUN=function(x) mean(x))
pw.Pel.DFA.Centroid <- sapply(Full.CDA$scores[Full.PeakNet$Group=="pw.Pel",2:3], FUN=function(x) mean(x))
af.Pec.DFA.Centroid <- sapply(Full.CDA$scores[Full.PeakNet$Group=="af.Pec",2:3], FUN=function(x) mean(x))
af.Pel.DFA.Centroid <- sapply(Full.CDA$scores[Full.PeakNet$Group=="af.Pel",2:3], FUN=function(x) mean(x))
pb.Pec.DFA.Centroid <- sapply(Full.CDA$scores[Full.PeakNet$Group=="pb.Pec",2:3], FUN=function(x) mean(x))

DFA.Centroids <- rbind(pw.Pec.DFA.Centroid,pw.Pel.DFA.Centroid, af.Pec.DFA.Centroid, af.Pel.DFA.Centroid, pb.Pec.DFA.Centroid)

quartz()
# Since the first column are the group names, canonical axis 1 = column 2
# Plotting axis 1 vs. axis 2
plot(Full.CDA$scores[,2], Full.CDA$scores[,3], type='n', xlab=paste('DF 1 (', round(Full.CDA$pct[1],2), "%)", sep=""), ylab=paste('DF 2 (', round(Full.CDA$pct[2],2), "%)", sep=""))
points(Full.CDA$scores[Full.PeakNet$Group=="pw.Pec",2], Full.CDA$scores[Full.PeakNet$Group=="pw.Pec",3], pch=21, bg='lightblue3')
points(Full.CDA$scores[Full.PeakNet$Group=="pw.Pel",2], Full.CDA$scores[Full.PeakNet$Group=="pw.Pel",3], pch=21, bg='steelblue3')
points(Full.CDA$scores[Full.PeakNet$Group=="af.Pec",2], Full.CDA$scores[Full.PeakNet$Group=="af.Pec",3], pch=21, bg='burlywood')
points(Full.CDA$scores[Full.PeakNet$Group=="af.Pel",2], Full.CDA$scores[Full.PeakNet$Group=="af.Pel",3], pch=21, bg='rosybrown4')
points(Full.CDA$scores[Full.PeakNet$Group=="pb.Pec",2], Full.CDA$scores[Full.PeakNet$Group=="pb.Pec",3], pch=25, bg='seagreen')
legend("topleft", title="Appendage Type", c("Aquatic FL", "Aquatic HL", "Terrestrial FL", "Terrestrial HL", "Terrestrial PF"),
       #pch=c(21,21,21,21,25), 
       fill=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))
# Adding centroids to the plot
points(pw.Pec.DFA.Centroid[1], pw.Pec.DFA.Centroid[2], pch=22, col="black", bg="lightblue3", cex=2)
points(pw.Pel.DFA.Centroid[1], pw.Pel.DFA.Centroid[2], pch=22, col="black", bg="steelblue3", cex=2)
points(af.Pec.DFA.Centroid[1], af.Pec.DFA.Centroid[2], pch=22, col="black", bg="burlywood", cex=2)
points(af.Pel.DFA.Centroid[1], af.Pel.DFA.Centroid[2], pch=22, col="black", bg="rosybrown4", cex=2)
points(pb.Pec.DFA.Centroid[1], pb.Pec.DFA.Centroid[2], pch=22, col="black", bg="seagreen", cex=2)

## Drawing convex hulls
Full.CDA.scores <- Full.CDA$scores
names(Full.CDA.scores)[1] <- "Group"

library(plyr) # required for ddply
find_hull <- function(Full.CDA.scores) Full.CDA.scores[chull(Full.CDA.scores[,2], Full.CDA.scores[,3]),]
hulls <- ddply(Full.CDA.scores, "Group", find_hull)

library(ggplot2)
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 5 Compare Data')
tiff(filename=paste("DFAPlot_NoOutlier_", SaveDate, ".tif", sep=""), width=600, height=600)
ggplot(data = Full.CDA.scores, aes(x = Can1, y = Can2, colour=Group, fill = Group)) +
  #scale_x_continuous(element_blank())+
  scale_x_continuous(paste('DF 1 (', round(Full.CDA$pct[1],2), "%)", sep=""))+
  scale_y_continuous(paste('DF 2 (', round(Full.CDA$pct[2],2), "%)", sep=""))+
  geom_point() + 
  geom_polygon(data = hulls, alpha = 0.5)+
  scale_color_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Terrestrial PF    ", "Semi-Aquatic FL  ", "Semi-Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("burlywood", "rosybrown4", "seagreen", "lightblue3", "steelblue3"))+
  scale_fill_manual(name="Appendage:", 
                     labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Terrestrial PF    ", "Semi-Aquatic FL  ", "Semi-Aquatic HL  "),
                     #values=c("ivory4", "ivory4", "ivory4"))+
                     values=c("burlywood", "rosybrown4", "seagreen", "lightblue3", "steelblue3"))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="bottom")
dev.off()

###  LINEAR DISCRIMINANT ANALYSIS
# Default CV is FALSE
# used to predict group membership
Full.LDA <- lda(Group~., data=Full.PeakNet[,c(1:7,13)], CV=FALSE)

# Running Cross-validation procedure
Full.CV.prior <- lda(Group~., data=Full.PeakNet[,c(1:7,13)], CV=FALSE)
# Gives you the number of trials classified in each group
Full.CV.prior$counts

Full.CV.posterior <- lda(Group~., data=Full.PeakNet[,c(1:7,13)], CV=TRUE)
# Use this code to get the crossvalidation results
round(Full.CV.posterior$posterior)
table(round(Full.CV.posterior$posterior))
# Gives you the number of trials that were classified in each group
table(Full.CV.posterior$class)
# Calculating the percentage of correct classifications
Full.PercentCorrect <- sum(diag(prop.table(table(Full.PeakNet$Group, Full.CV.posterior$class))))*100
Full.ct <- table(Full.PeakNet$Group, Full.CV.posterior$class)
diag(prop.table(Full.ct, 1))
# total percent correct
sum(diag(prop.table(Full.ct)))

# # Determining what trials were misclassified
# Full.allindices <- 1:nrow(Full.PeakNet)
# Full.misclassified <- Full.allindices[ predict(Full.LDA)$class != Full.PeakNet$Group]
# Full.right <- Full.PeakNet$Group[Full.misclassified]
# Full.wrong <- data.frame(predict(Full.LDA)$class[Full.misclassified])
# Full.MisclassID <- cbind(Full.right, Full.wrong)
# names(Full.MisclassID) <- c("Correct", "Misclassified")

# Generating table of what trials were misclassified
# The correct identities are the rows
Full.Misclass <- table(Full.PeakNet$Group, predict(Full.LDA)$class)


Species.LM <- manova(Species.Y~Species, data=Full.PeakNet)
Species.CDA <- candisc(Species.LM)   
#  Computing the MANOVA statistics
#  Currently set to use Wilks Lambda, but can do others
#  test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
summary(Species.LM, test=c("Wilks"))
detach(Full.PeakNet)


plot(Species.CDA)


# Since the first column are the group names, canonical axis 1 = column 2
# Plotting axis 1 vs. axis 2
plot(Species.CDA$scores[,2], Species.CDA$scores[,3], type='n', xlab='DF 1', ylab='DF 2')
points(Species.CDA$scores[Full.PeakNet$Species=="pw",2], Species.CDA$scores[Full.PeakNet$Species=="pw",3], pch=21, bg='blue')
points(Species.CDA$scores[Full.PeakNet$Species=="af",2], Species.CDA$scores[Full.PeakNet$Species=="af",3], pch=21, bg='tan')
points(Species.CDA$scores[Full.PeakNet$Species=="pb",2], Species.CDA$scores[Full.PeakNet$Species=="pb",3], pch=25, bg='green')
legend("bottomright", title="Species", c("Pw", "Af", "Pb"),
       #pch=c(21,21,21,21,25), 
       fill=c("blue", "tan", "green"))

#  LINEAR DISCRIMINANT ANALYSIS
# used to predict group membership
Species.LDA <- lda(Species~., data=Full.PeakNet[,c(1:7,10)], CV=FALSE)

# Running Cross-validation procedure
Species.CV.prior <- lda(Species~., data=Full.PeakNet[,c(1:7,10)], CV=FALSE)
# Gives you the number of trials classified in each group
Species.CV.prior$counts

Species.CV.posterior <- lda(Species~., data=Full.PeakNet[,c(1:7,10)], CV=TRUE)
# Use this code to get the crossvalidation results
round(Species.CV.posterior$posterior)
table(round(Species.CV.posterior$posterior))
# Gives you the number of trials that were classified in each group
table(Species.CV.posterior$class)
# Calculating the percentage of correct classifications
Species.PercentCorrect <- sum(diag(prop.table(table(Full.PeakNet$Species, Species.CV.posterior$class))))*100
Species.ct <- table(Full.PeakNet$Species, Species.CV.posterior$class)
diag(prop.table(Species.ct, 1))
# total percent correct
sum(diag(prop.table(Species.ct)))


##  Comparison of pectoral appendages
# Running analysis 
attach(Pec.PeakNet.Subset)
Pec.Y <- cbind(PercentStance, APAngle.Convert, MLAngle.Convert, InterpV.BW, InterpML.BW, InterpHz.BW, NetGRF.BW)
Pec.LM <- manova(Pec.Y~Species, data=Pec.PeakNet.Subset)
Pec.CDA <- candisc(Pec.LM)   
#  Computing the MANOVA statistics
#  Currently set to use Wilks Lambda, but can do others
#  test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
summary(Pec.LM, test=c("Wilks"))
detach(Pec.PeakNet.Subset)

Pec.PeakNet.Subset$Group <- paste(Pec.PeakNet.Subset$Species,".", Pec.PeakNet.Subset$Appendage, sep="")

plot(Pec.CDA)


# Since the first column are the group names, canonical axis 1 = column 2
# Plotting axis 1 vs. axis 2
plot(Pec.CDA$scores[,2], Pec.CDA$scores[,3], type='n', xlab='DF 1', ylab='DF 2')
points(Pec.CDA$scores[Pec.PeakNet.Subset$Group=="af.Pec",2], Pec.CDA$scores[Pec.PeakNet.Subset$Group=="af.Pec",3], pch=21, bg='tan')
points(Pec.CDA$scores[Pec.PeakNet.Subset$Group=="pw.Pec",2], Pec.CDA$scores[Pec.PeakNet.Subset$Group=="pw.Pec",3], pch=21, bg='blue')
points(Pec.CDA$scores[Pec.PeakNet.Subset$Group=="pb.Pec",2], Pec.CDA$scores[Pec.PeakNet.Subset$Group=="pb.Pec",3], pch=25, bg='green')
legend("bottomright", title="Species", c("Af", "Pw", "Pb"),
       #pch=c(21,21,21,21,25), 
       fill=c("tan", "blue", "green"))


#  LINEAR DISCRIMINANT ANALYSIS
# used to predict group membership
Pec.LDA <- lda(Group~., data=Pec.PeakNet.Subset[,c(1:7,12)], CV=FALSE)

# Running Cross-validation procedure
Pec.CV.prior <- lda(Group~., data=Pec.PeakNet.Subset[,c(1:7,12)], CV=FALSE)
# Gives you the number of trials classified in each group
Pec.CV.prior$counts

Pec.CV.posterior <- lda(Group~., data=Pec.PeakNet.Subset[,c(1:7,12)], CV=TRUE)
# Use this code to get the crossvalidation results
round(Pec.CV.posterior$posterior)
table(round(Pec.CV.posterior$posterior))
# Gives you the number of trials that were classified in each group
table(Pec.CV.posterior$class)
# Calculating the percentage of correct classifications
Pec.PercentCorrect <- sum(diag(prop.table(table(Pec.PeakNet.Subset$Group, Pec.CV.posterior$class))))*100
Pec.ct <- table(Pec.PeakNet.Subset$Group, Pec.CV.posterior$class)
diag(prop.table(Pec.ct, 1))
# total percent correct
sum(diag(prop.table(Pec.ct)))


### Plotting PCA
# prcomp is the preferred function from PCA (see documentation for more details)
# rotations are the loadings
# sdev are the square root of the eigenvalues
Y<-as.matrix(Full.PeakNet[,1:7])  
Y<-prcomp(Y, center=T, scale.=T)$x  #TO USE IF > 2 variables: for plotting PC axes at the end; scale=T makes it a correlation matrix
plot(Y[,c(1,2)],type="n",xlab="PC I", ylab="PC II",asp=1)
points(Y[Full.PeakNet$Species=="af" & Full.PeakNet$Appendage=="Pec",c(1,2)], pch=21, col="black", bg="darkgreen")
points(Y[Full.PeakNet$Species=="af" & Full.PeakNet$Appendage=="Pel",c(1,2)], pch=24, col="black", bg="green")
points(Y[Full.PeakNet$Species=="pw" & Full.PeakNet$Appendage=="Pec",c(1,2)], pch=21, col="black", bg="blue")
points(Y[Full.PeakNet$Species=="pw" & Full.PeakNet$Appendage=="Pel",c(1,2)], pch=24, col="black", bg="lightblue")
points(Y[Full.PeakNet$Species=="pb" & Full.PeakNet$Appendage=="Pec",c(1,2)], pch=21, col="black", bg="purple")

## I'm not sure if princomp is centering and scaling the data; might be better to stick with prcomp results from above
Full.PCA <- princomp(scale(Full.PeakNet[,1:7], center=T, scale=T), cor=TRUE, scores=TRUE)
Full.PCA.eigenvalue <- (Full.PCA$sdev)^2
Full.PCA.PerVar <- ((Full.PCA$sdev^2)/sum(Full.PCA$sdev^2))*100


# Since the first column are the group names, canonical axis 1 = column 2
## Plotting axis 1 vs. axis 2
plot(Full.PCA$scores[,2], Full.PCA$scores[,3], type='n', xlab=paste('PC 1 (', round(Full.PCA.PerVar[1],2), "%)", sep=""), ylab=paste('PC 2 (', round(Full.PCA.PerVar[2],2), "%)", sep=""))
points(Full.PCA$scores[Full.PeakNet$Group=="pw.Pec",2], Full.PCA$scores[Full.PeakNet$Group=="pw.Pec",3], pch=21, bg='lightblue3')
points(Full.PCA$scores[Full.PeakNet$Group=="pw.Pel",2], Full.PCA$scores[Full.PeakNet$Group=="pw.Pel",3], pch=21, bg='steelblue3')
points(Full.PCA$scores[Full.PeakNet$Group=="af.Pec",2], Full.PCA$scores[Full.PeakNet$Group=="af.Pec",3], pch=21, bg='burlywood')
points(Full.PCA$scores[Full.PeakNet$Group=="af.Pel",2], Full.PCA$scores[Full.PeakNet$Group=="af.Pel",3], pch=21, bg='rosybrown4')
points(Full.PCA$scores[Full.PeakNet$Group=="pb.Pec",2], Full.PCA$scores[Full.PeakNet$Group=="pb.Pec",3], pch=25, bg='seagreen')
legend("bottomright", title="Appendage Type", c("Pw FL", "Pw HL", "Af FL", "Af HL", "Pb PF"),
       #pch=c(21,21,21,21,25), 
       fill=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))


#############################  TESTING UNIVARIATE NORMALITY ####################

################################ RAW DATA ######################################
# Checking for normality for each indivdiual variable using Shapiro-Wilk test for normality
Pec.PeakNet.Subset.Norm <- sapply(Pec.PeakNet.Subset[,c(1:7)], FUN=function(x) shapiro.test(x))
Pec.PeakNet.Subset.NormPValues <- Pec.PeakNet.Subset.Norm[2,]
Pec.PeakNet.Subset.NumNormal <- sum(Pec.PeakNet.Subset.NormPValues > 0.05)
Pec.PeakNet.Subset.NormOutput <- data.frame(t(Pec.PeakNet.Subset.Norm[-4,])) # don't want to include the data.name info
Pec.PeakNet.Subset.NormOutput$Type <- 'Univariate'
Pec.PeakNet.Subset.NormOutput$Conclusion <- ifelse(Pec.PeakNet.Subset.NormOutput$p.value>0.05,'Normal','Not Normal')

Af.PeakNet.Subset.Norm <- sapply(Af.PeakNet.Subset[,c(1:7)], FUN=function(x) shapiro.test(x))
Af.PeakNet.Subset.NormPValues <- Af.PeakNet.Subset.Norm[2,]
Af.PeakNet.Subset.NumNormal <- sum(Af.PeakNet.Subset.NormPValues > 0.05)
Af.PeakNet.Subset.NormOutput <- data.frame(t(Af.PeakNet.Subset.Norm[-4,])) # don't want to include the data.name info
Af.PeakNet.Subset.NormOutput$Type <- 'Univariate'
Af.PeakNet.Subset.NormOutput$Conclusion <- ifelse(Af.PeakNet.Subset.NormOutput$p.value>0.05,'Normal','Not Normal')

Pw.PeakNet.Subset.Norm <- sapply(Pw.PeakNet.Subset[,c(1:7)], FUN=function(x) shapiro.test(x))
Pw.PeakNet.Subset.NormPValues <- Pw.PeakNet.Subset.Norm[2,]
Pw.PeakNet.Subset.NumNormal <- sum(Pw.PeakNet.Subset.NormPValues > 0.05)
Pw.PeakNet.Subset.NormOutput <- data.frame(t(Pw.PeakNet.Subset.Norm[-4,])) # don't want to include the data.name info
Pw.PeakNet.Subset.NormOutput$Type <- 'Univariate'
Pw.PeakNet.Subset.NormOutput$Conclusion <- ifelse(Pw.PeakNet.Subset.NormOutput$p.value>0.05,'Normal','Not Normal')

Pel.PeakNet.Subset.Norm <- sapply(Pel.PeakNet.Subset[,c(1:7)], FUN=function(x) shapiro.test(x))
Pel.PeakNet.Subset.NormPValues <- Pel.PeakNet.Subset.Norm[2,]
Pel.PeakNet.Subset.NumNormal <- sum(Pel.PeakNet.Subset.NormPValues > 0.05)
Pel.PeakNet.Subset.NormOutput <- data.frame(t(Pel.PeakNet.Subset.Norm[-4,])) # don't want to include the data.name info
Pel.PeakNet.Subset.NormOutput$Type <- 'Univariate'
Pel.PeakNet.Subset.NormOutput$Conclusion <- ifelse(Pel.PeakNet.Subset.NormOutput$p.value>0.05,'Normal','Not Normal')


###########################  TESTING VARIANCES ################################
# Variables with normal distribution will use Bartlett's test to test for similar variances betwen groups
# Non-parameteric alternative will be fligner.test

# Subset normal data
Pec.PeakNet.Subset.Norm <- rownames(Pec.PeakNet.Subset.NormOutput[which(Pec.PeakNet.Subset.NormOutput$Conclusion=="Normal"),])
Pec.PeakNet.Subset.NonNorm <- rownames(Pec.PeakNet.Subset.NormOutput[which(Pec.PeakNet.Subset.NormOutput$Conclusion=="Not Normal"),])

Pec.PeakNet.Subset.ParVar <- sapply(Pec.PeakNet.Subset[,1:7], FUN=function(x) bartlett.test(x ~ Pec.PeakNet.Subset$Species))
Pec.PeakNet.Subset.NonParVar <- sapply(Pec.PeakNet.Subset[,1:7], FUN=function(x) fligner.test(x ~ Pec.PeakNet.Subset$Species))

Af.PeakNet.Subset.ParVar <- sapply(Af.PeakNet.Subset[,1:7], FUN=function(x) bartlett.test(x ~ Af.PeakNet.Subset$Appendage))
Af.PeakNet.Subset.NonParVar <- sapply(Af.PeakNet.Subset[,1:7], FUN=function(x) fligner.test(x ~ Af.PeakNet.Subset$Appendage))

# For normal data
#PecPeakNet.ParVar <- sapply(Pec.PeakNet.Subset[,c(PecPeakNet.Norm)], FUN=function(x) bartlett.test(x ~ Pec.PeakNet.Subset$Species))

# For non-normal data
#PecPeakNet.NonParVar <- sapply(Pec.PeakNet.Subset[,c(PecPeakNet.NonNorm)], FUN=function(x) fligner.test(x ~ Pec.PeakNet.Subset$Species))

# Making dataframes of results
# Normal data
PecPeakNet.ParVar.Output <- data.frame(t(Pec.PeakNet.Subset.ParVar[-4,])) # don't want to include the data.name info; order is a little different from Fligner
PecPeakNet.ParVar.Output$Type <- 'Univariate'
PecPeakNet.ParVar.Output$Conclusion <- ifelse(PecPeakNet.ParVar.Output$p.value>0.05,'Equal Variance','Unequal Variance')

AfPeakNet.ParVar.Output <- data.frame(t(Af.PeakNet.Subset.ParVar[-4,])) # don't want to include the data.name info; order is a little different from Fligner
AfPeakNet.ParVar.Output$Type <- 'Univariate'
AfPeakNet.ParVar.Output$Conclusion <- ifelse(AfPeakNet.ParVar.Output$p.value>0.05,'Equal Variance','Unequal Variance')

# Non-normal data
PecPeakNet.NonParVar.Output <- data.frame(t(Pec.PeakNet.Subset.NonParVar[-5,])) # don't want to include the data.name info
PecPeakNet.NonParVar.Output$Type <- 'Univariate'
PecPeakNet.NonParVar.Output$Conclusion <- ifelse(PecPeakNet.NonParVar.Output$p.value>0.05,'Equal Variance','Unequal Variance')

AfPeakNet.NonParVar.Output <- data.frame(t(Af.PeakNet.Subset.NonParVar[-5,])) # don't want to include the data.name info
AfPeakNet.NonParVar.Output$Type <- 'Univariate'
AfPeakNet.NonParVar.Output$Conclusion <- ifelse(AfPeakNet.NonParVar.Output$p.value>0.05,'Equal Variance','Unequal Variance')


# Merging variance output together
PecPeakNet.Var <- rbind(PecPeakNet.ParVar.Output, PecPeakNet.NonParVar.Output)
AfPeakNet.Var <- rbind(AfPeakNet.ParVar.Output, AfPeakNet.NonParVar.Output)


# Saving the data on testing assumptions
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 5 Compare Data')
library(xlsx)
PecPeakNet.Filename <- paste("PecPeakNet_TestAss_", SaveDate, ".xlsx", sep="")
write.xlsx(Pec.PeakNet.Subset.NormOutput[order(rownames(Pec.PeakNet.Subset.NormOutput)),], file=PecPeakNet.Filename, sheetName="Uni Norm")
write.xlsx(PecPeakNet.Var[order(rownames(PecPeakNet.Var)),], file=PecPeakNet.Filename, sheetName="Variance", append=TRUE)

AfPeakNet.Filename <- paste("AfPeakNet_TestAss_", SaveDate, ".xlsx", sep="")
write.xlsx(Af.PeakNet.Subset.NormOutput[order(rownames(Af.PeakNet.Subset.NormOutput)),], file=AfPeakNet.Filename, sheetName="Uni Norm")
write.xlsx(AfPeakNet.Var[order(rownames(AfPeakNet.Var)),], file=AfPeakNet.Filename, sheetName="Variance", append=TRUE)


# Identifying which variables should be tested using parametric vs. non-parametric analyses
#if (substring(PecPeakNet.Var$method, 1,8)=='Bartlett' & PecPeakNet.Var$Conclusion == 'Equal Variance')
#  {PecPeakNet.Var$Decision <- 'Parametric'}


######################  COMPARING MEANS BETWEEN SPECIES ####################################

###  Linear mixed effects model 
# Recommended by Billy Bridges and Julia Sharp since I need to incorporate random effects, which aov() can't do
# Also, since we used the same individuals between the FL and HL comparison, we need to develop a different
# model from the FL vs. PF comparison 
# The lmer function in the lme4 package offers greater strength than aov() because it allows for the inclusion
# of random effects, and it does Type III error (aov() only does Type I)
# pvals.fnc from the languageR package is used to extract p-values for mixed models
# can tested for normality of residuals by typing: plot(fitted(PecNetGRFBW.LMER), residuals(PecNetGRFBW.LMER));
# should see a shot-gun splatter (i.e., no linear trends)
# use pMCMC values from pvals.fnc 'cause the Pr(>|t|) p-values are less conservative, and not as reliable
# why p-values aren't generated in lmer: https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html

## Talked to Saara, and it's okay that there is no p-value for the random effect since you're not trying to test
# whether it significantly affects your dependent variables.
# My including random effects into your model, you're trying to say "given the random effects of this factor,
# does this fixed effect still explain my dependent variable?"  (or something like that)
# Thus, the output for your fixed effect from the lmer() procedure with random effects specified is the 
# influence of the fixed effect with your random effect(s) already accounted for
# Thus, it doesn't really help you all that much to do a model selection procedure where you drop the random 
# effect term(s) to see whether the reduced model is the best fit for your data
# Also, typically (in ecology at least) you do either a model selection procedure OR hypothesis testing (not one
# and then the other).  That's because if you're choosing your best fit model, it's already assumed that the
# fixed effects you entered into the model significantly explain your data 
# In order to account for the multiple measures per individual, Saara actually recommended that I do model 
# selection, so that I can figure out whether I need to include the interaction between trial and invidiual as
# another random effect into my model.  THen you can use AIC to determine your best fit model (the best one is 
# the AIC value that is lowest, and there should be a difference of at least 2 units)
# With this, I can see whether I can go with the interaction random effect, or just stick with a more simple model
# However, the interaction between ind and group didn't really work because each trial has a unique ID 
# Talked to Mike Sears, and he said that we can stick with a fairly simplistic model with Species or Appendage Type
# as the fixed factor, and then Ind nested within Species or Appendage Type as the random effect
# Do that using the lme() function in nlme package, and then use Anova(x, Type=II) to use Type II error to 
# evaluate significance.
# Problem with a standard nested ANOVA is that it doesn't handle unbalanced designs very well, and it also 
# calculates your error differently

#####  THIS IS THE ACTUAL METHOD THAT WE USE FOR OUR ANALYSES
##  Evaluating significance
library(lme4)
library(sciplot) # for calculating standard errors
library(boot) # for bootstrapping
library(lsmeans) #for Tukey post hoc comparison
library(pbkrtest)
# use the p-values in pMCMC (posterior distribution from a Markov Chain Monte Carlo)
# The intercept is the baseline (usually the group name starting with a letter earlier in the alphabet)
# The way you interpret the output is the the pMCMC value will tell you whether that group (e.g., pw)
# is significantly different from the baseline or intercept


###  The MCMC sampling from before isn't acceptable anymore, so going to use alternative approaches
##  Bolker et al. (2008) described how likelihood ratio tests are not recommended for GLMMs
## Although I'm not using the results from the LRT, I'm running it just to see the results
## Tukey post hoc comparisons on the least squares means are then being use to look at pairwise comparisons
## Although REML is less biased, you need to have Maximum Likelihood to calculate your AIC and BIC values and to conduct the LRT

PercentStance.LMER <- lmer(PercentStance~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
PercentStance.LMER.Null <- lmer(PercentStance~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
PercentStance.LMER.Compare <- anova(PercentStance.LMER, PercentStance.LMER.Null)
PercentStance.Tukey <- lsmeans(PercentStance.LMER, pairwise~group, adjust="tukey")
PercentStance.Means <- aggregate(Full.PeakNet$PercentStance, by=list(Full.PeakNet$group), FUN=mean)[2]
PercentStance.SE <- aggregate(Full.PeakNet$PercentStance, by=list(Full.PeakNet$group), FUN=se)[2]
PercentStanceSummaryStats <- data.frame(unique(Full.PeakNet$group), PercentStance.Means, PercentStance.SE)
names(PercentStanceSummaryStats) <- c("group", "mean", "se")

APAngleConvert.LMER <- lmer(APAngle.Convert~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
APAngleConvert.LMER.Null <- lmer(APAngle.Convert~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
APAngleConvert.LMER.Compare <- anova(APAngleConvert.LMER, APAngleConvert.LMER.Null)
APAngleConvert.Tukey <- lsmeans(APAngleConvert.LMER, pairwise~group, adjust="tukey")
APAngleConvert.Means <- aggregate(Full.PeakNet$APAngle.Convert, by=list(Full.PeakNet$group), FUN=mean)[2]
APAngleConvert.SE <- aggregate(Full.PeakNet$APAngle.Convert, by=list(Full.PeakNet$group), FUN=se)[2]
APAngleConvertSummaryStats <- data.frame(unique(Full.PeakNet$group), APAngleConvert.Means, APAngleConvert.SE)
names(APAngleConvertSummaryStats) <- c("group", "mean", "se")

MLAngleConvert.LMER <- lmer(MLAngle.Convert~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
MLAngleConvert.LMER.Null <- lmer(MLAngle.Convert~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
MLAngleConvert.LMER.Compare <- anova(MLAngleConvert.LMER, MLAngleConvert.LMER.Null)
MLAngleConvert.Tukey <- lsmeans(MLAngleConvert.LMER, pairwise~group, adjust="tukey")
MLAngleConvert.Means <- aggregate(Full.PeakNet$MLAngle.Convert, by=list(Full.PeakNet$group), FUN=mean)[2]
MLAngleConvert.SE <- aggregate(Full.PeakNet$MLAngle.Convert, by=list(Full.PeakNet$group), FUN=se)[2]
MLAngleConvertSummaryStats <- data.frame(unique(Full.PeakNet$group), MLAngleConvert.Means, MLAngleConvert.SE)
names(MLAngleConvertSummaryStats) <- c("group", "mean", "se")

InterpVBW.LMER <- lmer(InterpV.BW~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
InterpVBW.LMER.Null <- lmer(InterpV.BW~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
InterpVBW.LMER.Compare <- anova(InterpVBW.LMER, InterpVBW.LMER.Null)
InterpVBW.Tukey <- lsmeans(InterpVBW.LMER, pairwise~group, adjust="tukey")
InterpVBW.Means <- aggregate(Full.PeakNet$InterpV.BW, by=list(Full.PeakNet$group), FUN=mean)[2]
InterpVBW.SE <- aggregate(Full.PeakNet$InterpV.BW, by=list(Full.PeakNet$group), FUN=se)[2]
InterpVBWSummaryStats <- data.frame(unique(Full.PeakNet$group), InterpVBW.Means, InterpVBW.SE)
names(InterpVBWSummaryStats) <- c("group", "mean", "se")

InterpMLBW.LMER <- lmer(InterpML.BW~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
InterpMLBW.LMER.Null <- lmer(InterpML.BW~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
InterpMLBW.LMER.Compare <- anova(InterpMLBW.LMER, InterpMLBW.LMER.Null)
InterpMLBW.Tukey <- lsmeans(InterpMLBW.LMER, pairwise~group, adjust="tukey")
InterpMLBW.Means <- aggregate(Full.PeakNet$InterpML.BW, by=list(Full.PeakNet$group), FUN=mean)[2]
InterpMLBW.SE <- aggregate(Full.PeakNet$InterpML.BW, by=list(Full.PeakNet$group), FUN=se)[2]
InterpMLBWSummaryStats <- data.frame(unique(Full.PeakNet$group), InterpMLBW.Means, InterpMLBW.SE)
names(InterpMLBWSummaryStats) <- c("group", "mean", "se")

InterpHzBW.LMER <- lmer(InterpHz.BW~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
InterpHzBW.LMER.Null <- lmer(InterpHz.BW~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
InterpHzBW.LMER.Compare <- anova(InterpHzBW.LMER, InterpHzBW.LMER.Null)
InterpHzBW.Tukey <- lsmeans(InterpHzBW.LMER, pairwise~group, adjust="tukey")
InterpHzBW.Means <- aggregate(Full.PeakNet$InterpHz.BW, by=list(Full.PeakNet$group), FUN=mean)[2]
InterpHzBW.SE <- aggregate(Full.PeakNet$InterpHz.BW, by=list(Full.PeakNet$group), FUN=se)[2]
InterpHzBWSummaryStats <- data.frame(unique(Full.PeakNet$group), InterpHzBW.Means, InterpHzBW.SE)
names(InterpHzBWSummaryStats) <- c("group", "mean", "se")

NetGRFBW.LMER <- lmer(NetGRF.BW~group+(1|Ind), Full.PeakNet, verbose=T, REML=F)
NetGRFBW.LMER.Null <- lmer(NetGRF.BW~1+(1|Ind), Full.PeakNet, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
NetGRFBW.LMER.Compare <- anova(NetGRFBW.LMER, NetGRFBW.LMER.Null)
NetGRFBW.Tukey <- lsmeans(NetGRFBW.LMER, pairwise~group, adjust="tukey")
NetGRFBW.Means <- aggregate(Full.PeakNet$NetGRF.BW, by=list(Full.PeakNet$group), FUN=mean)[2]
NetGRFBW.SE <- aggregate(Full.PeakNet$NetGRF.BW, by=list(Full.PeakNet$group), FUN=se)[2]
NetGRFBWSummaryStats <- data.frame(unique(Full.PeakNet$group), NetGRFBW.Means, NetGRFBW.SE)
names(NetGRFBWSummaryStats) <- c("group", "mean", "se")

# May be better to run DFA and then calculate the mahalanobis distance to do pairwise comparisons of the groups

# Testing homogeneity of variances/covariances
Full.PeakNet.Bartlett <-  bartlett.test(Full.PeakNet[,1:7], Species, data=Full.PeakNet)
Pec.PeakNet.Bartlett <-  bartlett.test(Pec.PeakNet.Subset[,1:7], Species, data=Pec.PeakNet.Subset)
Pel.PeakNet.Bartlett <-  bartlett.test(Pel.PeakNet.Subset[,1:7], Species, data=Pel.PeakNet.Subset)

# Testing multivariate normality (null = data are normal)
library(mvnormtest)
Full.PeakNet.MNorm <- mshapiro.test(as.matrix(t(Full.PeakNet[,1:7])))
Pec.PeakNet.MNorm <- mshapiro.test(as.matrix(t(Pec.PeakNet.Subset[,1:7])))
Pel.PeakNet.MNorm <- mshapiro.test(as.matrix(t(Pel.PeakNet.Subset[,1:7])))

# Graphical Assessment of Multivariate Normality
Full.x <- as.matrix(Full.PeakNet[,1:7]) # n x p numeric matrix
Full.center <- colMeans(Full.x) # centroid
Full.n <- nrow(Full.x); Full.p <- ncol(Full.x); Full.cov <- cov(Full.x);
Full.d <- mahalanobis(Full.x,Full.center,Full.cov) # distances
qqplot(qchisq(ppoints(Full.n),df=Full.p),Full.d,
       main="QQ Plot Assessing Multivariate Normality - Full", ylab="Mahalanobis D2")
abline(a=0,b=1)  # abline adds straight lines to a figure (a=slope, b=intercept)

Pec.x <- as.matrix(Pec.PeakNet.Subset[,1:7]) # n x p numeric matrix
Pec.center <- colMeans(Pec.x) # centroid
Pec.n <- nrow(Pec.x); Pec.p <- ncol(Pec.x); Pec.cov <- cov(Pec.x);
Pec.d <- mahalanobis(Pec.x,Pec.center,Pec.cov) # distances
qqplot(qchisq(ppoints(Pec.n),df=Pec.p),Pec.d,
       main="QQ Plot Assessing Multivariate Normality - Pec", ylab="Mahalanobis D2")
abline(a=0,b=1)  # abline adds straight lines to a figure (a=slope, b=intercept)

Pel.x <- as.matrix(Pel.PeakNet.Subset[,1:7]) # n x p numeric matrix
Pel.center <- colMeans(Pel.x) # centroid
Pel.n <- nrow(Pel.x); Pel.p <- ncol(Pel.x); Pel.cov <- cov(Pel.x);
Pel.d <- mahalanobis(Pel.x,Pel.center,Pel.cov) # distances
qqplot(qchisq(ppoints(Pel.n),df=Pel.p),Pel.d,
       main="QQ Plot Assessing Multivariate Normality - Pel", ylab="Mahalanobis D2")
abline(a=0,b=1)  # abline adds straight lines to a figure (a=slope, b=intercept)
