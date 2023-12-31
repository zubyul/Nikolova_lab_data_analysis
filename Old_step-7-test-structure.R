library(ggplot2)

setwd('/Users/Administrator/Downloads/YZ-20191122T154314Z-001/YZ')

### read in data

genetic1 <- read.csv("scores_frontal_cortex_predicted_expression.csv",header=T,stringsAsFactors=F)
clinical <- read.csv("HCP_behavioral_data.csv",header=T,stringsAsFactors=F)
imaging <- read.csv("unrestricted_amymiles_11_18_2019_13_10_13.csv",stringsAsFactors=FALSE)
HCP <- read.delim("HCP_sample_ref.txt", stringsAsFactors = FALSE)
### subset by variables of interest in one big dataframe

genetic <- merge(HCP[,c(4,5)],genetic1,by.x='SAMPLE_ID', by.y='IID')
colnames(genetic)[2] <- "SID"
data <- merge(genetic[,c(2,3)],clinical[,c(1,2,10,11,204)],by.x='SID',by.y='Subject')
data <- merge(data,imaging,by.x='SID',by.y='Subject')

regressors <- data[,1:7]
regressors$Gender <- factor(regressors$Gender) 

sample <- subset(regressors,Race=='White')
sample <- subset(regressors,Ethnicity=='Not Hispanic/Latino')
sample <- sample[,-c(4,5)]

volume <- data[,c(1,30:33,37,38,40,48:54)]
thickness <- data[,c(1,70:137)]
surfarea <- data[,c(1,138:205)]


### test main effects on volume

test.volume <- merge(sample,volume,by='SID')

tval.volume = array()
pval.volume = array()

for(i in names(test.volume[,6:19])){
  var <- paste(i)
  lm <- summary(glm(test.volume[[i]]~test.volume$male.weighted+test.volume$Age_in_Yrs+test.volume$Gender+test.volume$FS_IntraCranial_Vol))
  tval <- lm$coefficients[2,3]
  pval <- lm$coefficients[2,4]
  tval.volume[[i]] <- assign(var,round(tval,3))
  pval.volume[[i]] <- assign(var,round(pval,3))}

tval.volume <- tval.volume[2:15]
pval.volume <- pval.volume[2:15]
fdr.volume <- round(p.adjust(pval.volume,method='fdr'),3)
out.volume <- as.data.frame(cbind(tval.volume,pval.volume,fdr.volume))
out.volume$ROI <- names(test.volume[,6:19])
out.volume <- out.volume[,c(4,1,2,3)]
rownames(out.volume) <- 1:14
out.volume <- out.volume[order(out.volume$fdr.volume),]
sig.volume <- subset(out.volume,fdr.volume<0.05)

### test main effects on thickness

test.thickness <- merge(sample[,1:4],thickness,by='SID')

tval.thickness = array()
pval.thickness = array()

for(i in names(test.thickness[5:72])){
  var <- paste(i)
  lm <- summary(glm(test.thickness[[i]]~test.thickness$male.weighted+test.thickness$Age_in_Yrs+test.thickness$Gender))
  tval <- lm$coefficients[2,3]
  pval <- lm$coefficients[2,4]
  tval.thickness[[i]] <- assign(var,round(tval,3))
  pval.thickness[[i]] <- assign(var,round(pval,3))}

tval.thickness <- tval.thickness[2:69]
pval.thickness <- pval.thickness[2:69]
fdr.thickness <- round(p.adjust(pval.thickness,method='fdr'),3)
out.thickness <- as.data.frame(cbind(tval.thickness,pval.thickness,fdr.thickness))
out.thickness$ROI <- names(test.thickness[,5:72])
out.thickness <- out.thickness[,c(4,1,2,3)]
rownames(out.thickness) <- 1:68
out.thickness <- out.thickness[order(out.thickness$fdr.thickness),]
sig.thickness <- subset(out.thickness,fdr.thickness<0.05)

### test main effects on surfarea

test.surfarea <- merge(sample,surfarea,by='SID')

tval.surfarea = array()
pval.surfarea = array()

for(i in names(test.surfarea[,6:73])){
  var <- paste(i)
  lm <- summary(glm(test.surfarea[[i]]~test.surfarea$male.weighted
                    +test.surfarea$Age_in_Yrs+test.surfarea$Gender
                    +test.surfarea$FS_IntraCranial_Vol))
  tval <- lm$coefficients[2,3]
  pval <- lm$coefficients[2,4]
  tval.surfarea[[i]] <- assign(var,round(tval,3))
  pval.surfarea[[i]] <- assign(var,round(pval,3))}

tval.surfarea <- tval.surfarea[2:69]
pval.surfarea <- pval.surfarea[2:69]
fdr.surfarea <- round(p.adjust(pval.surfarea,method='fdr'),3)
out.surfarea <- as.data.frame(cbind(tval.surfarea,pval.surfarea,fdr.surfarea))
out.surfarea$ROI <- names(test.surfarea[,6:73])
out.surfarea <- out.surfarea[,c(4,1,2,3)]
rownames(out.surfarea) <- 1:68
out.surfarea <- out.surfarea[order(out.surfarea$fdr.surfarea),]
sig.surfarea <- subset(out.surfarea,fdr.surfarea<0.05)

write.csv(out.surfarea,'out_SA.csv')
write.csv(out.thickness,'out_thic.csv')
write.csv(out.volume,'out_vol.csv')

