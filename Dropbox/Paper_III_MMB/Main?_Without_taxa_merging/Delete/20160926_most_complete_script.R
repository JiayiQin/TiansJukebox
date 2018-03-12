###### Import data sheets ------------------------------------------------------

# MMB is for datasheet process
# OTU_ is for data analyse

#set the work routine
setwd("~/Dropbox/Paper_III_MMB/Without_taxa_merging/Delete")

###### Combine two OTU table according to the amplicon ---------------------------

MMB_combine <- read.csv("MMB.csv",header = T, sep =",",check.names=FALSE)

MMB_combine <- MMB_combine[order(MMB_combine$Taxa_result_CGG),]

write.csv(MMB_combine, "MMB_sort.csv")

TREAT<-read.csv("Data_structure.csv",header=TRUE,check.names=FALSE)

###### Organise taxa annotation  -- ROW ----------------------------------------
# MMB_combine is the table for further analysis
# 52 variables for RISO
# 156 variables for CGG

# Dim function to convert data from factor to character, *use paste(X)
FacToCha <- function(x) {if(is.factor(x)) 
{as.character(paste(x))} else {x}} 

MMB_combine$Taxa_result_CGG <- FacToCha(MMB_combine$Taxa_result_CGG)
MMB_combine$Taxa_result_riso <- FacToCha(MMB_combine$Taxa_result_riso)

# Manually check if the same amplicon got same taxonomy annotation in Excel
# Based on Taxa_result_CGG, copy taxa result 
TAXA <- MMB_combine$Taxa_result_CGG
for ( i in 1:672) {
  if (MMB_combine$Taxa_result_CGG[i] != "NA") 
  {TAXA[i] = MMB_combine$Taxa_result_CGG[i]}
  else
  {TAXA[i] = MMB_combine$Taxa_result_riso[i]}
}

MMB_taxa = cbind(TAXA,MMB_combine)

# IntToNum <- function(x) {if(is.integer(x)) 
# {as.numeric(paste(x))} else {x}} 
# MMB_taxa<- IntToNum(MMB_taxa)

# change the NA to value 0 for further analysis
MMB_taxa[is.na(MMB_taxa)]<- 0
# This command will lead to several warning message, 
# because there are several NAs as factor. No problem going further.
# This command will also change integer into numeric automatically.

# Minimize the table, remains the TAXA and OTU reads table
MMB_taxa_clean<- MMB_taxa[,-2:-7]
MMB_taxa_clean<- MMB_taxa_clean[,-54:-57] # extra column as character at the end of the table, need to delete
MMB_taxa_clean_onlynumbers<-MMB_taxa_clean[,2:209] # without taxa information
totalreads<-MMB_taxa_clean_onlynumbers[1,]
totalreads<-colSums(MMB_taxa_clean_onlynumbers)


###### Merge OTUs by taxa name ######################
# 2nd run merge -- COLUMN    #

# according to the original count
MMB_taxa_clean
MMB_merge <- MMB_taxa_clean[1,]

for (i in 2:672) #column
{
  if (MMB_taxa_clean$TAXA[i] %in% MMB_merge[,1]) # if the taxa already exists in rMMB_merge
  {
    position<- match(MMB_taxa_clean$TAXA[i], MMB_merge[,1]) # Find row number
    MMB_merge[position,2:209] <- MMB_merge[position,2:209] + MMB_taxa_clean[i,2:209] # Merge    
  }
  else
  {
    MMB_merge<- rbind(MMB_merge, MMB_taxa_clean[i,]) # not exist, add new line
  }  
}


## derep TruSeq replicate according to counts ####
MMB_merge
MMB_merge_derep <- MMB_merge[,1:105]
count2=0
count1=0
for (n in 1 : 52)
{ 
  r = 53 + 3 * (n-1) + 1 # locate the starting of calculation, colum
  for (m in 1:133)
  { # scan entir column
    T <- c(MMB_merge[m,r]>0, MMB_merge[m,r+1] >0, MMB_merge[m,r+2] >0) # logical calculation of the detection rate
    if (sum(T) > 1) # At least 2 out of 3 PCR got the OTU
    {
      count1 <- count1 + 1
      MMB_merge_derep[m,n+53] <- (MMB_merge[m,r]+MMB_merge[m,r+1]+MMB_merge[m,r+2])/3
    }
    else
    {
      MMB_merge_derep[m,n+53] <- 0
      count2 <- count2 +1
    }
  }
}

MMB_merge_derep_no0 <- MMB_merge_derep
MMB_merge_derep_no0 <- MMB_merge_derep[rowSums(MMB_merge_derep[,2:105]) > 0,]

write.csv(MMB_merge_derep_no0,"MMB_merge_derep_no0.csv")
#export the table and organise it in excel


################################ dataset ready for further organis ####################################

# The output should be MMB_merge_count to compatible with following analysis
# To do list for MMB_merge_derep_no0.csv
# 1. delete OTUs with the total reads lower than 5
# 2. give the the similarity between 90% - 97% a higher taxa annotation
# 3. delete column number 
# 4. delete rowSums
# 5. Organise the order for PM extraction

MMB_merge_count_pre<-read.csv("MMB_merge_derep_no0.csv")　
MMB_merge_count_pre<- MMB_merge_count_pre[,2:106]
sample<-replicate(104,"PB")

sample[1:24]<-replicate(24,"PB_NXT")
sample[25:44]<-replicate(20,"PM_NXT")
sample[45:52]<-replicate(8,"SP_NXT")
sample[53:76]<-replicate(24,"PB_TSQ")
sample[77:96]<-replicate(20,"PM_TSQ")
sample[97:104]<-replicate(8,"SP_TSQ")
a<-sample
sample<- replicate(105,"Taxa")
sample[2:105]<-a

names(MMB_merge_count_pre) = sample

MMB_merge_count <- MMB_merge_count_pre[1,]

for (i in 2:110) #column
{
  if (MMB_merge_count_pre$Taxa[i] %in% MMB_merge_count[,1]) # if the taxa already exists in rMMB_merge
  {
    position<- match(MMB_merge_count_pre$Taxa[i], MMB_merge_count [,1]) # Find row number
    MMB_merge_count [position,2:105] <- MMB_merge_count [position,2:105] + MMB_merge_count_pre[i,2:105] # Merge    
  }
  else
  {
    MMB_merge_count <- rbind(MMB_merge_count, MMB_merge_count_pre[i,]) # not exist, add new line
  }  
}

MMB_merge_count <- MMB_merge_count[rowSums(MMB_merge_count[,2:105])>5,]
#### calculate relative abundace for each OTUs ####
library(vegan)
MMB_merge_rltv<-decostand(MMB_merge_count[,2:105],"total",2)
# decostand function for standarization of ecology community data
# Method = total to calculate the relative abundance
# check<-colSums(MMB_taxa_rltv) # double check the calculation
# write.csv(cbind(TAXA,MMB_taxa_rltv),"MMB_taxa_relative.csv")

MMB_merge_rltv<-cbind(MMB_merge_count[,1],MMB_merge_rltv)


# Prepare dataset for further analysis ----------------------------------

#### readcounts, extract_method level ######
MMB_merge_count <- MMB_merge_count[rowSums(MMB_merge_count[,2:105])>5,]

PB_NXT_count<-MMB_merge_count[,2:25]
PM_NXT_count<-MMB_merge_count[,30:45] # exclude 10 g extraction
SP_NXT_count<-MMB_merge_count[,46:53]
PB_TSQ_count<-MMB_merge_count[,54:77]
PM_TSQ_count<-MMB_merge_count[,82:97] # exclude 10 g extraction
SP_TSQ_count<-MMB_merge_count[,98:105]

MMB_merge_count_extract<-rbind(rowSums(PB_NXT_count),
                                  rowSums(PM_NXT_count),
                                  rowSums(SP_NXT_count),
                                  rowSums(PB_TSQ_count),
                                  rowSums(PM_TSQ_count),
                                  rowSums(SP_TSQ_count))

MMB_merge_count_extract_rarefy<- MMB_merge_count_extract

# rarefy data set 
for ( i in 1: 6){
  for (j in 1:87)
  {
    if (MMB_merge_count_extract[i,j] > 0) {MMB_merge_count_extract_rarefy[i,j]=1} else {MMB_merge_count_extract_rarefy[i,j]=0}
  }
}
extract<-c("PB_NXT","PM_NXT","SP_NXT","PB_TSQ","PM_TSQ","SP_TSQ")
row.names(MMB_merge_count_extract_rarefy) <- extract
#### readcounts, soil core level ######
##### NXT
## PB_NXT_count
PB_NXT_count_derep =matrix(nrow=87,ncol=8)  #creat an empty matrix for PB
for ( i in 1:8)
{
  m<-(i-1)*3+1
  PB_NXT_count_derep[,i]<- rowSums(PB_NXT_count[,m:(m+2)])
}

## PM_NXT_count
PM_NXT_count_derep =matrix(nrow=87,ncol=8)  #creat an empty matrix for PM
for ( i in 1:8)
{
  m<-(i-1)*2+1
  PM_NXT_count_derep[,i]<- rowSums(PM_NXT_count[,m:(m+1)])
}

#### TSQ
##PB_TSQ_count
PB_TSQ_count_derep =matrix(nrow=87,ncol=8)  #creat an empty matrix for PB
for ( i in 1:8)
{
  m<-(i-1)*3+1
  PB_TSQ_count_derep[,i]<- rowSums(PB_TSQ_count[,m:(m+2)])
}

## PM_TSQ_count
PM_TSQ_count_derep =matrix(nrow=87,ncol=8)  #creat an empty matrix for PM
for ( i in 1:8)
{
  m<-(i-1)*2+1
  PM_TSQ_count_derep[,i]<- rowSums(PM_TSQ_count[,m:(m+1)])
}

PB_NXT_count_derep
PM_NXT_count_derep
SP_NXT_count<-SP_NXT_count[1:87,]
PB_TSQ_count_derep
PM_TSQ_count_derep
SP_TSQ_count<-SP_TSQ_count[1:87,]

MMB_merge_count_soilcore <- cbind(MMB_merge_count$Taxa,PB_NXT_count_derep,PM_NXT_count_derep,SP_NXT_count,PB_TSQ_count_derep,PM_TSQ_count_derep,SP_TSQ_count)

write.csv(MMB_merge_count_soilcore,"check_the_list.csv")
sample_soilcore<-replicate(48,"PB")
sample_soilcore<-c(replicate(8,"PB_NXT"),replicate(8,"PM_NXT"),replicate(8,"SP_NXT"),replicate(8,"PB_TSQ"),replicate(8,"PM_TSQ"),replicate(8,"SP_TSQ"))

################################ dataset ready to use ####################################

# network analysis --------------------------------------------------------
# prepare dataset 
library(igraph)

MMB_merge_count_extract_rarefy
rowSums(MMB_merge_count_extract_rarefy)

Network_rarefy<-as.matrix(MMB_merge_count_extract_rarefy)
network2graph <- graph_from_incidence_matrix(Network_rarefy)
network2graph.bp <- bipartite.projection(network2graph)
plot(network2graph.bp$proj1)
plot(network2graph.bp$proj2)

V(network2graph)$color <- c("steel blue", "orange")[V(network2graph)$type+1]
V(network2graph)$shape <- c("square", "circle")[V(network2graph)$type+1]
#V(network2graph)$label <- ""
#V(network2graph)$label[V(network2graph)$type==F] <- nodes2$media[V(network2graph)$type==F]
V(network2graph)$label.cex=.4
V(network2graph)$label.font=2

plot(network2graph, 
     vertex.size=7, 
     layout=layout_as_bipartite) # layout spicify the type of network


#reference from the package example
plot(net2, vertex.shape="none", vertex.label=NA,
     vertex.label.color=V(net2)$color, vertex.label.font=2.5,
     vertex.label.cex=.6, edge.color="gray70", edge.width=2, layout=layout_as_bipartite)



#### Venn diagram ####

library(gplots)
#prepare dataset
nr <- nrow(MMB_merge_count)
PBL_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PB_NXT_count[i,]) > 0) {PBL_NXT[i]= MMB_merge_count[i,1]}}
PBL_NXT <- PBL_NXT[PBL_NXT!="NA"]

PML_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PM_NXT_count[i,]) > 0) {PML_NXT[i]= MMB_merge_count[i,1]}}
PML_NXT <- PML_NXT[PML_NXT!="NA"]

SPL_NXT = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(SP_NXT_count[i,]) > 0) {SPL_NXT[i]= MMB_merge_count[i,1]}}
SPL_NXT <- SPL_NXT[SPL_NXT!="NA"]

PBL_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PB_TSQ_count[i,]) > 0) {PBL_TSQ[i]= MMB_merge_count[i,1]}}
PBL_TSQ <- PBL_TSQ[PBL_TSQ!="NA"]

PML_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(PM_TSQ_count[i,]) > 0) {PML_TSQ[i]= MMB_merge_count[i,1]}}
PML_TSQ <- PML_TSQ[PML_TSQ!="NA"]

SPL_TSQ = replicate(nr, "NA")
for (i in 1:nr)
{if (rowSums(SP_TSQ_count[i,]) > 0) {SPL_TSQ[i]= MMB_merge_count[i,1]}}
SPL_TSQ <- SPL_TSQ[SPL_TSQ!="NA"]

##### draw graph #####

input_1 <- list (SPL_NXT,SPL_TSQ)
venn(input_1)

input_2 <- list (PBL_TSQ,PML_TSQ,SPL_TSQ)
venn(input_2)

input_3 <-list (PBL_NXT,PML_NXT,SPL_NXT)
venn(input_3)

input_4<-list(PBL_NXT,PML_NXT,SPL_NXT,PBL_TSQ,PML_TSQ,SPL_TSQ)

library(VennDiagram)
# the input is based on the calculation from Venn function
# According to the paper described the package, it is impossible to get
# scaled venn diagram from this function for 3 sets data.
# Alternatively, I employed a website based function, BioVenn from Tim 
venn.plot <- venn.diagram(
  x = input_2,
  filename = "Venn_3set_complex.tiff",
  scaled = TRUE,
  euler.d=TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  category = c("Phophate buffer", "PowerMax","DNA soup"),
  fill = c("blue", "red", "orange"),
  lty = "blank",
  cex = 1,
  cat.cex = 1,
  cat.col = c("blue", "red", "orange"),
  cat.dist =0.05,
  rotation.degree = 45,
  main = "Venn Diagram"
  #sub = "Featuring: rotation and external lines",
  #main.cex = 2,
  #sub.cex = 1
)

grid.newpage()
grid.draw(venn.plot)

#### individual based rarefaction curve #####################

# The rarefaction curve produced from this part is based the individuals (abundance data)
library(vegan)

# dataset
# subset for each method
# use counts for each method
tPB_NXT_count_derep<-t(PB_NXT_count_derep)
tPM_NXT_count_derep<-t(PM_NXT_count_derep) 
tSP_NXT_count_derep<-t(SP_NXT_count)
tPB_TSQ_count_derep<-t(PB_TSQ_count_derep) 
tPM_TSQ_count_derep<-t(PM_TSQ_count_derep) 
tSP_TSQ_count_derep<-t(SP_TSQ_count)

## PB_NXT accumulation 
S_PB_NXT <- specnumber(tPB_NXT_count_derep) 
raremax_PB_NXT<-min(rowSums(tPB_NXT_count_derep))
Srare_PB_NXT <- rarefy(tPB_NXT_count_derep, raremax_PB_NXT)

plot(S_PB_NXT, Srare_PB_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## PM_NXT accumulation 
S_PM_NXT <- specnumber(tPM_NXT_count_derep) 
raremax_PM_NXT<-min(rowSums(tPM_NXT_count_derep))
Srare_PM_NXT <- rarefy(tPM_NXT_count_derep, raremax_PM_NXT)

plot(S_PM_NXT, Srare_PM_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## SP_NXT accumulation
S_SP_NXT <- specnumber(tSP_NXT_count) 
raremax_SP_NXT<-min(rowSums(tSP_NXT_count))
Srare_SP_NXT <- rarefy(tSP_NXT_count, raremax_SP_NXT)

plot(S_SP_NXT, Srare_SP_NXT, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

abline(0, 1)

## PB_TSQ accumulation 
S_PB_TSQ <- specnumber(round(tPB_TSQ_count_derep)) 
raremax_PB_TSQ<-min(rowSums(round(tPB_TSQ_count_derep)))
Srare_PB_TSQ <- rarefy(round(tPB_TSQ_count_derep), raremax_PB_TSQ)

plot(S_PB_TSQ, Srare_PB_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="blue")

abline(0, 1,col="blue")

## PM_TSQ accumulation 
S_PM_TSQ <- specnumber(round(tPM_TSQ_count_derep)) 
raremax_PM_TSQ<-min(rowSums(round(tPM_TSQ_count_derep)))
Srare_PM_TSQ <- rarefy(round(tPM_TSQ_count_derep), raremax_PM_TSQ)

plot(S_PM_TSQ, Srare_PM_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="red")

abline(0, 1,col="red")

## SP_TSQ accumulation
S_SP_TSQ <- specnumber(round(tSP_TSQ_count)) 
raremax_SP_TSQ<-min(rowSums(round(tSP_TSQ_count)))
Srare_SP_TSQ <- rarefy(round(tSP_TSQ_count), raremax_SP_TSQ)

plot(S_SP_TSQ, Srare_SP_TSQ, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",col="orange")

abline(0, 1,col="orange")

# export individual based rarefaction curve
par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)

par(mfrow=c(2,3))
rarecurve(tPB_NXT_count_derep, 
          sample = raremax_PB_NXT, 
          col = "blue", 
          xlim=c(0,40000),
          ylim=c(0,60),
          xaxt="n",
          label = FALSE,
          cex = 0.6)
mtext("Nextera", font=2,side=2,line=2.5,cex =1)
mtext("Phosphate Buffer",font=2,side=3,line=1.5,cex =1)
rarecurve(tPM_NXT_count_derep, 
          sample = raremax_PM_NXT, 
          col = "red", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          xaxt="n",
          yaxt="n",
          cex = 0.6)
mtext("Powermax",font=2,side=3,line=1.5,cex =1)
rarecurve(tSP_NXT_count, 
          sample = raremax_SP_NXT, 
          col = "orange", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          xaxt="n",
          yaxt="n",
          cex = 0.6)
mtext("DNA soup",font=2,side=3,line=1.5,cex =1)

rarecurve(round(tPB_TSQ_count_derep), 
          sample = raremax_PB_TSQ, 
          col = "blue", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          cex = 0.6)
mtext("TruSeq", font=2,side=2,line=2.5,cex =1)
rarecurve(round(tPM_TSQ_count_derep), 
          sample = raremax_PM_TSQ, 
          col = "red", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          yaxt="n",
          cex = 0.6)

rarecurve(round(tSP_TSQ_count), 
          sample = raremax_SP_TSQ, 
          col = "orange", 
          xlim=c(0,40000),
          ylim=c(0,60),
          label = FALSE,
          yaxt="n",
          cex = 0.6)

#### Sample based rarefaction curve #############
# Use vegan based package:rareNMtests to maximize the possibility of drawing rarefaction curve

library("rareNMtests")

# Ecological null model test using sample-based rarefaction curves

# for species richness (q = 0)

## sample based rarefaction curve for entire dataset #####

## sample based rarefaction with readcount in soil core level #########
#dataset
MMB_merge_count_soilcore
sample_soilcore

sample_rare_soilcore<-cbind(sample_soilcore,t(MMB_merge_count_soilcore_rafy))


sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[,-1], by=sample_rare_soilcore[,1], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix PB compare NXT and TSQ p=0.11
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[1:8,-1],sample_rare_soilcore[25:32,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix PM compare NXT and TSQ p=0.075
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[9:16,-1],sample_rare_soilcore[33:40,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# fix SP compare NXT and TSQ p=0.255
sbecoq_sample_rare_soilcore <- EcoTest.sample(rbind(sample_rare_soilcore[17:24,-1],sample_rare_soilcore[41:48,-1]), by=c(Library[1:8],Library[25:32]), MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

#Compare within NXT p=0.25
sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[1:24,-1], by=Extract[1:24], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

#Compare within TSQ p=0.07
sbecoq_sample_rare_soilcore <- EcoTest.sample(sample_rare_soilcore[25:48,-1], by=Extract[25:48], MARGIN=1)
plot(sbecoq_sample_rare_soilcore)

# library no difference regarding to the rarefaction curve
# extraction no difference regarding to the rarefaction cure as well

###### accumulation curve ######################
library(vegan)
# Accumulation curve

# Function specaccum finds species accumulation curves or the number of species for a certain number of sampled sites or individuals.

#dataset
PB_NXT_count_derep
PM_NXT_count_derep
SP_NXT_count
PB_TSQ_count_derep
PM_TSQ_count_derep
SP_TSQ_count

#alpha("blue", 0.1)
#adjustcolor( "blue", alpha.f = 0.2)
sp1 <- specaccum(t(PB_NXT_count_derep))
sp2 <- specaccum(t(PM_NXT_count_derep))
sp3 <- specaccum(t(SP_NXT_count))
sp4 <- specaccum(t(PB_TSQ_count_derep))
sp5 <- specaccum(t(PM_TSQ_count_derep))
sp6 <- specaccum(t(SP_TSQ_count))


#### draw accumulation curve -------------------------------------------------
par(mfrow=c(1,1))
plot(sp1,  
     xlim=c(1,8),
     ylim=c(15,80),
     main="Accumulation curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp2, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     col="red", 
     lwd=2, 
     ci.lty=0, 
     ci.type="poly", 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp3, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)
#boxplot(sp6, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
par(mfrow=c(1,1))
plot(sp4,  
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     #axes=F,
     ci.type="poly", 
     col="blue", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "lightblue", alpha.f = 0.5) #adjust the transparency 
)
#boxplot(sp2, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp5, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="",#main title of the graph
     xlab = "",#Lable for x axis
     ylab = "",# Lable for y axis
     axes=F,
     ci.type="poly", 
     col="red", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "pink", alpha.f = 0.5)
)
#boxplot(sp4, col="yellow", add=TRUE, pch="+")

par(new=TRUE)
plot(sp6, 
     xlim=c(1,8),
     ylim=c(15,80),
     main="Accumulation curve",#main title of the graph
     xlab = "Number of samples",#Lable for x axis
     ylab = "Number of species",# Lable for y axis
     ci.type="poly", 
     col="orange", 
     lwd=2, 
     lty=2,
     ci.lty=0, 
     ci.col=adjustcolor( "yellow", alpha.f = 0.5)
)


legend(4.5,25, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black", 
       title="Nextera",
       lty = c(1, 1, 1),
       bty="n",lwd=2,cex=1)
legend(4.5,25, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       #title="TruSeq",
       text.col = "black", 
       lty = c(2, 2, 2),
       bty="n",lwd=2,cex=1)

#boxplot(sp6, col="yellow", add=TRUE, pch="+")

# done


############################## relative abundance dataset ##################
# alpha diversity ---------------------------------------------------------

# dataset
tPB_NXT_count_derep<-t(PB_NXT_count_derep)
tPM_NXT_count_derep<-t(PM_NXT_count_derep)
tSP_NXT_count 

tPB_TSQ_count_derep<-t(PB_TSQ_count_derep) 
tPM_TSQ_count_derep <-t(PM_TSQ_count_derep) 
tSP_TSQ_count 

tPB_NXT_count_derep_rltv<-decostand(tPB_NXT_count_derep,"total",2)
tPM_NXT_count_derep_rltv<-decostand(tPM_NXT_count_derep,"total",2)
tSP_NXT_count_rltv<-decostand(tSP_NXT_count,"total",2)

tPB_TSQ_count_derep_rltv<-decostand(tPB_TSQ_count_derep,"total",2)
tPM_TSQ_count_derep_rltv<-decostand(tPM_TSQ_count_derep,"total",2)
tSP_TSQ_count_rltv<-decostand(tSP_TSQ_count,"total",2)


#function example
#diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
#fisher.alpha(x, MARGIN = 1, se = FALSE, ...)

PB_NXT_shanon<-diversity(tPB_NXT_count_derep_rltv, index = "shannon", MARGIN = 1, base = exp(1))
PM_NXT_shanon<-diversity(tPM_NXT_count_derep_rltv, index = "shannon", MARGIN = 1, base = exp(1))
SP_NXT_shanon<-diversity(tSP_NXT_count_rltv, index = "shannon", MARGIN = 1, base = exp(1))
PB_TSQ_shanon<-diversity(tPB_TSQ_count_derep_rltv, index = "shannon", MARGIN = 1, base = exp(1))
PM_TSQ_shanon<-diversity(tPM_TSQ_count_derep_rltv, index = "shannon", MARGIN = 1, base = exp(1))
SP_TSQ_shanon<-diversity(tSP_TSQ_count_rltv, index = "shannon", MARGIN = 1, base = exp(1))


shanon <- rep(0,48)
shanon[1:8] <- PB_NXT_shanon
shanon[9:16] <- PM_NXT_shanon        
shanon[17:24] <- SP_NXT_shanon
shanon[25:32] <- PB_TSQ_shanon
shanon[33:40] <- PM_TSQ_shanon
shanon[41:48] <- SP_TSQ_shanon

Extract<-c(rep("PB",8),
         rep("PM",8),
         rep("SP",8),
         rep("PB",8),
         rep("PM",8),
         rep("SP",8))
Library<-c(rep("NXT",24),
           rep("TSQ",24))
Sampleset<-c(rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2),
             rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2),
             rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2),
             rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2),
             rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2),
             rep("1", 2),rep("2", 2),rep("3", 2),rep("4", 2)
             )
Soilcore<-c(rep(seq(1,8,1),6))

shanon <- data.frame(shanon)
shanon<-cbind(Extract,Library,Sampleset,Soilcore,shanon)

NormalDistribution <- function(i){
  
  st<- shapiro.test(i)
  print (st)
  # Shapiro-Wilk normality test, p与查表中的W alpha 比较. p< W的可能性>0.05,即在0.05水平，p >0.05 就认为数据不符合正态分布 p值小于0.05，数据为正态分布
  
  kt<- ks.test(i, "pnorm", mean = mean(i), sd = sqrt(var(i))) 
  # This test is used for big data set>5000
  #Kolmogorov-Smirnov检验需要三个输入变量，及数据本身、均值及标准差
  # p与查表中的D alpha 比较. p > D 的可能性>0.05, 即，在0.05水平，p >0.05 就认为数据符合正态分布
  print(kt)
}

NormalDistribution(shanon[1:24,5])
NormalDistribution(shanon[25:48,5]) # maybe not fitted for normal distribution
NormalDistribution(shanon[,5])

# aov test for all the treatments 
# We dont consider Soilcore here, as Soilcore has no replicates and each of them are completely different
lmer_shanon_1 <- lmer(shanon ~ Extract*Library+(1|Sampleset),data=shanon)
summary(lmer_shanon_1)

lmer_shanon_2<-lm(shanon ~ Extract*Library,data=shanon)

anova(lmer_shanon_1,lmer_shanon_2)

summary(lmer_shanon_1) # no effect
anova(lmer_shanon_1)

# aov test for extracts within same library preparation method
summary(lm(shanon~Extract,data=shanon[1:24,])) # no effect
summary(lm(shanon~Extract,data=shanon[25:48,])) # no effect
anova(lm(shanon~Extract,data=shanon[1:24,])) # no effect
anova(lm(shanon~Extract,data=shanon[25:48,]))# no effect

# aov test for library preparation method within same extration method
PB_shanon<-subset(shanon,Extract=="PB")
PM_shanon<-subset(shanon,Extract=="PM")
SP_shanon<-subset(shanon,Extract=="SP")

summary(lm(shanon~Library,data=PB_shanon)) # no effect
summary(lm(shanon~Library,data=PM_shanon)) # no effect
summary(lm(shanon~Library,data=SP_shanon)) # no effect
anova(lm(shanon~Library,data=PB_shanon)) # no effect
anova(lm(shanon~Library,data=PM_shanon)) # no effect
anova(lm(shanon~Library,data=SP_shanon)) # no effect

LMER_TSQ_1<-lmer(shanon~Extract+(1|Sampleset),data=shanon[25:48,])
LMER_TSQ_2<-lm(shanon~Extract,data=shanon[25:48,])
anova(LMER_TSQ_1,LMER_TSQ_2)

anova(LMER_TSQ_2)
summary(LMER_TSQ_2)
plot(LMER_TSQ_2)


########### beta diversity ###################
MMB_merge_count_soilcore_rltv<-decostand(MMB_merge_count_soilcore,"total",2)
names(MMB_merge_count_soilcore_rltv)<-paste(Extract,sep = "_", Library)

plot(capscale(betad~Extract+Library+Sampleset+Soilcore,dist="Jaccard"))

NumToChr <- function(x) {if(is.numeric(x)) 
{as.character(paste(x))} else {x}} 

Soilcore<- NumToChr(Soilcore)

betad<-betadiver(t(MMB_merge_count_soilcore_rltv),"z")
adonis(betad~Sampleset,method="bray",
       #strata=Extract, 
       perm=200) # effect!!!
adonis(betad~Extract,method="bray",strata=Library, perm=200) #effect
adonis(betad~Library,method="bray",strata=Extract, perm=200) #effect


mod<-betadisper(betad,Sampleset)

plot(mod)

boxplot(mod)



betad<-betadiver(t(MMB_merge_count_soilcore_rltv[,1:24]),"z")
adonis(betad~Extract[1:24],method="bray", perm=200)
mod<-betadisper(betad,Extract[1:24])
plot(mod)
boxplot(mod)

betad<-betadiver(t(MMB_merge_count_soilcore_rltv[,25:48]),"z")
adonis(betad~Extract[25:48],method="Jaccard", perm=999)
mod<-betadisper(betad,Extract[25:48])
plot(mod)
boxplot(mod)

###################################### absense/Presence dataset ############################

######## alpha diversity ###############

# dataset
tPB_NXT_count_derep<-t(PB_NXT_count_derep)
tPM_NXT_count_derep<-t(PM_NXT_count_derep)
tSP_NXT_count<-t(SP_NXT_count)

tPB_TSQ_count_derep<-t(PB_TSQ_count_derep) 
tPM_TSQ_count_derep <-t(PM_TSQ_count_derep) 
tSP_TSQ_count <-t(SP_TSQ_count)

tPB_NXT_count_derep_rafy<-decostand(tPB_NXT_count_derep,"pa",2)
tPM_NXT_count_derep_rafy<-decostand(tPM_NXT_count_derep,"pa",2)
tSP_NXT_count_rafy<-decostand(tSP_NXT_count,"pa",2)

tPB_TSQ_count_derep_rafy<-decostand(tPB_TSQ_count_derep,"pa",2)
tPM_TSQ_count_derep_rafy<-decostand(tPM_TSQ_count_derep,"pa",2)
tSP_TSQ_count_rafy<-decostand(tSP_TSQ_count,"pa",2)

####### SHANON #######
# function example
# diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
# fisher.alpha(x, MARGIN = 1, se = FALSE, ...)
# shanon diversity is not suitable for incidence dataset, as it considered abundance input!

PB_NXT_shanon_rafy<-diversity(tPB_NXT_count_derep_rafy, index = "shannon", MARGIN = 1, base = exp(1))
PM_NXT_shanon_rafy<-diversity(tPM_NXT_count_derep_rafy, index = "shannon", MARGIN = 1, base = exp(1))
SP_NXT_shanon_rafy<-diversity(tSP_NXT_count_rafy, index = "shannon", MARGIN = 1, base = exp(1))
PB_TSQ_shanon_rafy<-diversity(tPB_TSQ_count_derep_rafy, index = "shannon", MARGIN = 1, base = exp(1))
PM_TSQ_shanon_rafy<-diversity(tPM_TSQ_count_derep_rafy, index = "shannon", MARGIN = 1, base = exp(1))
SP_TSQ_shanon_rafy<-diversity(tSP_TSQ_count_rafy, index = "shannon", MARGIN = 1, base = exp(1))


shanon_rafy <- rep(0,48)
shanon_rafy[1:8] <- PB_NXT_shanon_rafy
shanon_rafy[9:16] <- PM_NXT_shanon_rafy        
shanon_rafy[17:24] <- SP_NXT_shanon_rafy
shanon_rafy[25:32] <- PB_TSQ_shanon_rafy
shanon_rafy[33:40] <- PM_TSQ_shanon_rafy
shanon_rafy[41:48] <- SP_TSQ_shanon_rafy

shanon_rafy <- data.frame(shanon_rafy)
shanon_rafy <-cbind(Extract,Library,Sampleset,Soilcore,shanon_rafy)

NormalDistribution(shanon_rafy[1:24,5]) # maybe not fitted
NormalDistribution(shanon_rafy[25:48,5]) # maybe not fitted
NormalDistribution(shanon_rafy[,5]) # maybe not fitted


# aov test for all the treatments
lmer_shanon_rafy_1 <- lmer(shanon_rafy~Extract*Library+(1|Sampleset),data=shanon_rafy)
summary(lmer_shanon_rafy_1)

lmer_shanon_rafy_2<-lm(shanon_rafy~Extract*Library,data=shanon_rafy)
lmer_shanon_rafy_2

anova(lmer_shanon_rafy_1,lmer_shanon_rafy_2)

summary(lmer_shanon_rafy_2)

anova(lmer_shanon_rafy_2)


# aov test for extracts within same library preparation method
summary(lm(shanon_rafy~Extract,data=shanon_rafy[1:24,])) # no effect
summary(lm(shanon_rafy~Extract,data=shanon_rafy[25:48,])) # effect difficult to explain !!!!
anova(lm(shanon_rafy~Extract,data=shanon_rafy[1:24,])) # no effect
anova(lm(shanon_rafy~Extract,data=shanon_rafy[25:48,])) # effect difficult to explain !!!!

LMER_TSQ_rafy_1<-lmer(shanon_rafy~Extract+(1|Sampleset),data=shanon_rafy[25:48,])
LMER_TSQ_rafy_2<-lm(shanon_rafy~Extract,data=shanon_rafy[25:48,])
anova(LMER_TSQ_rafy_1,LMER_TSQ_rafy_2)

anova(LMER_TSQ_rafy_2)
summary(LMER_TSQ_rafy_2)
plot(LMER_TSQ_rafy_2)






# aov test for library preparation method within same extration method
PB_shanon_rafy<-subset(shanon_rafy,Extract=="PB")
PM_shanon_rafy<-subset(shanon_rafy,Extract=="PM")
SP_shanon_rafy<-subset(shanon_rafy,Extract=="SP")

summary(lm(shanon_rafy~Library,data=PB_shanon_rafy)) # no
summary(lm(shanon_rafy~Library,data=PM_shanon_rafy)) # no
summary(lm(shanon_rafy~Library,data=SP_shanon_rafy)) # no
anova(lm(shanon_rafy~Library,data=PB_shanon_rafy)) # no
anova(lm(shanon_rafy~Library,data=PM_shanon_rafy)) # no
anova(lm(shanon_rafy~Library,data=SP_shanon_rafy)) # no

#### Chao2/ICE/Jackkinfe #####

library(fossil)

Chao2_PB_TSQ<-chao2(tPB_TSQ_count_derep_rafy,taxa.row=F)
Chao2_PM_TSQ<-chao2(tPM_TSQ_count_derep_rafy,taxa.row=F)
Chao2_SP_TSQ<-chao2(tSP_TSQ_count_rafy,taxa.row=F)

alpha_fossil_PB_NXT<-spp.est(t(tPB_NXT_count_derep_rafy),abund=F)
alpha_fossil_PM_NXT<-spp.est(t(tPM_NXT_count_derep_rafy),abund=F)
alpha_fossil_SP_NXT<-spp.est(t(tSP_NXT_count_rafy),abund=F) 
alpha_fossil_PB_TSQ<-spp.est(t(tPB_TSQ_count_derep_rafy),abund=F)
alpha_fossil_PM_TSQ<-spp.est(t(tPM_TSQ_count_derep_rafy),abund=F)
alpha_fossil_SP_TSQ<-spp.est(t(tSP_TSQ_count_rafy),abund=F)

Chao2 <- rep(0,48)
Chao2[1:8] <- alpha_fossil_PB_NXT[,5]
Chao2[9:16] <- alpha_fossil_PM_NXT[,5]    
Chao2[17:24] <- alpha_fossil_SP_NXT[,5]
Chao2[25:32] <- alpha_fossil_PB_TSQ[,5]
Chao2[33:40] <- alpha_fossil_PM_TSQ[,5]
Chao2[41:48] <- alpha_fossil_SP_TSQ[,5]

alpha_fossil_Chao2 <-shanon_rafy
alpha_fossil_Chao2[,5]<-Chao2

ICE <- rep(0,48)
ICE[1:8] <- alpha_fossil_PB_NXT[,8]
ICE[9:16] <- alpha_fossil_PM_NXT[,8]    
ICE[17:24] <- alpha_fossil_SP_NXT[,8]
ICE[25:32] <- alpha_fossil_PB_TSQ[,8]
ICE[33:40] <- alpha_fossil_PM_TSQ[,8]
ICE[41:48] <- alpha_fossil_SP_TSQ[,8]

alpha_fossil_ICE <-shanon_rafy
alpha_fossil_ICE[,5]<-ICE

Jack1 <- rep(0,48)
Jack1[1:8] <- alpha_fossil_PB_NXT[,11]
Jack1[9:16] <- alpha_fossil_PM_NXT[,11]    
Jack1[17:24] <- alpha_fossil_SP_NXT[,11]
Jack1[25:32] <- alpha_fossil_PB_TSQ[,11]
Jack1[33:40] <- alpha_fossil_PM_TSQ[,11]
Jack1[41:48] <- alpha_fossil_SP_TSQ[,11]

alpha_fossil_Jack1 <-shanon_rafy
alpha_fossil_Jack1[,5]<-Jack1


# aov test
# Anova test is not propiate for incidence based diversity indexes

#### graph for CHAO2/ICE/Jackkinfe2 NXT####
par(cex = 0.8)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)

par(mfrow=c(2,3))
plot(alpha_fossil_PB_NXT[,5]~alpha_fossil_PB_NXT[,1],
     ylim=c(10,100),
     xlim=c(1,8),
     xaxt="n",
     col="blue",
     type="l") #Chao2
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,6],rev(alpha_fossil_PB_NXT[,7])),border=NA,col=rgb(0, 0, 1,0.5))#blue
mtext("Nextera", font=2,side=2,line=2.5,cex =1)
mtext("Chao2",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_NXT[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     type="l") #Chao2
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,6],rev(alpha_fossil_PM_NXT[,7])),border=NA,col=rgb(1, 0, 0,0.5))#red

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     type="l") #Chao2
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,6],rev(alpha_fossil_SP_NXT[,7])),border=NA,col=rgb(1, 1, 0,0.5))#orange

plot(alpha_fossil_PB_NXT[,8],
     ylim=c(10,100),
     xaxt="n",
     yaxt="n",
     col="blue",
     type="l") #ICE
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,9],rev(alpha_fossil_PB_NXT[,10])),border=NA,col=rgb(0, 0, 1,0.5))#red
mtext("ICE",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_NXT[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     type="l") #ICE
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,9],rev(alpha_fossil_PM_NXT[,10])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     type="l") #ICE
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,9],rev(alpha_fossil_SP_NXT[,10])),border=NA,col=rgb(1, 1, 0,0.5))

plot(alpha_fossil_PB_NXT[,11],
     ylim=c(10,100),
     col="blue",
     xaxt="n",
     yaxt="n",
     type="l") #Jacknife
polygon(c(alpha_fossil_PB_NXT[,1],rev(alpha_fossil_PB_NXT[,1])),
        c(alpha_fossil_PB_NXT[,12],rev(alpha_fossil_PB_NXT[,13])),border=NA,col=rgb(0, 0, 1,0.5))
mtext("Jackknife1",font=2,side=3,line=1.5,cex =1)

par(new=TRUE)  
plot(alpha_fossil_PM_NXT[,11],
     ylim=c(10,100),
     axes=F,
     xlab="",
     xaxt="n",
     yaxt="n",
     col="red",
     type="l") #Jacknife
polygon(c(alpha_fossil_PM_NXT[,1],rev(alpha_fossil_PM_NXT[,1])),
        c(alpha_fossil_PM_NXT[,12],rev(alpha_fossil_PM_NXT[,13])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_NXT[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     yaxt="n",
     xlab="",
     col="orange",
     type="l") #Jacknife
polygon(c(alpha_fossil_SP_NXT[,1],rev(alpha_fossil_SP_NXT[,1])),
        c(alpha_fossil_SP_NXT[,12],rev(alpha_fossil_SP_NXT[,13])),border=NA,col=rgb(1, 1, 0,0.5))

legend(3,20, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black",
       lty=c(1, 1, 1),
       bty="n",cex=1)

####graph for CHAO2/ICE/Jackkinfe2 TSQ #################

plot(alpha_fossil_PB_TSQ[,5]~alpha_fossil_PB_TSQ[,1],
     ylim=c(10,100),
     xlim=c(1,8),
     xlab="Chao2",
     col="blue",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,6],rev(alpha_fossil_PB_TSQ[,7])),border=NA,col=rgb(0, 0, 1,0.5))#blue
mtext("TruSeq", font=2,side=2,line=2.5,cex =1)

par(new=TRUE)
plot(alpha_fossil_PM_TSQ[,5],
     ylim=c(10,100),
     axes=F,
     xlab="",
     col="red",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,6],rev(alpha_fossil_PM_TSQ[,7])),border=NA,col=rgb(1, 0, 0,0.5))#red

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,5],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     lty=2,
     type="l") #Chao2
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,6],rev(alpha_fossil_SP_TSQ[,7])),border=NA,col=rgb(1, 1, 0,0.5))#orange

plot(alpha_fossil_PB_TSQ[,8],
     ylim=c(10,100),
     xlab="ICE",
     yaxt="n",
     col="blue",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,9],rev(alpha_fossil_PB_TSQ[,10])),border=NA,col=rgb(0, 0, 1,0.5))#red

par(new=TRUE)
plot(alpha_fossil_PM_TSQ[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,9],rev(alpha_fossil_PM_TSQ[,10])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,8],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="orange",
     lty=2,
     type="l") #ICE
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,9],rev(alpha_fossil_SP_TSQ[,10])),border=NA,col=rgb(1, 1, 0,0.5))

plot(alpha_fossil_PB_TSQ[,11],
     ylim=c(10,100),
     col="blue",
     yaxt="n",
     lty=2,
     type="l") #Jacknife
polygon(c(alpha_fossil_PB_TSQ[,1],rev(alpha_fossil_PB_TSQ[,1])),
        c(alpha_fossil_PB_TSQ[,12],rev(alpha_fossil_PB_TSQ[,13])),border=NA,col=rgb(0, 0, 1,0.5))

par(new=TRUE)  
plot(alpha_fossil_PM_TSQ[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     col="red",
     lty=2,
     type="l") #Jacknife
polygon(c(alpha_fossil_PM_TSQ[,1],rev(alpha_fossil_PM_TSQ[,1])),
        c(alpha_fossil_PM_TSQ[,12],rev(alpha_fossil_PM_TSQ[,13])),border=NA,col=rgb(1, 0, 0,0.5))

par(new=TRUE)
plot(alpha_fossil_SP_TSQ[,11],
     ylim=c(10,100),
     axes=F,
     xaxt="n",
     xlab="",
     lty=2,
     col="orange",
     type="l") #Jacknife
polygon(c(alpha_fossil_SP_TSQ[,1],rev(alpha_fossil_SP_TSQ[,1])),
        c(alpha_fossil_SP_TSQ[,12],rev(alpha_fossil_SP_TSQ[,13])),border=NA,col=rgb(1, 1, 0,0.5))

legend(3,20, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black",
       lty=c(2, 2, 2),
       bty="n",cex=1)



########### beta diversity ###################
MMB_merge_count_soilcore_rafy<-decostand(MMB_merge_count_soilcore,"pa",2)

NumToChr <- function(x) {if(is.numeric(x)) 
{as.character(paste(x))} else {x}} 

Soilcore<- NumToChr(Soilcore)

betad_rafy<-betadiver(t(MMB_merge_count_soilcore_rafy),"z")

beta.dist<-vegdist(MMB_merge_count_soilcore_rafy,"binomial",binary=TRUE) # define the dissimilarity indexes
adonis(betad_rafy~Extract, method="Jacaard", strata=Library,perm=200)
adonis(betad_rafy~Library, method="Jacaard", strata=Extract,perm=200)
#The function also finds indices for presence/ absence data by setting binary = TRUE

Extract_Lib <- paste(Extract,Library)
par(mfrow=c(1,1))
par(mar = c(4, 4, 2, 0), oma = c(1, 1, 1, 0.5))
mod_rafy<-betadisper(betad_rafy,Extract_Lib)
# method from vegan tutorial
anova(mod_rafy)
TukeyHSD(mod_rafy)
# how to extract % for each axes https://stat.ethz.ch/pipermail/r-sig-ecology/2011-June/002183.html
eig<-eigenvals(mod_rafy)
axes_percent<-eig/sum(eig)
plot(mod_rafy,
     col= c("blue","blue","red","red","orange","orange"),
     main="Multivariate Dispersion",
     pch = c(15,22,17,24,18,23),
     segments = FALSE, # control grey line
     hull = FALSE, ellipse = FALSE,label = F,lwd=0.5)
ordiellipse(mod_rafy, Extract_Lib, 
            kind="se", 
            conf=0.95, 
            lwd=1, 
            lty = c(1,2,1,2,1,2),
            col= c("blue","blue","red","red","orange","orange")
            )

legend(-0.5,0.3, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black", 
       title="Nextera",
       pch = c(15, 17, 18),
       lty=0,
       bty="n",lwd=1,cex=0.75)
legend(-0.5,0.19, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       title="TruSeq",
       text.col = "black", 
       pch = c(22, 24, 23),
       lty=0,
       bty="n",lwd=1,cex=0.75)


betad_NXT<-betadiver(t(MMB_merge_count_soilcore_rafy[,1:24]),"z")
adonis(betad_NXT~Extract[1:24],method="beta.dist", perm=999)
mod_NXT<-betadisper(betad_NXT,Extract[1:24])
plot(mod)
boxplot(mod)

betad_TSQ<-betadiver(t(MMB_merge_count_soilcore_rafy[,25:48]),"z")
adonis(betad_TSQ~Extract[25:48],method="beta.dist", perm=200)
mod_TSQ<-betadisper(betad_TSQ,Extract[25:48])
par(mfrow=c(1,2))
plot(mod_NXT,
     col= c("blue","red","orange"),
     main="Multivariate Dispersion for Nextera",
     pch = c(15,17,18),
     hull = FALSE, ellipse = TRUE,label = F,lwd=1
)
plot(mod_TSQ,
     col= c("blue","red","orange"),
     main="Multivariate Dispersion for TruSeq",
     pch = c(15,17,18),
     hull = FALSE, ellipse = TRUE,label = F,lwd=1
     )

legend(-0.5,0.25, c("Phosphate buffer", "Powermax", "DNA soup"), 
       col = c( "blue","red","orange" ),
       text.col = "black",
       
       pch = c(15,17,18),
       bty="n",cex=1)

boxplot(mod)










 
#### PCA       #########

library("ade4")
names(OTU_merge_minrep)
MMB_merge_rltv

OTU_pca <-MMB_merge_count_soilcore_rltv
OTU_pca <-t(OTU_pca)
row.names(OTU_pca)<-paste(row.names(OTU_pca),"_",rep(seq(1,8,1),6))

boxplot(OTU_pca,las=2)
boxplot(data.frame(scale(OTU_pca)),las=2)

#OTU_transpca<-log(OTU_pca+1)

#Perform PCA
row1<-c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8)) #row
pca1<-dudi.pca(OTU_pca,
               row.w=row1,
               scale = TRUE,
               scannf = F,
               nf=3)
#How many axes shall we keep? This is used to judge the important components.

barplot(pca1$eig)
explain<-pca1$eig/sum(pca1$eig)
barplot(explain)
#Correlation circle
s.corcircle(pca1$co)

#Interpretation
#plot()
scatter(pca1)
gcol = c("red","blue","green","orange","purple","grey")
Frow1<-as.factor(c(rep("Phosphate_buffer_NXT",8),
                   rep("Powermax_NXT",8),
                   rep("DNA_soup_NXT",8),
                   rep("Phosphate_buffer_TSQ",8),
                   rep("Powermax_TSQ",8),
                   rep("DNA_soup_TSQ",8)))
par(mfrow=c(2,2))
s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 1, 
        yax = 2,
        addaxes=T,
        grid=F,
        include.origin = T,
        origin=c(0,0),# the fixed point for the origin axes
        clabel = 1)
s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 2, 
        yax = 3,
        grid=F,
        addaxes=T,
        origin=c(0,0),# the fixed point for the origin axes
        clabel = 1)

s.class(pca1$li,
        fac=Frow1,
        col = gcol, 
        xax = 1, 
        yax = 3,
        addaxes=T,
        origin=c(0,0),# the fixed point for the origin axes
        grid=F,
        clabel = 1)

s.label(pca1$li,clab=0.5,add.p=T)
s.arrow(pca1$c1,add.p=T)

plot3d(pca1$li,col=gcol,type="p",radius=0.05)

####BCA between site PCA ########
pca_trybca <- dudi.pca(OTU_pca, scan = FALSE, nf = 3)
bet1 <- bca(pca_trybca,Frow1, scan = FALSE, nf = 3)
par(mfrow=c(2,2))
s.class(bet1$ls, 
        Frow1, 
        #sub = "Between sites PCA", 
        xax = 1, 
        yax = 2,
        csub = 1.75)
par(new=F)
s.class(bet1$ls, 
        Frow1, 
        #sub = "Between sites PCA", 
        xax = 3, 
        yax = 2,
        csub = 1.75)
par(new=F)
s.class(bet1$ls, 
        Frow1, 
        #sub = "Between sites PCA", 
        xax = 1, 
        yax = 3,
        csub = 1.75)

par(mfrow=c(1,1))
plot(bet1)

##### PERMANOVA ####
library(vegan)
TREAT <- read.csv("/Users/JiayiQin/Dropbox/Paper_III_MMB/Data_structure.csv")
PERM <- MMB_merge_count[,2:105]
IntToNum <- function(x) {if(is.integer(x)) 
            {as.numeric(paste(x))} else {x}} 
PERM<- IntToNum(PERM)
PERM<-t(PERM)
adonis2(PERM~ Method+Library, data = TREAT)
