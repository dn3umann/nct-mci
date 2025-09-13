
# ---- CONFIGURATION ---- #

#Set up relative file paths (see README for file structure)
base.dir=file.path("Project Folder")
matlab.dir=file.path(base.dir,"MATLAB")
data.dir=file.path(base.dir,"Data")
results.dir=file.path(base.dir,"Results")
figures.dir=file.path(results.dir,"Figures")
  
#Load R packages
library(R.matlab)
library(stringr)
library(viridis)
library("lattice")
library(robustbase)
library(matrixStats)
library(vioplot)
library(pracma)
library(gridExtra)
library(fmsb)
library(RColorBrewer)
library(ggplot2)
library(viridisLite)
library(ggpubr)
library(ppcor)

#Load time series data and find number of subjects
file.names=list.files(file.path(data.dir,"TimeSeries"))
length(file.names) #total number of subjects with time series data



# ---- DATA VALIDATION ---- #

#Ensure that there is a complete set of time series data for each subject
for(i in 1:length(file.names)){
  x=read.table(file.names[i],header=F,sep="",dec=".")
  if (all(x==0)) { #find the files where the time series is all zeroes
    print(paste((file.names[i]),"has no data.")) }
  else {
    if (any(x==0)) { #find the files where the time series has zeroes
      print(paste((file.names[i]),"has bad data.")) }}
} #If any files are returned, remove them from the TimeSeries folder



# ---- DEMOGRAPHICS ---- #

#Demographic.Data.RData should have the following columns: ID, Age, Sex, Status, T1, one column per region

#Import data
load(file.path(data.dir,"Demographic.Data.RData"))
dem.data=as.data.frame(data)
dem.data=dem.data[order(dem.data$T1,decreasing=FALSE),]

#Find age distribution
print(paste("Mean Age (All Subjects):",mean(dem.data$Age)))
print(paste("Mean Age (HC):",mean(subset(dem.data,Status=="HC")$Age)))
print(paste("Mean Age (MCI):",mean(subset(dem.data,Status=="MCI")$Age)))

#Find sex distribution; females = 1, males = 0
print(paste("Percent Female (All Subjects):",
            (100*(sum(dem.data$Sex==1)/nrow(dem.data)))))
print(paste("Percent Female (HC):",
            (100*(sum(subset(dem.data,Status=="HC")$Sex==1))/nrow(subset(dem.data,Status=="HC")))))
print(paste("Percent Female (MCI):",
            (100*(sum(subset(dem.data,Status=="MCI")$Sex==1))/nrow(subset(dem.data,Status=="MCI")))))

#Find ICV distribution
print(paste("Mean ICV (All Subjects):",mean(dem.data$ICV)))
print(paste("Mean ICV (HC):",mean(subset(dem.data,Status=="HC")$ICV)))
print(paste("Mean ICV (MCI):",mean(subset(dem.data,Status=="MCI")$ICV)))

#Wilcoxon rank-sum to compare age and ICV between HC and MCI
wilcox.test(Age~Status,data=dem.data)
wilcox.test(ICV~Status,data=dem.data)

#Chi-square test to compare sex between HC and MCI
chisq.test(dem.data$Sex,dem.data$Status)



# ----TIME SERIES Z-SCORING ---- #

time.series.matrix=NULL #initialize the matrix

#Z-score all time series data at the subject-level
for(i in 1:length(file.names)){
  matrix.entry=as.matrix(read.table(file.names[i],header=F,sep="",dec=".")) #TS data for subject i
  mean1=mean(matrix.entry) #mean for subject i
  sd1=sd(matrix.entry) #SD for subject i
  time.series.matrix.sub.level.zscore=(matrix.entry-mean1)/sd1 #z-score for subject i
  time.series.matrix=rbind(time.series.matrix,t(time.series.matrix.sub.level.zscore)) #matrix containing z-scores for all subjects
}

#Print dimensions of the matrix to confirm correct number of regions, then save as a MATLAB file
print('Dimensions of Matrix');dim(time.series.matrix)
writeMat(con=file.path(matlab.dir,"TimeSeriesZScored.mat"),labpcexport=time.series.matrix)



# ---- EXECUTE BLOCKS 1-5 IN MATLAB ---- #



# ---- TIME SERIES LENGTH ---- #

#TS length differs between subjects and must be accounted for when calculating transition probability

#Create a matrix containing the TS lengths for each subject
lengthTS=NULL

for(i in 1:length(file.names)){
  subjectTS=read.table(file.names[i],header=F,sep="",dec=".") #TS data for subject i
  lengthTS=rbind(lengthTS,ncol(subjectTS)) #matrix containing TS lengths for all subjects
}

#Export matrix to MATLAB file
writeMat(con=file.path(matlab.dir,"LengthTS.mat"),labpcexport=lengthTS)



# ---- EXECUTE BLOCKS 6-7 IN MATLAB ---- #



# ---- IMPORT TRANSITION PROBABILITIES TO R ---- #

#Import transition probabilities and confirm that length matches the number of subjects
transition.prob=readMat(file.path(matlab.dir,"Transition.Prob.mat"))$tran.prob
length(transition.prob)
n=554 #n = number of subjects

#Create subject-wise transition probabilities
transition.prob.persubject=NULL
for(i in 1:n){
  transition.prob.persubject=rbind(transition.prob.persubject,as.numeric(transition.prob[[i]][[1]]))
}

#Dimensions should be n x k^2
dim(transition.prob.persubject)



# ---- LOAD THE STRUCTURAL CONNECTOME AND EXAMINE DATA ---- #

#SC.mat should be the structural connectome for the applicable parcellation (i.e., Schaefer 200) 

#Assign average SC data (diagonals equal to zero) to a matrix
SC.matrix=readMat(file.path(data.dir,"SC.mat"))$labpcexport

#Confirm matrix dimensions and diagonal values
dim(SC.matrix) #dimensions should be number of regions x number of regions
SC.matrix[1:5,1:5] #view first 5 rows and first 5 columns to confirm that the diagonals are set to zero

#Save the matrix as a MATLAB file
writeMat(con=file.path(matlab.dir,"SC.Matrix.diag0.mat"),labpcexport=SC.matrix)



# ---- IDENTIFY BRAIN STATES AND CREATE RADIAL PLOTS ---- #

#sch232_to_yeo.csv matches the Schaefer 232 parcellation to the Yeo 7 networks

#Load the region-to-network table and centroids
sch.232.yeo=as.numeric(read.table(file.path(data.dir,"sch232_to_yeo.csv"),sep=",",header=FALSE))
sch.200.yeo=sch.232.yeo[1:200] #omit the last 32 regions (subcortex) because they're not used
centroids.k=readMat(file.path(matlab.dir,"Centroids.k.mat"))$labpcexport
dim(centroids.k6) #dimensions should be number of subjects x k

#Create and store files for the radial plots
radial.plot.files=c("State.1.Radial.pdf",
                    "State.2.Radial.pdf",
                    "State.3.Radial.pdf",
                    "State.4.Radial.pdf",
                    "State.5.Radial.pdf",
                    "State.6.Radial.pdf")
radial.plots=file.path(figures.dir,radial.plot.files)

#Find which state corresponds to each cluster
k=6 #number of clusters
for(clusternumber in 1:k){cat(clusternumber)
  centroids.cluster.k=centroids.k[,clusternumber] #isolate centroids for cluster k
  pos.centroids.cluster.k=replace(centroids.cluster.k,which(centroids.cluster.k<0),0) #isolate positive centroids
  neg.centroids.cluster.k=replace(centroids.cluster.k,which(centroids.cluster.k>0),0) #isolate negative centroids
  
  positive.centroid.cluster.k=NULL
  negative.centroid.cluster.k=NULL
  
  for(i in c(1,2,3,4,5,6,7)){ #assign each centroid to a network
    positive.centroid.cluster.k=rbind(positive.centroid.cluster.k,
                                      mean(pos.centroids.cluster.k[which(sch.200.yeo==i)]))
    negative.centroid.cluster.k=rbind(negative.centroid.cluster.k,
                                      mean(neg.centroids.cluster.k[which(sch.200.yeo==i)]))
  }
  
  #Find the highest amplitude positive centroid and highest amplitude negative centroid
  cbind(positive.centroid.cluster.k/max(positive.centroid.cluster.k),
        negative.centroid.cluster.k/(max(abs(negative.centroid.cluster.k))))
  
  #Find the highest amplitude centroid among both the positive and negative centroids
  cbind(positive.centroid.cluster.k,negative.centroid.cluster.k)
  maxvalue=max(max(positive.centroid.cluster.k),max(abs(negative.centroid.cluster.k)));maxvalue
  
  #Normalize centroids so that each is scaled between 0 and 1
  positive.state=positive.centroid.cluster.k/maxvalue
  negative.state=negative.centroid.cluster.k/maxvalue
  
  maxvalue.chart=1 #upper bound of the plot
  
  #Create a data frame with the following parameters
    #columns: correspond to each network
    #row 1: defines the upper bound of the chart (constant for all clusters; essentially the scale)
    #row 2: defines lower bound of the chart (zero)
    #row 3: average value of the positive centroids for each network
    #row 4: average value of the negative centroids for each network
  data.state=as.data.frame(rbind(as.numeric(positive.state), #row 3
                           as.numeric(abs(negative.state))),ncol=7) #row 4 
  colnames(data.state)=c("Visual","Somatomotor","Dorsal Attention","Salient Ventral Attention",
                         "Limbic","Control","Default") #columns
  data.state=rbind(rep(maxvalue.chart,k),rep(0,k),data.state) #row 1 and row 2
  
  
  #Determine which network corresponds to each cluster
  most.positive=colnames(data.state)[which.max(data[3,])] #find highest amplitude positive network
  most.negative=colnames(data.state)[which.max(data[4,])] #find highest amplitude negative network
  if(max(data[3,])>max(data.state[4,])){ #highest overall amplitude network among both positive and negative
    highest.value=paste(most.positive,"(+)")
  } else {
    highest.value=paste(most.negative,"(-)")
  }
  
  #Print out the names of the networks for each state
  network.results=NULL
  network.results=c(results,paste("\nCluster",clusternumber,"Network:",highest.value,
                                  "\nCluster",clusternumber,"High Amplitude Network:",most.positive,
                                  "\nCluster",clusternumber,"Low Amplitude Network:",most.negative,"\n",collapse="\n"))
  cat(network.results)
  
  
  #Create a radial plot for each state
  pdf(file=radial.plot.files[[clusternumber]],width=18,height=12)
  par(mar=c(5,5,5,5));par(mfrow=c(1,1))
  newradarchart(data=data.state,
                color=c("#FFC300","purple"), #choose colors for high and low amplitude networks
                vlabels=c("VIS","SMN","DAN","VAN","LIM","FPN","DMN"),vlcex=3, #list names of each network
                caxislabels=c(paste(0,"  "),
                              paste(0.25,"  "),
                              paste(0.5,"  "),
                              paste(0.75,"  "),
                              paste(1,"  ")),
                calcex=2)
  dev.off()
} 



# ---- EXECUTE BLOCKS 8-9 IN MATLAB ---- #



# ---- FIND BEST T AND CREATE A GRAPH FOR CORRELATION ---- #

#Load data from MATLAB
E.full.T=readMat(file.path(results.dir,"E.Full.T.mat"))$E.Full.T

#Run a loop to find the correlation coefficients and p-values for each T
cor.energy=NULL
pvals.energy=NULL

for(t in 1:20){ #there are 20 different values of T being tested
  cor.entry=cor.test(rank(colMeans(transition.prob.persubject)),
                     rank(E.full.T[t,]),method="spearman")
  cor.energy=rbind(cor.energy,cor.entry$estimate) #save each correlation coefficient
  pvals.energy=rbind(pvals.energy,cor.entry$p.value) #save p-values for each coefficient
}

#Print correlation coefficients and p-values
print(cat("Correlation Coefficients:",head(cor.energy)))
print(cat("p-values:",head(pvals.energy)))

#Find which T corresponds to the highest amplitude of correlation
best.t=seq(0.001,10,0.5)[which.max(abs(cor.energy))]
print(paste("Best T:",best.t))

#Graph the amplitudes of correlation for each value of T
pdf(file=file.path(figures.dir,"T.Correlation.Graph.pdf"),width=14,height=8)
plot(seq(0.001,10,0.5),abs(cor.energy),type="b",
     ylab="Amplitude of the correlation",xaxt="none",xlab="t")
axis(1,seq(0.001,10,0.5)) #range of T values
dev.off()



# ---- EXECUTE BLOCKS 10-11 IN MATLAB ---- #



# ---- COMPUTE ENTROPY ---- #

n=554 #number of subjects
regions=200 #number of regions

#Initialize values for the loop
entropy=matrix(NA,ncol=regions,nrow=n)
k=1
for(i in 1:n){cat(k)
  timeseries=t(read.table(file.names[i],header=F,sep="",dec="."));dim(timeseries)
  for(region in 1:(dim(timeseries)[2])){
    timeseries1=timeseries[,region]
    entropy[k,region]=sample_entropy(timeseries1[-c(1:10)],edim=2) #calculate per subject, per region entropy
  }
  k=k+1
}
dim(entropy)

#Save entropy data
save(entropy,file=file.path(results.dir,"Entropy.RData"))



# ---- IMPORT AND PAIR PET DATA ---- #

#Import PET data
raw.plaque=read.csv(file.path(data.dir,"Plaque.csv"))
names(raw.plaque)[1]="ID" #rename first column to "ID"
raw.tau=read.csv(file.path(data.dir,"Tau.csv"))
names(raw.tau)[1]="ID" #rename new first column to "ID"
dem.data=as.data.frame(dem.data)

#Only use data for subjects that were included for TE
tau=raw.tau[raw.tau$ID %in% dem.data$ID,] #only use subjects that were used for the other analyses
tau=tau[match(dem.data$ID,tau$ID),]
plaque=raw.plaque[raw.plaque$ID %in% dem.data$ID,] #only use subjects that were used for the other analyses
plaque=plaque[match(dem.data$ID,plaque$ID),]
dim(tau);dim(plaque) #should match number of subjects and numbers of regions

#Save the cleaned PET data
save(tau,file=file.path(data.dir,"Tau.Clean.RData"))
save(plaque,file=file.path(data.dir,"Plaque.Clean.RData"))



# ---- CREATE STATE-WISE TRANSITION ENERGY HEAT MAP ---- #

#Import data and initialize values for the loop
minimumTE=readMat(file.path(results.dir,"State.Minimum.TE.mat"))$minimumTE
lm.entry=NULL #variable to store lm results for i
lm.state.TE.estimate=NULL #variable to store each lm estimate
lm.state.TE.pvalue=NULL #variable to store each lm p-value

for(i in 1:36){ #compare TE for 6x6 clusters
  lm.entry=summary(lm(minimumcontrolenergy[,i]~Status+Sex+Age+ICV,data=dem.data))
  lm.state.TE.estimate=rbind(lm.state.TE.estimate,lm.entry$coefficients[2,1])
  lm.state.TE.pvalue=rbind(lm.state.TE.pvalue,lm.entry$coefficients[2,4])
}

#Correct the p-values using BH
adjusted.lm.state.TE.pvalue=p.adjust(lm.state.TE.pvalue,n=36,method="BH");adjusted.lm.state.TE.pvalue


#Input the networks from line 235 in the order they were given
networks=c("DAN (+)","DMN (+)","DMN (-)",
           "DAN (-)", "VIS (-)","VIS (+)")

#Input the preferred order of networks to be displayed on the plot
ordered.networks=c("DAN (+)","DAN (-)","DMN (+)",
                   "DMN (-)","VIS (+)","VIS (-)")

#Create a matrix with the lm estimates
matrix=matrix(lm.state.TE.estimate,6,6)

#Match networks to ordered networks (i.e., DAN- listed fourth in networks but second in ordered networks)
ordered.matrix=matrix[c(1,4,2,3,6,5),c(1,4,2,3,6,5)]

#Create and save the state-wise TE heat map using the linear model estimates
png(file=file.path(figures.dir,"Statewise.Comparison.png"),width=800,height=600)
plot=levelplot(t(ordered.matrix),
               col.regions=colorRampPalette(rev(brewer.pal(11,"RdBu"))),
               colorkey=list(space="left"), #location of color key
               at=seq(-71.6,71.6,length.out=100),
               main="Comparison of State-Wise Transition Energy using lm Estimates", #positive means MCI > HC
               xlab="",ylab="",
               scales=list(x=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2),
                           y=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2)));plot #omit scales argument if networks are not yet known


#Create and plot TE for HC only
average.state.TE.HC=colMeans(minimumTE[which(dem.data$Status=="HC"),]) #isolate HC data only
matrix.HC=matrix(average.state.TE.HC,6,6) #create matrix for HC data only
ordered.matrix.HC=matrix[c(1,4,2,3,6,5),c(1,4,2,3,6,5)] #order should match ordered.matrix
png(file=file.path(figures.dir,"Statewise.HC.png"),width=800,height=600)
plot=levelplot(t(ordered.matrix.HC),
               col.regions=colorRampPalette(plasma(200)), #color palette and how many sub-colors
               colorkey=list(space="left"), #location of color key
               main="State-Wise Transition Energy for HC",
               xlab="",ylab="",
               scales=list(x=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2),
                           y=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2)));plot


#Create and plot TE for MCI only
average.state.TE.MCI=colMeans(minimumTE[which(dem.data$Status=="MCI"),]) #isolate HC data only
matrix.MCI=matrix(average.state.TE.MCI,6,6) #create matrix for HC data only
ordered.matrix.HC=matrix[c(1,4,2,3,6,5),c(1,4,2,3,6,5)] #order should match ordered.matrix
png(file=file.path(figures.dir,"Statewise.MCI.png"),width=800,height=600)
plot=levelplot(t(ordered.matrix.MCI),
               col.regions=colorRampPalette(plasma(200)), #color palette and how many sub-colors
               colorkey=list(space="left"), #location of color key
               main="State-Wise Transition Energy for HC",
               xlab="",ylab="",
               scales=list(x=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2),
                           y=list(at=seq(1,6),labels=ordered.networks,cex=1.5,font=2)));plot



# ---- ANALYZE REGIONAL TRANSITION ENERGIES ---- #

#Import and initialize data for the loop
regional.TE=readMat(file.path(results.dir,"Regional.Minimum.TE.mat"))$regionalTE[,1,]
lm.estimates=NULL #variable to store lm estimates for each region
lm.pvalue=NULL #variable to store p-value for each region

#Assign reference values to status and sex
dem.data$Status=relevel(dem.data$Status,ref="HC") #estimate > 0 indicates TE higher for MCI (MCI - HC)
dem.data$Sex=relevel(dem.data$Sex,ref="1") #females are coded as 1; estimate > 0 indicates higher TE for males (males - females)

#Compare TE for each region using the linear model adjusted for sex, age, and ICV
for(i in 1:200) { #adjust bounds to match number of regions (and columns storing the data)
  lm.entry=summary(lm(regional.TE[,i]~Status+Sex+Age+ICV,data=dem.data))
  lm.estimates=rbind(lm.estimates,lm.entry$coefficients[2,1]) #isolate coefficients
  lm.pvalue=rbind(lm.pvalue,lm.entry$coefficients[2,4]) #isolate p-value for status
}

#Create MATLAB files that will be used to make the brain plots for comparisons
lm.estimates.status=lm.estimates[,1]
writeMat(con=file.path(data.dir,"TE.Estimates.Status.mat"),labpcexport=lm.estimates.status)

#Create MATLAB files containing the regional TE for HC and MCI separately
regional.HC=colMeans(regional.TE[which(dem.data$Status=="HC"),])
writeMat(con=file.path(data.dir,"Regional.TE.HC.mat"),labpcexport=regional.HC)
regional.MCI=colMeans(regional.TE[which(dem.data$Status=="MCI"),])
writeMat(con=file.path(data.dir,"Regional.TE.MCI.mat"),labpcexport=regional.MCI)

#Identify regions with significantly different TE between HC and MCI after BH correction
significant.regions=which(p.adjust(lm.pvalue,n=200,method="BH")<0.05)
print(paste("Number of Regions with Significantly Different TE between MCI and HC: ",
            length(significant.regions)))

#Find regions where TE is greater in MCI and display them in order of magnitude
mci.greater=significant.regions[which(lm.estimates[significant.regions]>0)]
print(paste("Number of Regions where TE is Significantly Greater in MCI than in HC:",
            length(mci.greater)))
print(paste("Regions (decreasing strength): ",toString(mci.greater[order(lm.estimates[mci.greater],decreasing=TRUE)])))

#Find regions where TE is greater in HC and display them in order of magnitude
hc.greater=significant.regions[which(lm.estimates[significant.regions]<0)]
print(paste("Number of Regions where TE is Significantly Greater in HC than in MCI:",
            length(hc.greater)))
print(paste("Regions (decreasing strength): ",toString(hc.greater[order(lm.estimates[hc.greater],decreasing=FALSE)])))



# ---- ANALYZE REGIONAL ENTROPIES ---- #

#Import and initialize data for the loop
regional.entropy=load(file.path(results.dir,"Entropy.RData"))
lm.estimates=NULL #variable to store lm estimates for each region
lm.pvalue=NULL #variable to store p-value for each region

#Compare entropy for each region using the linear model adjusted for sex, age, and ICV
for(i in 1:200) { #adjust bounds to match number of regions (and columns storing the data)
  lm.entry=summary(lm(regional.entropy[,i]~Status+Sex+Age+ICV,data=dem.data))
  lm.estimates=rbind(lm.estimates,lm.entry$coefficients[2,1]) #isolate coefficients
  lm.pvalue=rbind(lm.pvalue,lm.entry$coefficients[2,4]) #isolate p-value for status
}

#Create MATLAB files that will be used to make the brain plots for comparisons
lm.estimates.status=lm.estimates[,1]
writeMat(con=file.path(data.dir,"Entropy.Estimates.Status.mat"),labpcexport=lm.estimates.status)

#Create MATLAB files containing the regional entropy for HC and MCI separately
regional.HC=colMeans(regional.entropy[which(dem.data$Status=="HC"),])
writeMat(con=file.path(data.dir,"Regional.Entropy.HC.mat"),labpcexport=regional.HC)
regional.MCI=colMeans(regional.entropy[which(dem.data$Status=="MCI"),])
writeMat(con=file.path(data.dir,"Regional.Entropy.MCI.mat"),labpcexport=regional.MCI)

#Identify regions with significantly different entropy between HC and MCI after BH correction
significant.regions=which(p.adjust(lm.pvalue,n=200,method="BH")<0.05)
print(paste("Number of Regions with Significantly Different Entropy between MCI and HC: ",
            length(significant.regions)))

#Find regions where entropy is greater in MCI and display them in order of magnitude
mci.greater=significant.regions[which(lm.estimates[significant.regions]>0)]
print(paste("Number of Regions where Entropy is Significantly Greater in MCI than in HC:",
            length(mci.greater)))
print(paste("Regions (decreasing strength): ",toString(mci.greater[order(lm.estimates[mci.greater],decreasing=TRUE)])))

#Find regions where entropy is greater in HC and display them in order of magnitude
hc.greater=significant.regions[which(lm.estimates[significant.regions]<0)]
print(paste("Number of Regions where Entropy is Significantly Greater in HC than in MCI:",
            length(hc.greater)))
print(paste("Regions (decreasing strength): ",toString(hc.greater[order(lm.estimates[hc.greater],decreasing=FALSE)])))



# ---- ANALYZE REGIONAL TAU LEVELS ---- #

#Import and initialize data for the loop
regional.tau=load(file.path(data.dir,"Tau.Clean.RData"))
lm.estimates=NULL #variable to store lm estimates for each region
lm.pvalue=NULL #variable to store p-value for each region

#Compare tau for each region using the linear model adjusted for sex, age, and ICV
for(i in 2:201) { #adjust bounds to match number of regions (and columns storing the data)
  lm.entry=summary(lm(regional.tau[,i]~Status+Sex+Age+ICV,data=dem.data))
  lm.estimates=rbind(lm.estimates,lm.entry$coefficients[2,1]) #isolate coefficients
  lm.pvalue=rbind(lm.pvalue,lm.entry$coefficients[2,4]) #isolate p-value for status
}

#Create MATLAB files that will be used to make the brain plots for comparisons
lm.estimates.status=lm.estimates[,1]
writeMat(con=file.path(data.dir,"Tau.Estimates.Status.mat"),labpcexport=lm.estimates.status)

#Create MATLAB files containing the regional tau for HC and MCI separately
regional.HC=colMeans(regional.tau[which(dem.data$Status=="HC"),])
writeMat(con=file.path(data.dir,"Regional.Tau.HC.mat"),labpcexport=regional.HC)
regional.MCI=colMeans(regional.tau[which(dem.data$Status=="MCI"),])
writeMat(con=file.path(data.dir,"Regional.Tau.MCI.mat"),labpcexport=regional.MCI)

#Identify regions with significantly different tau between HC and MCI after BH correction
significant.regions=which(p.adjust(lm.pvalue,n=200,method="BH")<0.05)
print(paste("Number of Regions with Significantly Different Tau between MCI and HC: ",
            length(significant.regions)))

#Find regions where tau is greater in MCI and display them in order of magnitude
mci.greater=significant.regions[which(lm.estimates[significant.regions]>0)]
print(paste("Number of Regions where Tau is Significantly Greater in MCI than in HC:",
            length(mci.greater)))
print(paste("Regions (decreasing strength): ",toString(mci.greater[order(lm.estimates[mci.greater],decreasing=TRUE)])))

#Find regions where tau is greater in HC and display them in order of magnitude
hc.greater=significant.regions[which(lm.estimates[significant.regions]<0)]
print(paste("Number of Regions where Tau is Significantly Greater in HC than in MCI:",
            length(hc.greater)))
print(paste("Regions (decreasing strength): ",toString(hc.greater[order(lm.estimates[hc.greater],decreasing=FALSE)])))



# ---- ANALYZE REGIONAL AMYLOID BETA PLAQUE LEVELS ---- #

#Import and initialize data for the loop
regional.plaque=load(file.path(data.dir,"Plaque.Clean.RData"))
lm.estimates=NULL #variable to store lm estimates for each region
lm.pvalue=NULL #variable to store p-value for each region

#Compare plaque for each region using the linear model adjusted for sex, age, and ICV
for(i in 2:201) { #adjust bounds to match number of regions (and columns storing the data)
  lm.entry=summary(lm(regional.plaque[,i]~Status+Sex+Age+ICV,data=dem.data))
  lm.estimates=rbind(lm.estimates,lm.entry$coefficients[2,1]) #isolate coefficients
  lm.pvalue=rbind(lm.pvalue,lm.entry$coefficients[2,4]) #isolate p-value for status
}

#Create MATLAB files that will be used to make the brain plots for comparisons
lm.estimates.status=lm.estimates[,1]
writeMat(con=file.path(data.dir,"Plaque.Estimates.Status.mat"),labpcexport=lm.estimates.status)

#Create MATLAB files containing the regional plaque for HC and MCI separately
regional.HC=colMeans(regional.plaque[which(dem.data$Status=="HC"),])
writeMat(con=file.path(data.dir,"Regional.Plaque.HC.mat"),labpcexport=regional.HC)
regional.MCI=colMeans(regional.plaque[which(dem.data$Status=="MCI"),])
writeMat(con=file.path(data.dir,"Regional.Plaque.MCI.mat"),labpcexport=regional.MCI)

#Identify regions with significantly different plaque between HC and MCI after BH correction
significant.regions=which(p.adjust(lm.pvalue,n=200,method="BH")<0.05)
print(paste("Number of Regions with Significantly Different Plaque between MCI and HC: ",
            length(significant.regions)))

#Find regions where plaque is greater in MCI and display them in order of magnitude
mci.greater=significant.regions[which(lm.estimates[significant.regions]>0)]
print(paste("Number of Regions where Plaque is Significantly Greater in MCI than in HC:",
            length(mci.greater)))
print(paste("Regions (decreasing strength): ",toString(mci.greater[order(lm.estimates[mci.greater],decreasing=TRUE)])))

#Find regions where plaque is greater in HC and display them in order of magnitude
hc.greater=significant.regions[which(lm.estimates[significant.regions]<0)]
print(paste("Number of Regions where Plaque is Significantly Greater in HC than in MCI:",
            length(hc.greater)))
print(paste("Regions (decreasing strength): ",toString(hc.greater[order(lm.estimates[hc.greater],decreasing=FALSE)])))



# ---- EXECUTE BLOCKS 12-15 IN MATLAB ---- #



# ---- TEST REGIONAL ASSOCIATION BETWEEN TE AND ENTROPY ---- #

#Initialize variables for TE vs. entropy association
lm.entry.TE.entropy=NULL;lm.results.TE.entropy=NULL
cor.entry.TE.entropy=NULL;cor.results.TE.entropy=NULL

#Use lm and Spearman correlation to examine association between TE and entropy
for(i in 1:200) { #bounds of i should match the number of region
  lm.entry.TE.entropy=summary(lm(regional.TE[,i]~entropy[,i]+
                                   dem.data$Sex+dem.data$Age+dem.data$ICV))
  lm.results.TE.entropy=rbind(lm.results.TE.entropy,lm.entry.TE.entropy$coefficients[2,1:4]) #isolate coefficients
  cor.entry.TE.entropy=cor.test(as.numeric(regional.TE[,i]),as.numeric(entropy[,i]),method="spearman")
  cor.results.TE.entropy=rbind(cor.results.TE.entropy,data.frame(cor.entry.TE.entropy$estimate,cor.entry.TE.entropy$p.value))
}

#Save lm and correlation estimates as MATLAB files
cor.estimates.TE.entropy=cor.results.TE.entropy[,1]
writeMat(con=file.path(results.dir,"TE.Entropy.Cor.Estimates.mat"),labpcexport=cor.estimates.TE.entropy)
lm.estimates.TE.entropy=lm.results.TE.entropy[,1]
writeMat(con=file.path(results.dir,"TE.Entropy.lm.Estimates.mat"),labpcexport=lm.estimates.TE.entropy)

#Display results of the lm association, including BH-adjusted significance
print(paste("Number of Regions with Significant Association between TE and Entropy:",length(which(p.adjust(lm.results.TE.entropy[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with any Positive Association between TE and Entropy:",length(which(p.adjust(lm.results.TE.entropy[,4],n=200,method="BH")<0.05))>0))
print(paste("Number of Regions with any Negative Association between TE and Entropy:",length(which(p.adjust(lm.results.TE.entropy[,4],n=200,method="BH")<0.05))<0))
print(paste("Regions with Significant Association between TE and Entropy:",toString(which(p.adjust(lm.results.TE.entropy[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with Significant Correlation between TE and Entropy:",length(which(p.adjust(cor.results.TE.entropy[,2],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Correlation between TE and Entropy:",toString(which(p.adjust(cor.results.TE.entropy[,2],n=200,method="BH")<0.05))))



# ---- EXECUTE BLOCKS 16-17 IN MATLAB ---- #



# ---- TEST REGIONAL CORRELATIONS AND ASSOCIATIONS BETWEEN PET METRICS, TE, AND ENTROPY ---- #

#Isolate data for MCI only
regional.TE.MCI=regional.TE[which(dem.data$Status=="MCI"),]
entropy.MCI=entropy[which(dem.data$Status=="MCI"),]
plaque.MCI=plaque[which(dem.data$Status=="MCI"),]
tau.MCI=tau[which(dem.data$Status=="MCI"),]

#Initialize variables for TE vs. plaque association
lm.entry.TE.plaque=NULL;lm.results.TE.plaque=NULL
cor.entry.TE.plaque=NULL;cor.results.TE.plaque=NULL

#Initialize variables for TE vs. tau association
lm.entry.TE.tau=NULL;lm.results.TE.tau=NULL
cor.entry.TE.tau=NULL;cor.results.TE.tau=NULL

#Initialize variables for entropy vs. plaque association
lm.entry.entropy.plaque=NULL;lm.results.entropy.plaque=NULL
cor.entry.entropy.plaque=NULL;cor.results.entropy.plaque=NULL

#Initialize variables for entropy vs. tau association
lm.entry.entropy.tau=NULL;lm.results.entropy.tau=NULL
cor.entry.entropy.tau=NULL;cor.results.entropy.tau=NULL

#Use lm and Spearman correlation to examine association for MCI only
for(i in 1:200) { #bounds of i should match the number of region
  lm.entry.TE.plaque=summary(lm(plaque.MCI[,i+1]~regional.TE.MCI[,i]+
                                  dem.data[which(dem.data$Status=="MCI"),]$Sex+
                                  dem.data[which(dem.data$Status=="MCI"),]$Age+
                                  dem.data[which(dem.data$Status=="MCI"),]$ICV))
  lm.results.TE.plaque=rbind(lm.results.TE.plaque,lm.entry.TE.plaque$coefficients[2,1:4]) #isolate coefficients
  cor.entry.TE.plaque=cor.test(as.numeric(plaque.MCI[,i+1]),as.numeric(regional.TE.MCI[,i]),method="spearman")
  cor.results.TE.plaque=rbind(cor.results.TE.plaque,data.frame(cor.entry.TE.plaque$estimate,cor.entry.TE.plaque$p.value))
  
  
  lm.entry.TE.tau=summary(lm(tau.MCI[,i+1]~regional.TE.MCI[,i]+
                               dem.data[which(dem.data$Status=="MCI"),]$Sex+
                               dem.data[which(dem.data$Status=="MCI"),]$Age+
                               dem.data[which(dem.data$Status=="MCI"),]$ICV))
  lm.results.TE.tau=rbind(lm.results.TE.tau,lm.entry.TE.tau$coefficients[2,1:4])
  cor.entry.TE.tau=cor.test(as.numeric(tau.MCI[,i+1]),as.numeric(regional.TE.MCI[,i]),method="spearman")
  cor.results.TE.tau=rbind(cor.results.TE.tau,data.frame(cor.entry.TE.tau$estimate,cor.entry.TE.tau$p.value))
  
  
  lm.entry.entropy.plaque=summary(lm(plaque.MCI[,i+1]~entropy.MCI[,i]+
                                       dem.data[which(dem.data$Status=="MCI"),]$Sex+
                                       dem.data[which(dem.data$Status=="MCI"),]$Age+
                                       dem.data[which(dem.data$Status=="MCI"),]$ICV))
  lm.results.entropy.plaque=rbind(lm.results.entropy.plaque,lm.entry.entropy.plaque$coefficients[2,1:4])
  cor.entry.entropy.plaque=cor.test(as.numeric(plaque.MCI[,i+1]),as.numeric(entropy.MCI[,i]),method="spearman")
  cor.results.entropy.plaque=rbind(cor.results.entropy.plaque,data.frame(cor.entry.entropy.plaque$estimate,cor.entry.entropy.plaque$p.value))
  
  
  lm.entry.entropy.tau=summary(lm(tau.MCI[,i+1]~entropy.MCI[,i]+
                                    dem.data[which(dem.data$Status=="MCI"),]$Sex+
                                    dem.data[which(dem.data$Status=="MCI"),]$Age+
                                    dem.data[which(dem.data$Status=="MCI"),]$ICV))
  lm.results.entropy.tau=rbind(lm.results.entropy.tau,lm.entry.entropy.tau$coefficients[2,1:4])
  cor.entry.entropy.tau=cor.test(as.numeric(tau.MCI[,i+1]),as.numeric(entropy.MCI[,i]),method="spearman")
  cor.results.entropy.tau=rbind(cor.results.entropy.tau,data.frame(cor.entry.entropy.tau$estimate,cor.entry.entropy.tau$p.value))
  
}

#Save lm and correlation estimates as MATLAB files
cor.estimates.TE.plaque=cor.results.TE.plaque[,1]
writeMat(con=file.path(results.dir,"TE.Plaque.Cor.Estimates.mat"),labpcexport=cor.estimates.TE.plaque)
lm.estimates.TE.plaque=lm.results.TE.plaque[,1]
writeMat(con=file.path(results.dir,"TE.Plaque.lm.Estimates.mat"),labpcexport=lm.estimates.TE.plaque)

cor.estimates.TE.tau=cor.results.TE.tau[,1]
writeMat(con=file.path(results.dir,"TE.Tau.Cor.Estimates.mat"),labpcexport=cor.estimates.TE.tau)
lm.estimates.TE.tau=lm.results.TE.tau[,1]
writeMat(con=file.path(results.dir,"TE.Tau.lm.Estimates.mat"),labpcexport=lm.estimates.TE.tau)

cor.estimates.entropy.plaque=cor.results.entropy.plaque[,1]
writeMat(con=file.path(results.dir,"Entropy.Plaque.Cor.Estimates.mat"),labpcexport=cor.estimates.entropy.plaque)
lm.estimates.entropy.plaque=lm.results.entropy.plaque[,1]
writeMat(con=file.path(results.dir,"Entropy.Plaque.lm.Estimates.mat"),labpcexport=lm.estimates.entropy.plaque)

cor.estimates.entropy.tau=cor.results.entropy.tau[,1]
writeMat(con=file.path(results.dir,"Entropy.Tau.Cor.Estimates.mat"),labpcexport=cor.estimates.entropy.tau)
lm.estimates.entropy.tau=lm.results.entropy.tau[,1]
writeMat(con=file.path(results.dir,"Entropy.Tau.lm.Estimates.mat"),labpcexport=lm.estimates.entropy.tau)

#Display results of the lm association, including BH-adjusted significance
print(paste("Number of Regions with Significant Association between TE and Plaque:",length(which(p.adjust(lm.results.TE.plaque[,4],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Association between TE and Plaque:",toString(which(p.adjust(lm.results.TE.plaque[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with Significant Correlation between TE and Plaque:",length(which(p.adjust(cor.results.TE.plaque[,2],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Correlation between TE and Plaque:",toString(which(p.adjust(cor.results.TE.plaque[,2],n=200,method="BH")<0.05))))

print(paste("Number of Regions with Significant Association between TE and Tau:",length(which(p.adjust(lm.results.TE.tau[,4],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Association between TE and Tau:",toString(which(p.adjust(lm.results.TE.tau[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with Significant Correlation between TE and Tau:",length(which(p.adjust(cor.results.TE.tau[,2],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Correlation between TE and Tau:",toString(which(p.adjust(cor.results.TE.tau[,2],n=200,method="BH")<0.05))))

print(paste("Number of Regions with Significant Association between Entropy and Plaque:",length(which(p.adjust(lm.results.entropy.plaque[,4],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Association between Entropy and Plaque:",toString(which(p.adjust(lm.results.entropy.plaque[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with Significant Correlation between Entropy and Plaque:",length(which(p.adjust(cor.results.entropy.plaque[,2],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Correlation between Entropy and Plaque:",toString(which(p.adjust(cor.results.entropy.plaque[,2],n=200,method="BH")<0.05))))

print(paste("Number of Regions with Significant Association between Entropy and Tau:",length(which(p.adjust(lm.results.entropy.tau[,4],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Association between Entropy and Tau:",toString(which(p.adjust(lm.results.entropy.tau[,4],n=200,method="BH")<0.05))))
print(paste("Number of Regions with Significant Correlation between Entropy and Tau:",length(which(p.adjust(cor.results.entropy.tau[,2],n=200,method="BH")<0.05))))
print(paste("Regions with Significant Correlation between Entropy and Tau:",toString(which(p.adjust(cor.results.entropy.tau[,2],n=200,method="BH")<0.05))))



# ---- EXECUTE LINES 649-732 IN MATLAB ---- #



# ---- CREATE BAR PLOTS TO SHOW TE IN EACH NETWORK ---- #

#Initialize data and variables
regional.TE.HC=regional.TE[which(dem.data$Status=="HC"),] #isolate regional TE for HC
regional.TE.MCI=regional.TE[which(dem.data$Status=="MCI"),] #isolate regional TE for MCI
networknames=c("VIS","SOM","DAN","VAN","LIM","FPN","DMN") #list the Yeo 7 network names
lm.p.values=NULL;lm.estimates=NULL #variable to store lm estimates and p-values
TE.individual.network=NULL;TE.individual.network.HC=NULL;TE.individual.network.MCI=NULL #variable to store TE for network i
TE.network.HC=NULL;TE.network.MCI=NULL #variable to store TE for all networks


for (i in 1:7){ #bounds of i match how each network is indexed in the sch.200.yeo file
  TE.individual.network.HC=rowMeans(regional.TE.HC[,which(sch.200.yeo==i)]) #find the average TE in HC for all regions in network i
  TE.network.HC=rbind(TE.network.HC,mean(TE.individual.network.HC)) #add network i to the variable storing all results for HC
  TE.individual.network.MCI=rowMeans(regional.TE.MCI[,which(sch.200.yeo==i)]) #find the average TE in MCI for all regions in network i
  TE.network.MCI=rbind(TE.network.MCI,mean(TE.individual.network.MCI)) #add network i to the variable storing all results for MCI
  TE.individual.network=rowMeans(regional.TE[,which(sch.200.yeo==i)]) #find the average TE for all subjects in regions of network i

  lm.per.network=summary(lm(TE.individual.network~Status+Sex+Age+ICV,data=dem.data)) #apply lm to compare TE between HC and MCI
  lm.estimates=rbind(lm.estimates,lm.per.network$coefficients[2,1]) #isolate lm estimates
  lm.p.values=rbind(lm.p.values,lm.per.network$coefficients[2,4]) #isolate lm p-values
}

#Adjust the lm p-values using BH and find which networks are significant
adjusted.lm.p.values=p.adjust(lm.p.values,n=7,method="BH")
networknames[which(adjusted.lm.p.values<0.05)]

#Create bar plot for lm estimates
data.plot=data.frame(Network=c("VIS","SOM","DAN","VAN","LIM","FPN","DMN"),Value=lm.estimates)
png(file=file.path(figures.dir,"TE.Network.Bar.Plot.Comparison.png"),width=600,height=1000)
ggplot(data.plot,aes(x=Network,y=Value,fill=Value)) +
  geom_bar(stat="identity") +
  scale_fill_gradientn(colors=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)) +
  coord_flip() +
  theme_minimal() +
  labs(title="",x="Network",y="Estimates") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  theme(legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),legend.text=element_text(size=30),
        legend.title=element_text(size=30))
dev.off()


#Create bar plot for TE in HC
data.plot=data.frame(Network=c("VIS","SOM","DAN","VAN","LIM","FPN","DMN"),Value=TE.network.HC)
png(file=file.path(figures.dir,"TE.Network.Bar.Plot.HC.png"),width=600,height=1000)

ggplot(data.plot,aes(x=Network,y=Value)) +
  geom_bar(stat="identity",fill=plasma.colors[6]) +
  coord_flip() +
  theme_minimal() +
  labs(title="",x="Network",y="TE") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  theme(legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),legend.text=element_text(size=30),
        legend.title=element_text(size=30))
dev.off()


#Create bar plot for TE in MCI
data.plot=data.frame(Network=c("VIS","SOM","DAN","VAN","LIM","FPN","DMN"),Value=TE.network.HC)
png(file=file.path(figures.dir,"TE.Network.Bar.Plot.MCI.png"),width=600,height=1000)

ggplot(data.plot,aes(x=Network,y=Value)) +
  geom_bar(stat="identity",fill=plasma.colors[6]) +
  coord_flip() +
  theme_minimal() +
  labs(title="",x="Network",y="TE") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  theme(legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),legend.text=element_text(size=30),
        legend.title=element_text(size=30))
dev.off()



# ---- CREATE BAR PLOTS FOR THE ASSOCIATION BETWEEN NETWORK-LEVEL TE AND ENTROPY ---- #

#Initialize values for the loop
TE.individual.network=NULL;entropy.individual.network=NULL #TE and entropy for network i
lm.per.network #variable to store the lm result for network i
lm.estimates=NULL;lm.p.values=NULL #variables to store the lm estimates and p-values for all networks


for (i in 1:7){ #bounds of i match how each network is indexed in the sch.200.yeo file
  TE.individual.network=rowMeans(regional.TE[,which(sch.200.yeo==i)]) #isolate TE for network i
  entropy.individual.network=rowMeans(entropy[,which(sch.200.yeo==i)]) #isolate entropy for network i
  lm.per.network=summary(lm(TE.individual.network~entropy.individual.network+Sex+Age+ICV,data=dem.data))
  lm.estimates=rbind(lm.estimates,lm.per.network$coefficients[2,1]) #isolate lm estimate
  lm.p.values=rbind(lm.p.values,lm.per.network$coefficients[2,4]) #isolate lm p-value
}

#Adjust the lm p-values using BH and find which networks are significant
adjusted.lm.p.values=p.adjust(lm.p.values,n=7,method="BH")
networknames[which(adjusted.lm.p.values<0.05)]

data.plot=data.frame(Network=c("VIS","SOM","DAN","VAN","LIM","FPN","DMN"),Value=lm.estimates)
png(file=file.path(figures.dir,"TE.Entropy.Network.Bar.Plot.png"),width=600,height=1000)

ggplot(data.plot,aes(x=Network,y=Value)) +
  geom_bar(stat="identity",fill=plasma.colors[6]) +
  coord_flip() +
  theme_minimal() +
  labs(title="",x="Network",y="Association") +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  theme(legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),legend.text=element_text(size=30),
        legend.title=element_text(size=30))
dev.off()



# ---- ANALYZE GLOBAL TE, ENTROPY, AND PET DATA ---- #

#Compute global entropy, global tau, and global amyloid beta plaque levels
global.TE=readMat(file.path(results.dir,"Global.Minimum.TE.mat"))$globalTE
global.entropy=rowMeans(regional.entropy)
global.plaque=rowMeans(regional.plaque[,2:201],na.rm=TRUE)
global.tau=rowMeans(regional.tau[,2:201],na.rm=TRUE)

#Apply the linear model to compare each metric between HC and MCI
summary(lm(global.TE~Status+Sex+Age+ICV,data=dem.data))
summary(lm(global.entropy~Status+Sex+Age+ICV,data=dem.data))
summary(lm(global.plaque~Status+Sex+Age+ICV,data=dem.data))
summary(lm(global.tau~Status+Sex+Age+ICV,data=dem.data))



# ---- CREATE VIOLIN PLOTS TO SHOW GLOBAL DIFFERENCES ---- #

#Choose colors for HC and MCI
colors=c("#00AFBB","#08EFFF")

#Create and save the violin plot for global TE
data.plot=data.frame(group=dem.data$Status,y=global.TE)
pdf(file=file.path(figures.dir,"Global.TE.Violin.pdf"))

ggplot(data=data.plot,aes(x=group,y=y,fill=group))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="",y="Global Transition Energy")+
  scale_fill_manual(values=rep(colors,1))+
  scale_x_discrete(labels=c("HC","MCI"))+
  theme(legend.position="none",axis.text.x=element_text(size=6))


#Create and save the violin plot for global entropy
data.plot=data.frame(group=dem.data$Status,y=global.entropy)
pdf(file=file.path(figures.dir,"Global.Entropy.Violin.pdf"))

ggplot(data=data.plot,aes(x=group,y=y,fill=group))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="",y="Global Transition Energy")+
  scale_fill_manual(values=rep(colors,1))+
  scale_x_discrete(labels=c("HC","MCI"))+
  theme(legend.position="none",axis.text.x=element_text(size=6))


#Create and save the violin plot for global tau
data.plot=data.frame(group=dem.data$Status,y=global.tau)
pdf(file=file.path(figures.dir,"Global.Tau.Violin.pdf"))

ggplot(data=data.plot,aes(x=group,y=y,fill=group))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="",y="Global Tau")+
  scale_fill_manual(values=rep(colors,1))+
  scale_x_discrete(labels=c("HC","MCI"))+
  theme(legend.position="none",axis.text.x=element_text(size=6))


#Create and save the violin plot for global plaque
data.plot=data.frame(group=dem.data$Status,y=global.plaque)
pdf(file=file.path(figures.dir,"Global.Plaque.Violin.pdf"))

ggplot(data=data.plot,aes(x=group,y=y,fill=group))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1,fill="white")+
  labs(title="",x="",y="Global Plaque")+
  scale_fill_manual(values=rep(colors,1))+
  scale_x_discrete(labels=c("HC","MCI"))+
  theme(legend.position="none",axis.text.x=element_text(size=6))



# ---- CREATE A SCATTERPLOT TO SHOW THE CORRELATION BETWEEN TE AND TRANSITION PROBABILITY ---- #

#Initialize values and run loop
TE.TP.cor=NULL
for(i in 1:554){ #bounds of i should match number of subjects
  TE.TP.cor=rbind(TE.TP.cor,cor.test(minimumTE[i,],transition.prob.persubject[i,]))
}

#Find the correlation for global TE and transition probability
global.cor=cor.test(rank(colMeans(transition.prob.persubject)),rank(colMeans(minimumTE)))

#Create and save the scatterplot
pdf(file=file.path(figures.dir,"Correlation.TE.TP.pdf",width=8,height=8))
par(mfrow=c(1,1))
plot(rank(colMeans(transition.prob.persubject)),rank(colMeans(minimumTE)),
     xlab="Transition Energy (rank)",ylab="Global Transition Energy (rank)",
     main=paste("Pearson Correlation =",round(global.cor$estimate,7),"\n p = ",round(global.cor$p.value,7)),
     col=c("black","blue","red")[z],pch=19)
abline(lm(rank(colMeans(minimumTE))~rank(colMeans(transition.prob.persubject))),col="red",lwd=2)



# ---- CREATE A SCATTERPLOT TO SHOW THE RELATIONSHIP BETWEEN TE AND ENTROPY ---- #

#Find the correlation for global TE and entropy
global.cor=cor.test(rank(colMeans(transition.prob.persubject)),rank(colMeans(minimumTE)))
correlation.TE.entropy=cor.test(global.TE,rowMeans(entropy),method='spearman');correlation.TE.entropy

#Find the lm association between TE and entropy
global.lm.TE.entropy=summary(lm(as.numeric(global.TE)~global.entropy+Sex+Age+ICV,data=dem.data));global.lm.TE.entropy

#Create a data frame for global metrics
data.plot=data.frame(status=dem.data$Status,
                     age=dem.data$Age,
                     sex=dem.data$Sex,
                     icv=dem.data$ICV,
                     TE=global.TE,
                     entropy=rowMeans(entropy),
                     tau=global.tau,
                     plaque=global.plaque)

#Create and save the scatterplot
png(file=file.path(figures.dir,"Scatter.TE.Entropy.png",width=1000,height=600))

colors=c("#B12B8F","#00AFBB") #first color is for HC; second color is for MCI
ggplot(data=data.plot,aes(x=TE,y=entropy)) +
  geom_point(aes(color=status),size=3) + #point colors based on status
  geom_smooth(method="lm",se=TRUE,aes(fill=status,color=status),size=1.5) + #regression calculated based on status
  geom_smooth(method="lm",se=FALSE,colour="black",linetype="dashed",size=2) +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=20),
        plot.title=element_text(size=20)) +
  labs(title="",x="Global Transition Energy",y="Global Entropy") +
  scale_color_manual(values=c("HC"=colors[1],"MCI"=colors[2])) +
  scale_fill_manual(values=c("HC"=colors[1],"MCI"=colors[2]))

#Show results for HC and MCI separately
lm.TE.entropy.HC=summary(lm(as.numeric(TE)~entropy+sex+age+icv,data=data.plot[data.plot$status=="HC",]));lm.TE.entropy.HC
lm.TE.entropy.MCI=summary(lm(as.numeric(TE)~entropy+sex+age+icv,data=data.plot[data.plot$status=="MCI",]));lm.TE.entropy.MCI



# ---- CREATE A SCATTERPLOT FOR TE AND AMYLOID BETA PLAQUE ---- #

#Create a data frame for MCI and HC
data.mci=data.frame(status="MCI",
                    age=dem.data$Age[which(dem.data$Status=="MCI")],
                    sex=dem.data$Sex[which(dem.data$Status=="MCI")],
                    icv=dem.data$ICV[which(dem.data$Status=="MCI")],
                    TE=global.TE[which(dem.data$Status=="MCI")],
                    entropy=(rowMeans(entropy))[which(dem.data$Status=="MCI")],
                    tau=global.tau[which(dem.data$Status=="MCI")],
                    plaque=global.plaque[which(dem.data$Status=="MCI")])
data.hc=data.frame(status="HC",
                   age=dem.data$Age[which(dem.data$Status=="HC")],
                   sex=dem.data$Sex[which(dem.data$Status=="HC")],
                   icv=dem.data$ICV[which(dem.data$Status=="HC")],
                   TE=global.TE[which(dem.data$Status=="HC")],
                   entropy=(rowMeans(entropy))[which(dem.data$Status=="HC")],
                   tau=global.tau[which(dem.data$Status=="HC")],
                   plaque=global.plaque[which(dem.data$Status=="HC")])

#Find the correlation and lm association between global TE and global plaque
correlation.TE.plaque=cor.test(global.TE,global.plaque,method='spearman');correlation.TE.plaque
global.lm.TE.plaque=summary(lm(as.numeric(global.TE)~global.plaque+Sex+Age+ICV,data=dem.data));global.lm.TE.plaque

lm.global.plaque=lm(plaque~TE+age+sex+icv,data=data.plot)
lm.global.plaque.estimates=summary(lm.global.plaque)$coefficients["TE","Estimate"]
lm.global.plaque.p.value=summary(lm.global.plaque)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate All: ",lm.global.plaque.estimates," | p-value: ",lm.global.plaque.p.value))

lm.global.plaque.mci=lm(plaque~TE+age+sex+icv,data=data.mci)
lm.global.plaque.estimates.mci=summary(lm.global.plaque.mci)$coefficients["TE","Estimate"]
lm.global.plaque.p.value.mci=summary(lm.global.plaque.mci)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate MCI: ",lm.global.plaque.estimates.mci," | p-value: ",lm.global.plaque.p.value.mci))

lm.global.plaque.hc=lm(plaque~TE+age+sex+icv,data=data.hc)
lm.global.plaque.estimates.hc=summary(lm.global.plaque.hc)$coefficients["TE","Estimate"]
lm.global.plaque.p.value.hc=summary(lm.global.plaque.hc)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate HC: ",lm.global.plaque.estimates.hc," | p-value: ",lm.global.plaque.p.value.hc))

#Create and save the scatterplot
png(file=file.path(figures.dir,"Scatter.TE.Plaque.png",width=1000,height=600))
colors=c("#B12B8F","#00AFBB") #first color is for HC; second color is for MCI
ggplot(data=rbind(data.mci,data.hc),
       aes(x=TE,y=plaque)) +
  xlim(min(data$TE),max(data$TE)) + #set x-axis limits
  geom_point(aes(color=status),size=3) + #point colors by status
  geom_smooth(method="lm",se=TRUE,aes(fill=status,color=status),size=1.5) + #regression lines
  theme(axis.title=element_text(size=20), #text size
        axis.text=element_text(size=20),
        plot.title=element_text(size=20)) +
  labs(title="",x="Global Transition Energy",y="Global Amyloid Beta Plaque Level") +
  scale_color_manual(values=c("HC"=colors[1],"MCI"=colors[2])) +
  scale_fill_manual(values=c("HC"=colors[1],"MCI"=colors[2]))



# ---- CREATE A SCATTERPLOT FOR TE AND TAU ---- #

#Find the correlation and lm association between global TE and global tau
correlation.TE.tau=cor.test(global.TE,global.tau,method='spearman');correlation.TE.tau
global.lm.TE.tau=summary(lm(as.numeric(global.TE)~global.tau+Sex+Age+ICV,data=dem.data));global.lm.TE.tau

lm.global.tau=lm(tau~TE+age+sex+icv,data=data.plot)
lm.global.tau.estimates=summary(lm.global.tau)$coefficients["TE","Estimate"]
lm.global.tau.p.value=summary(lm.global.tau)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate All: ",lm.global.tau.estimates," | p-value: ",lm.global.tau.p.value))

lm.global.tau.mci=lm(tau~TE+age+sex+icv,data=data.mci)
lm.global.tau.estimates.mci=summary(lm.global.tau.mci)$coefficients["TE","Estimate"]
lm.global.tau.p.value.mci=summary(lm.global.tau.mci)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate MCI: ",lm.global.tau.estimates.mci," | p-value: ",lm.global.tau.p.value.mci))

lm.global.tau.hc=lm(tau~TE+age+sex+icv,data=data.hc)
lm.global.tau.estimates.hc=summary(lm.global.tau.hc)$coefficients["TE","Estimate"]
lm.global.tau.p.value.hc=summary(lm.global.tau.hc)$coefficients["TE","Pr(>|t|)"]
print(paste0("Estimate HC: ",lm.global.tau.estimates.hc," | p-value: ",lm.global.tau.p.value.hc))

#Create and save the scatterplot
png(file=file.path(figures.dir,"Scatter.TE.Tau.png",width=1000,height=600))
colors=c("#B12B8F","#00AFBB") #first color is for HC; second color is for MCI
ggplot(data=rbind(data.mci,data.hc),
       aes(x=TE,y=tau)) +
  xlim(min(data$TE),max(data$TE)) + #set x-axis limits
  geom_point(aes(color=status),size=3) + #point colors by status
  geom_smooth(method="lm",se=TRUE,aes(fill=status,color=status),size=1.5) + #regression lines
  theme(axis.title=element_text(size=20), #text size
        axis.text=element_text(size=20),
        plot.title=element_text(size=20)) +
  labs(title="",x="Global Transition Energy",y="Global Tau Level") +
  scale_color_manual(values=c("HC"=colors[1],"MCI"=colors[2])) +
  scale_fill_manual(values=c("HC"=colors[1],"MCI"=colors[2]))



# ---- CREATE A SCATTERPLOT FOR ENTROPY AND AMYLOID BETA PLAQUE ---- #

#Find the correlation and lm association between global entropy and global plaque
correlation.entropy.plaque=cor.test(rowMeans(entropy),global.plaque,method='spearman');correlation.entropy.plaque
global.lm.entropy.plaque=summary(lm(as.numeric(rowMeans(entropy))~global.plaque+Sex+Age+ICV,data=dem.data));global.lm.entropy.plaque

lm.global.plaque=lm(plaque~entropy+age+sex+icv,data=data.plot)
lm.global.plaque.estimates=summary(lm.global.plaque)$coefficients["entropy","Estimate"]
lm.global.plaque.p.value=summary(lm.global.plaque)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate All: ",lm.global.plaque.estimates," | p-value: ",lm.global.plaque.p.value))

lm.global.plaque.mci=lm(plaque~entropy+age+sex+icv,data=data.mci)
lm.global.plaque.estimates.mci=summary(lm.global.plaque.mci)$coefficients["entropy","Estimate"]
lm.global.plaque.p.value.mci=summary(lm.global.plaque.mci)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate MCI: ",lm.global.plaque.estimates.mci," | p-value: ",lm.global.plaque.p.value.mci))

lm.global.plaque.hc=lm(plaque~entropy+age+sex+icv,data=data.hc)
lm.global.plaque.estimates.hc=summary(lm.global.plaque.hc)$coefficients["entropy","Estimate"]
lm.global.plaque.p.value.hc=summary(lm.global.plaque.hc)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate HC: ",lm.global.plaque.estimates.hc," | p-value: ",lm.global.plaque.p.value.hc))

#Create and save the scatterplot
png(file=file.path(figures.dir,"Scatter.Entropy.Plaque.png",width=1000,height=600))
colors=c("#B12B8F","#00AFBB") #first color is for HC; second color is for MCI
ggplot(data=rbind(data.mci,data.hc),
       aes(x=entropy,y=plaque)) +
  xlim(min(data$TE),max(data$TE)) + #set x-axis limits
  geom_point(aes(color=status),size=3) + #point colors by status
  geom_smooth(method="lm",se=TRUE,aes(fill=status,color=status),size=1.5) + #regression lines
  theme(axis.title=element_text(size=20), #text size
        axis.text=element_text(size=20),
        plot.title=element_text(size=20)) +
  labs(title="",x="Global Entropy",y="Global Amyloid Beta Plaque Level") +
  scale_color_manual(values=c("HC"=colors[1],"MCI"=colors[2])) +
  scale_fill_manual(values=c("HC"=colors[1],"MCI"=colors[2]))



# ---- CREATE A SCATTERPLOT FOR ENTROPY AND TAU ---- #

#Find the correlation and lm association between global entropy and global tau
correlation.entropy.tau=cor.test(rowMeans(entropy),global.tau,method='spearman');correlation.entropy.tau
global.lm.entropy.tau=summary(lm(as.numeric(rowMeans(entropy))~global.tau+Sex+Age+ICV,data=dem.data));global.lm.entropy.tau

lm.global.tau=lm(tau~entropy+age+sex+icv,data=data.plot)
lm.global.tau.estimates=summary(lm.global.tau)$coefficients["entropy","Estimate"]
lm.global.tau.p.value=summary(lm.global.tau)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate All: ",lm.global.tau.estimates," | p-value: ",lm.global.tau.p.value))

lm.global.tau.mci=lm(tau~entropy+age+sex+icv,data=data.mci)
lm.global.tau.estimates.mci=summary(lm.global.tau.mci)$coefficients["entropy","Estimate"]
lm.global.tau.p.value.mci=summary(lm.global.tau.mci)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate MCI: ",lm.global.tau.estimates.mci," | p-value: ",lm.global.tau.p.value.mci))

lm.global.tau.hc=lm(tau~entropy+age+sex+icv,data=data.hc)
lm.global.tau.estimates.hc=summary(lm.global.tau.hc)$coefficients["entropy","Estimate"]
lm.global.tau.p.value.hc=summary(lm.global.tau.hc)$coefficients["entropy","Pr(>|t|)"]
print(paste0("Estimate HC: ",lm.global.tau.estimates.hc," | p-value: ",lm.global.tau.p.value.hc))

#Create and save the scatterplot
png(file=file.path(figures.dir,"Scatter.Entropy.Tau.png",width=1000,height=600))
colors=c("#B12B8F","#00AFBB") #first color is for HC; second color is for MCI
ggplot(data=rbind(data.mci,data.hc),
       aes(x=entropy,y=tau)) +
  xlim(min(data$TE),max(data$TE)) + #set x-axis limits
  geom_point(aes(color=status),size=3) + #point colors by status
  geom_smooth(method="lm",se=TRUE,aes(fill=status,color=status),size=1.5) + #regression lines
  theme(axis.title=element_text(size=20), #text size
        axis.text=element_text(size=20),
        plot.title=element_text(size=20)) +
  labs(title="",x="Global Entropy",y="Global Tau Level") +
  scale_color_manual(values=c("HC"=colors[1],"MCI"=colors[2])) +
  scale_fill_manual(values=c("HC"=colors[1],"MCI"=colors[2]))

