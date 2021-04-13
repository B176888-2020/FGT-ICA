
#' AffyWorkflow
#' @description This
#'
#' @param dirInput The directory containing the CEL files and RDS file.
#' @param dirOutput The directory for your outputs and reports.
#'
#' @examples
#' 
#' 
#' @return the visualisation results and reports for Affymetrix Microarray
#' 
#' @references The Affymetrix Microarray Analysis Basic (Skeleton) Workflow script in course material


#Problems (features that might be improved)
# TODO: interactivity and feedback information
# TODO: duplication and layout
# TODO: Enhancement: SimpleAffy, PMA calls, targets file use etc and other missing features
# TODO: plot should be saved



# Packages
## Array Statistics pkgs
library(limma)
library(affy)
library(annotate)
library(mouse4302.db)  # load chip-specific annotation

## Visualisation pkgs
library(ggplot2)
library(affyPLM)  # Chip pseudo images
library(scatterplot3d)
library(pheatmap)  # Heat map


# The main function for Affymetrix Microarray workflow
affyWorkflow <- function(dirInput = "./data/CEL_data/", dirOutput = "./", 
                         sampleGroup = NULL,
                         microArray = TRUE, chipPseudo = FALSE){
  #---------- Input info and confirmation ----------#
  # Confirm the input arguments for users
  print("-------- AffyWorkflow(development version) --------")
  print("-------- Basic Input Information --------")
  print(paste("Input dir:", as.character(dirInput)))
  print(paste("Output dir:", as.character(dirOutput)))
  
  # Load and Read CEL files
  if (is.null(dirInput)) {
    print("Attention: Using the template data from CEL folder.")
  } else if (substr(dirInput, nchar(dirInput), nchar(dirInput)) != "/" || 
           substr(dirOutput, nchar(dirOutput), nchar(dirOutput)) != "/") {
    stop("Input Error: Please check if your inputs are empty or miss / in the end of dir.")
  } else {
    print("Process: Loading CEL data...")
    
  }
  mydata <- ReadAffy(celfile.path = dirInput)
  
  # View a summary of the CEL data
  ## Print the summary
  print(mydata)
  ## List the sample's name and feature
  print(paste("Sample Names:", paste(rownames(pData(mydata)), collapse = ",")))
  
  
  
  #---------- Build Quality Control Plots ----------#
  print("--------Quality Control Plots --------")
  # Separate data for further investigation
  dfPhenoData <- mydata@phenoData@data
  if (!is.null(sampleGroup) && length(sampleGroup) == length(mydata@phenoData@data$sample)){
    dfPhenoData$sample <- sampleGroup
    mydata@phenoData@data$sample <- sampleGroup 
  } else {
    print("Warning: the element number of sampleGroup is not the same as CEL data. Continue without sampleGroup...")
    sampleGroup <- NULL
  }
  
  # Plot microarray pictures and chip pseudo images  to show large inconsistencies
  if (microArray == TRUE || chipPseudo == TRUE){
    print("Process: Print the microarray pictures...")
    if (length(dfPhenoData$sample) <= 6){
      # Plot microarray pictures
      if (microArray == TRUE) {
        op = par(mfrow = c(2,3))
        for (i in 1:length(dfPhenoData$sample)){
          image(mydata[, i], main = rownames(dfPhenoData)[i])
        }
      }
      # Plot chip pseudo images (affyPLM pkg needed)
      if (chipPseudo == TRUE){
        chipP = fitPLM(mydata)
        op = par(mfrow = c(2,3))
        for (i in 1:length(dfPhenoData$sample)){
          image(chipP, which = i, type='resids', main = rownames(dfPhenoData)[i])
        }
      }
    } else {
      print("Skip: The number of samples is larger than 6. Skip microarray pictures")
    }
  } else {
    print("Skip: Skip microarray and chip pseudo images as arguments.")
  }
  
  # Plot histogram for log intensity and density f
  print("Process: Print the histogram...")
  op = par(mfrow = c(1,1))  # reset the arrangement of the pictures
  hist(mydata, xlab = 'Log intensity', ylab = 'Density')
  legend("topright", rownames(dfPhenoData), col = rownames(dfPhenoData))
  
  # Plot box plot with different colours per sample group
  print("Process: Print the boxplot...")
  if (is.null(sampleGroup)) {
    # Normal box plot
    boxplot(mydata, col = 'white', xaxt = "n")  
  } else {
    # Select colors for different groups
    colours <- grDevices::colors()[10:9+length(unique(sampleGroup))]
    lsColours <- sampleGroup
    for (i in 1:length(unique(sampleGroup))) {
      group <- unique(sampleGroup)[i]
      lsColours[which(lsColours == group)] <- colours[i]
    }
    # Box plot with groups
    boxplot(mydata, col = lsColours, xaxt = "n")  
  }
  # Add axis text with angel
  text(seq_along(mydata), par("usr")[3] - 0.5, labels = rownames(dfPhenoData), 
       srt = 30, adj = 1, xpd = TRUE, cex = 0.8)  # angle, position and size
  
  
  #---------- Normalise the data using RMA ----------#
  print("--------Normalisation(RMA) --------")
  print("Process: RMA...")
  eset <- rma(mydata)
  # Summary of RMA results
  eset
  # To obtain a matrix of the expression values, use exprs() 
  values <- exprs(eset)
  
  
  #---------- Plot Normalised Data ----------#
  # Boxplot to observe the results of normalisation
  # Notice differences with the boxplot from the raw data
  if (is.null(sampleGroup)){
    # Normal box plot
    boxplot(values, col = 'white', xaxt = "n") 
  } else {
    # Box plot with groups
    boxplot(values, col = lsColours, xaxt = "n")  
  }
  # Add axis text with angel
  text(seq_along(mydata), par("usr")[3] - 0.5, labels = rownames(dfPhenoData), 
       srt = 30, adj = 1, xpd = TRUE, cex = 0.8)  # angle, position and size
  
  
  
  
  
}







## Plot Normalised Data


# Boxplot to observe the results of normalisation
# Notice differences with the boxplot from the raw data
boxplot(values, col=colours,las=2)

# MA plot of the samples 1 and 4
mva.pairs(values[, c(1,4)])
# The same plot for the non-normalised raw data
# Note that the mva.pairs call below only plots a few of the  #samples – you may wish to plot them all but this is slow
mva.pairs(pm(mydata)[, c(1,4)])

## Plot Heatmap

# To facilitate interpretation, let’s replace the columns # # header,currently
# displaying the filename, to show the name of each sample 
# (if you have a targets file)
colnames(values) <- rownames(pData(adf))
# Performs hierarchical clustering with average linkage based on
# Pearson’s Correlation Coefficient
hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
plot(hc)
## Perform PCA

pca <- prcomp(t(values), scale=T)
# Plot the PCA results

s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(1))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(values),pos = 3,offset = 0.5)

## Perform fold filtering

#obtaining a matrix of expression values
exprsvals <- exprs(eset)
#RMA outputs log2 data while MAS5 outputs linear data
#To convert from log…
exprsvals10 <-2^exprsvals
#check conversion
exprsvals[1:10,]
#converted
exprsvals10[1:10,]

#More fold filtering
#check order of sample names
mysamples <- sampleNames(eset)
#display the list
mysamples
#it is useful to obtain a vector of ProbeIDs here
probesets <- probeNames(mydata)
#display the first 10 ProbeSets
probesets[1:10]

#Build final fold table
#Calculate the means
#Note mean of the log is not the same as the log of the mean!!
Plac.mean <- apply(exprsvals10[,c("GSM3946100_01_Plac.CEL", "GSM3946101_02_Plac.CEL","GSM3946102_03_Plac.CEL")],1,mean)
E2.mean <- apply(exprsvals10[,c("GSM3946103_04_E2.CEL", "GSM3946104_05_E2.CEL","GSM3946105_06_E2.CEL")],1,mean)
#calculate some fold changes
E2c_Plac <-E2.mean/Plac.mean
#build a summary table to hold all the data
all.data <- cbind(Plac.mean, E2.mean, E2c_Plac)
#check the column names
colnames(all.data)
#write the table of means as an output
write.table(all.data, file="group_means.txt", quote=F, sep="\t",col.names=NA)

## Beginning statistical analysis

#Check original sample order
sampleNames(eset)
#Rename the samples
#sampleNames(eset) <-c("ESC.1","ESC.2","ESC.3","iPS2.2","iPS2.3","NSC.1","NSC.2","iPS2.1","iPS4.1","iPS4.2","iPS.3")
#Check the samples have renamed
sampleNames(eset)


##Building annotation for differential gene identification
#establish annotation for MOE430v2
#which annotation do we need
#modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-results-with-gene.html

eset@annotation


#packages in the annotation package
ls("package:mouse4302.db")

#build an annotation table
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters 
#assign as feature data of the current Eset
fData(eset) <- tmp

## Statistical analysis using Limma

#Build the design matrix
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
design
colnames(design) <- c("Plac","E2")
#Check it makes sense
sampleNames(eset)
#output the design matrix
design

#This instructs Limma which comparisons to make
contrastmatrix <- makeContrasts(Plac-E2, levels=design)
contrastmatrix

#issue these commands to fit the model
#and make the contrasts
fit <- lmFit(eset, design)

fit2 <- contrasts.fit(fit, contrastmatrix)

#this last part essentially moderates the t-statistic using 
#the borrowed variance approach described in class
fit2 <- eBayes(fit2)

#get the results
topTable(fit2,coef=1,adjust="fdr")
myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(eset))
write.table(myresults,"myresults.txt")

#make a venn diagram
clas <- classifyTestsF(fit2)
vennDiagram(clas)


## Carry out Functional Enrichment analysis

Mm.H <- readRDS("./data/Mm.h.all.v7.1.entrez.rds") 

#Check that you have the required objects
ls()

#Show the full contents of the annotation package
ls("package:mouse4302.db")

#Show the annotation keys in this database
keytypes(mouse4302.db) 

sampleNames(eset)
## Process annotation for functional enrichment

#Here we select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(eset), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
#View the top of the table
head(res)
#find the index of each row of the expression set in the #annotation object res
idx <- match(rownames(eset), res$PROBEID)
#Use the index to set the phenotypic data in the ExpressionSet
fData(eset) <- res[idx, ]
head(fData(eset), 10)
#Find all rows that don’t have an EntrezID and remove then
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,]


## Functional Enrichment Analysis

#convert to indexes
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)
#Pick the most suitable enrichment analysis tool to find #enrichment signatures in the data and run this tool So:-

#I just run mroast here as an example- justify the selection of this method!

#if you want to run mroast
results <-mroast(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run camera
#results <-camera(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#if you want to run romer
#results <-romer(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#View the results
results
#Use help for other parameters. Note we might decide to use #exactly the same model as our differential gene analysis for #the enrichment analysis- in this case we can extract it from #the fit
#sv <- squeezeVar(fit$sigma^2,df=fit$df.residual)

write.table(results,"enrichment.txt",sep="\t")
#You can then examine the results in “enrichment.txt”.  It is a text file.  It can be downloaded to view in a spreadsheet such as Excel.

## Session Information

sessionInfo()


