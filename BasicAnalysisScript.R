
#' AffyWorkflow
#' @description A function for paried experiment affy array analysis. 
#' This development version only allows one pair of groups with multiple samples
#' and further development may extend the scale for multiple combination of pairs.
#'
#' @param dirInput character, the directory containing the CEL files and RDS file.
#' @param dirOutput character,the directory for your outputs and reports.
#' @param sampleGroup list, NULL, specify group type for each sample
#' @param microArray boolean, plot the microArray picture or not
#' @param chipPseudo boolean, plot the chip pseudo picture or not
#' @param func_enrich character, set the functional enrichment method, e.g. "camera", "mroast", "romer"
#' @param functionality character, set the biological process from MSigDB that you want to inspect.
#'
#' @examples
#' affyWorkflow()
#' affyWorkflow(dirInput = "./data/", dirOutput = "./", sampleGroup = c(rep("Plac", 3), rep("E2", 3)), microArray = TRUE, chipPseudo = TRUE, func_enrich = "camera", functionality = "HALLMARK_MYOGENESIS")
#' 
#' @return the visualisation results and reports for Affymetrix Microarray
#' 
#' @references The Affymetrix Microarray Analysis Basic (Skeleton) Workflow script in course material


# Packages
## Array Statistics pkgs
library(limma)
library(affy)
library(annotate)
library(mouse4302.db)  # load chip-specific annotation

# Data frame combination pkg
library(dplyr)

## Visualisation pkgs
library(ggplot2)
library(affyPLM)  # Chip pseudo images
library(scatterplot3d)
library(pheatmap)  # Heat map


# The main function for Affymetrix Microarray workflow
affyWorkflow <- function(dirInput = "./data/", dirOutput = "./", 
                         sampleGroup = NULL,
                         microArray = TRUE, chipPseudo = FALSE,
                         func_enrich = "camera",
                         functionality = "HALLMARK_MYOGENESIS"){
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
  print("-------- Quality Control Plots --------")
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
        png(paste(dirOutput, "all_mp.png", sep = ""), 
            width     = 4,
            height    = 3,
            units     = "in",
            res       = 600,
            pointsize = 4)
        op <- par(mfrow = c(2,3))
        for (i in 1:length(dfPhenoData$sample)){
          image(mydata[, i], main = rownames(dfPhenoData)[i])
        }
        dev.off()
      }
      # Plot chip pseudo images (affyPLM pkg needed)
      if (chipPseudo == TRUE){
        chipP = fitPLM(mydata)
        png(paste(dirOutput, "all_resids.png", sep = ""), 
            width     = 4,
            height    = 3,
            units     = "in",
            res       = 600,
            pointsize = 4)
        op <- par(mfrow = c(2,3))
        for (i in 1:length(dfPhenoData$sample)){
          image(chipP, which = i, type = 'resids', main = rownames(dfPhenoData)[i])
        }
        dev.off()
      }
    } else {
      print("Skip: The number of samples is larger than 6. Skip microarray pictures")
    }
  } else {
    print("Skip: Skip microarray and chip pseudo images as arguments.")
  }
  
  # Plot histogram for log intensity and density f
  print("Process: Print the histogram...")
  png(paste(dirOutput, "all_hist.png", sep = ""), 
      width     = 4,
      height    = 3,
      units     = "in",
      res       = 600,
      pointsize = 4)
  op = par(mfrow = c(1,1))  # reset the arrangement of the pictures
  hist(mydata, xlab = 'Log intensity', ylab = 'Density')
  legend("topright", rownames(dfPhenoData), col = rownames(dfPhenoData))
  dev.off()
  
  # Plot box plot with different colours per sample group
  print("Process: Print the boxplot...")
  if (is.null(sampleGroup)) {
    # Normal box plot
    png(paste(dirOutput, "all_boxplot.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
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
    png(paste(dirOutput, "all_boxplot.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
    boxplot(mydata, col = lsColours, xaxt = "n")
  }
  # Add axis text with angel
  text(seq_along(mydata), par("usr")[3] - 0.5, labels = rownames(dfPhenoData), 
       srt = 20, adj = 1, xpd = TRUE, cex = 0.8)  # angle, position and size
  dev.off()
  
  # The same plot for the non-normalised raw data
  # Note that the mva.pairs call below only plots a few of the  #samples – you may wish to plot them all but this is slow
  print("Process: Print the MA plot...")
  for (i in 1:(length(unique(sampleGroup))+1)){
    png(paste(dirOutput, rownames(dfPhenoData)[i], rownames(dfPhenoData)[i+(length(unique(sampleGroup))+1)], "_maplot.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
    pairSample <- i + length(sampleGroup)/length(unique(sampleGroup))
    mva.pairs(pm(mydata)[, c(i, pairSample)])
    dev.off()
  }
  
  calls <-mas5calls(mydata)
  calls <- exprs(calls)
  absent <- colSums(calls == 'A')
  present <-colSums(calls=='P')
  present_absent <- cbind(rownames(dfPhenoData), absent, present)
  write.table(present_absent, paste(dirOutput, "present_absent.txt", sep = ""))
  print(paste("Output: Present_absent Table - saved in", dirOutput))
  
  
  #---------- Normalise the data using RMA ----------#
  print("-------- Normalisation(RMA) --------")
  print("Process: RMA...")
  eset <- rma(mydata)
  # Summary of RMA results
  eset
  # To obtain a matrix of the expression values, use exprs() 
  values <- exprs(eset)
  
  
  #---------- Plot Normalised Data ----------#
  # Boxplot to observe the results of normalisation
  # Notice differences with the boxplot from the raw data
  print("Process: Print the boxplot after normalisation...")
  if (is.null(sampleGroup)){
    # Normal box plot
    png(paste(dirOutput, "all_boxplot_normalised.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
    boxplot(values, col = 'white', xaxt = "n") 
  } else {
    # Box plot with groups
    png(paste(dirOutput, "all_boxplot_normalised.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
    boxplot(values, col = lsColours, xaxt = "n")  
  }
  # Add axis text with angel
  text(seq_along(mydata), par("usr")[3] - 0.5, labels = rownames(dfPhenoData), 
       srt = 20, adj = 1, xpd = TRUE, cex = 0.8)  # angle, position and size
  dev.off()
  
  # MA plot of each pair of samples
  print("Process: Print the MVA plot for normalised data...")
  for (i in 1:(length(unique(sampleGroup))+1)){
    png(paste(dirOutput, rownames(dfPhenoData)[i], rownames(dfPhenoData)[i+(length(unique(sampleGroup))+1)], "_maplot_normalised.png", sep = ""), 
        width     = 4,
        height    = 3,
        units     = "in",
        res       = 600,
        pointsize = 4)
    pairSample <- i + length(sampleGroup)/length(unique(sampleGroup))
    mva.pairs(values[, c(i, pairSample)])
    dev.off()
  }
  png(paste(dirOutput, "all_maplot_normalised.png", sep = ""), 
      width     = 4,
      height    = 3,
      units     = "in",
      res       = 600,
      pointsize = 4)
  mva.pairs(values)
  dev.off()
  
  
  #---------- Person Correlation Coefficient ----------#
  print("-------- Correlation and Difference --------")
  # Performs hierarchical clustering with average linkage based on Pearson’s Correlation Coefficient
  print("Process: Perform hierarchical clustering...")
  hc <- hclust(as.dist(1-cor(values, method="pearson")), method="average")
  # Visualise the hc results
  print("Process: Plot the hierarchical tree...")
  png(paste(dirOutput, "all_cluster_dendrogram.png", sep = ""), 
      width     = 4,
      height    = 3,
      units     = "in",
      res       = 600,
      pointsize = 4)
  plot(hc)
  dev.off()
  
  
  #---------- PCA ----------#
  print("Process: Perform PCA...")
  pca <- prcomp(t(values), scale=T)
  # Plot the PCA results
  print("Process: Plot PCA 3D distribution...")
  png(paste(dirOutput, "all_pca.png", sep = ""), 
      width     = 4,
      height    = 3,
      units     = "in",
      res       = 600,
      pointsize = 4)
  s3d<-scatterplot3d(pca$x[,1:3], pch=19, 
                     color = rainbow(length(dfPhenoData$sample)))
  s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
  text(s3d.coords$x, s3d.coords$y, labels = colnames(values), 
       pos = 3,offset = 0.5, cex = 0.8)
  dev.off()
  
  
  #---------- Perform fold change ----------#
  print("Process: Perform fold change...")
  #obtaining a matrix of expression values
  exprsvals <- exprs(eset)
  #RMA outputs log2 data while MAS5 outputs linear data
  #To convert from log…
  exprsvals10 <-2^exprsvals
  #check conversion
  exprsvals[1:10,]
  #converted
  exprsvals10[1:10,]
  
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
  if (is.null(sampleGroup)){
    print("Skip: fold table need sample group information.")
  } else {
    controlGroup <- rownames(dfPhenoData)[which(dfPhenoData$sample == unique(sampleGroup)[1])]
    exprGroup <- rownames(dfPhenoData)[which(dfPhenoData$sample == unique(sampleGroup)[2])]
    control.mean <- apply(exprsvals10[,controlGroup],1,mean)
    expr.mean <- apply(exprsvals10[,exprGroup],1,mean)
    #calculate some fold changes
    expr_control <-expr.mean/control.mean
    expr_diff <- abs(expr_control-1)
    #build a summary table to hold all the data
    all.data <- cbind(control.mean, expr.mean, expr_control)
    dfAllData <- as.data.frame(all.data)[order(as.data.frame(all.data)$expr_control, decreasing = TRUE), ]
    # check the column names
    colnames(dfAllData)
    # write the table of means as an output
    write.table(all.data, file=paste(dirOutput, "group_means.txt", sep = ""), 
                quote = F, sep = "\t", col.names=NA)
    write.table(dfAllData, file=paste(dirOutput, "group_means_diff.txt", sep = ""), 
                quote = F, sep = "\t", col.names=NA)
    print(paste("Output: Fold table - saved in", dirOutput))
  }
  
  #---------- Statistical analysis ----------#
  print("-------- Differential Analysis --------")
  # Check original sample order and annotation information
  sampleNames(eset)
  eset@annotation
  #packages in the annotation package
  ls("package:mouse4302.db")
  ## Building annotation for differential gene identification
  # establish annotation for MOE430v2 modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-results-with-gene.html
  #build an annotation table
  ID <- featureNames(eset)
  Symbol <- getSYMBOL(ID, "mouse4302.db")
  Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
  tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
  tmp[tmp=="NA"] <- NA #fix padding with NA characters 
  #assign as feature data of the current Eset
  fData(eset) <- tmp
  
  #Build the design matrix
  design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
  design
  colnames(design) <- c(unique(sampleGroup)[1],unique(sampleGroup)[2])
  #Check it makes sense
  sampleNames(eset)
  #output the design matrix
  design
  
  #This instructs Limma which comparisons to make
  eval(parse(text=paste("contrastmatrix <- makeContrasts(", 
                        unique(sampleGroup)[1], "-", unique(sampleGroup)[2], 
                        ", levels=design)", sep = "")))
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
  myresults <- topTable(fit2,coef=1, adjust="fdr", number=nrow(eset))
  myresults_LFC <- topTable(fit2,coef=1, adjust="fdr", number=nrow(eset), sort.by = "logFC")
  write.table(myresults, paste(dirOutput, "myresults.txt", sep = ""))
  write.table(myresults_LFC, paste(dirOutput, "myresults_LFC.txt", sep = ""))
  print(paste("Output: Statistics Table - saved in", dirOutput))
  
  #make a venn diagram
  clas <- classifyTestsF(fit2)
  png(paste(dirOutput, "all_venn.png", sep = ""), 
      width     = 4,
      height    = 3,
      units     = "in",
      res       = 600,
      pointsize = 4)
  vennDiagram(clas)
  dev.off()
  
  # pHeatmap
  myIndexes <-rownames(myresults[1:40,])
  myIndexes_LFC <-rownames(myresults_LFC[1:40,])
  plotme <-exprs(eset)[myIndexes,]
  plotme_LFC <-exprs(eset)[myIndexes_LFC,]
  png(paste(dirOutput, "heatmap.png", sep = ""), 
      width     = 8,
      height    = 6,
      units     = "in",
      res       = 300)
  pheatmap(plotme,scale="row")
  dev.off()
  png(paste(dirOutput, "heatmap_LFC.png", sep = ""), 
      width     = 8,
      height    = 6,
      units     = "in",
      res       = 300)
  pheatmap(plotme_LFC,scale="row")
  dev.off()
  
  #Use help for other parameters. Note we might decide to use #exactly the same model as our differential gene analysis for #the enrichment analysis- in this case we can extract it from #the fit
  sv <- squeezeVar(fit$sigma^2, df=fit$df.residual)
  write.table(sv, paste(dirOutput, "sv.txt", sep = ""))
  print(paste("Output: Bayes variances Table - saved in", dirOutput))
  
  #---------- Functional Enrichment analysis ----------#
  print("-------- Functional Enrichment analysis --------")
  Mm.H <- readRDS(paste(dirInput, "Mm.h.all.v7.1.entrez.rds", sep = "")) 
  
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
  
  #convert to indexes
  H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)
  #Pick the most suitable enrichment analysis tool to find #enrichment signatures in the data and run this tool So:-
  
  #I just run mroast here as an example- justify the selection of this method!
  
  #if you want to run mroast
  if (func_enrich == "mroast") {
    results <-mroast(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
  } else if (func_enrich == "camera") {
    #if you want to run camera
    results <-camera(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
  } else if (func_enrich == "romer") {
    #if you want to run romer
    results <-romer(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
  } else {
    stop("Variable Error: Please check if the func_enrich variable is properly settled or quoted.")
  }

  #View the results
  results
  write.table(results,paste(dirOutput, func_enrich, "_enrichment.txt", sep = ""), sep="\t")
  
  # index table for trace back
  dfHIndex <- data.frame()
  for (lsName in names(H.indices)) {
    df <- eval(parse(text=paste("H.indices$", lsName, sep = "")))
    df <- cbind(ENTREZID=df, path=rep(lsName, length(df)))
    dfHIndex <- rbind(dfHIndex, df)
  }
  dfFinal <- left_join(fData(eset_t), dfHIndex, by = c("ENTREZID" = "ENTREZID"))
  dfFinal <- left_join(myresults_LFC, dfFinal, by = c("ID" = "PROBEID"))
  write.table(dfFinal,paste(dirOutput, func_enrich, "_index.txt", sep = ""), sep="\t", row.names = FALSE)
  
  # Generate tables for specific functionality
  dfSPFunc <- dfFinal[which(dfFinal$path == functionality),]
  #dfSPFunc <- dfSPFunc[order(dfSPFunc$P.Value), ]
  write.table(dfSPFunc, paste(dirOutput, func_enrich, functionality, "_index.txt", sep = ""), sep="\t", row.names = FALSE)
  
  # Match top ten with table of the specific functionality
  df10 <- left_join(myresults[1:10, ], dfFinal, by = c("ID" = "ID"))
  df10_LFC <- left_join(myresults_LFC[1:10, ], dfFinal, by = c("ID" = "ID"))
  
  write.table(df10_LFC, paste(dirOutput, func_enrich, "_top10.txt", sep = ""), sep="\t", row.names = FALSE)
  
  print(paste("Output: Functional Enrichment Tables - saved in", dirOutput))
  #You can then examine the results in “enrichment.txt”.  It is a text file.  It can be downloaded to view in a spreadsheet such as Excel.
  
}


## Session Information
sessionInfo()


