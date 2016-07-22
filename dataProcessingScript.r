## Processing of CDR3 data from PTCL samples
## July 22, 2016
## Scott Brown

###################
#### Libraries ####
###################

if(!"entropy" %in% installed.packages()){
        install.packages("entropy")
    }
library("entropy")

if(!"ggplot2" %in% installed.packages()){
        install.packages("ggplot2")
    }
library("ggplot2")
    

###############
#### Files ####
###############

WORKINGDIR <- "/path/to/working/directory/"

shallowCDR3File <- "combinedShallowMiTCRresultFile.tsv"
## expected format (based on MiTCR output):
# library   sample      chain       aaSeq       nucSeq                  abundance
# lib_id    sampleID    alpha/beta  CAVPVTF     TGTGCTGTTCCTGTTACCTTT   5

shallowSampleFile <- "combinedShallowMiTCRsamplesFile.tsv"
## expected format:
# library   sample      readLength  seqDepth    fastq
# lib_id    sampleID    125         2869630     /path/to/sequence_file.tar.gz

deepCDR3File <- "combinedDeepMiTCRresultFile.tsv"
deepSampleFile <- "combinedDeepMiTCRsamplesFile.tsv"

abTCR_TPM_shallow_file <- "abTCR_TPM_shallow.tsv"
## expected format:
# sample    ENST##########.#    ENST#########.#     ...
# sampleID  76.11               15.2                ...

gdTCR_TPM_shallow_file <- "gdTCR_TPM_shallow.tsv"
    
    
    
    
###################
#### Functions ####
###################

##    read in and do preliminary data cleaning
readAndCleanCDR3 <- function(dataFile){
    cdr3 <- read.table(dataFile, header=T, sep="\t", stringsAsFactors=F)
    names(cdr3) <- c("lib","sample","chain","aaSeq","nucSeq","abundance")
    
    print(paste("There are", nrow(cdr3), "total CDR3s"))
    cdr3 <- cdr3[!grepl("~",cdr3$aaSeq),]
    cdr3 <- cdr3[!grepl("[*]",cdr3$aaSeq),]
    print(paste("There are", nrow(cdr3), "productive CDR3s"))
    
    ## There are cases where then same CDR3 aaSeq was determined to be both alpha and beta (due to ambiguous alignemnts)
    ## Need to choose the chain which has more read support.
    print("Fixing cases where single CDR3 called as both alpha and beta.")
    cdr3 <- do.call(rbind,by(cdr3,cdr3[,c("sample","aaSeq")], function(x){
        df <- as.data.frame(x)
        df$abundance <- as.numeric(df$abundance)
        
        ans <- df[df$abundance == max(df$abundance),]
        if(nrow(ans) > 1){
            ## assign unknown
            ans <- ans[1,]
            ans$chain <- "Ambiguous"
            return(ans)
        }else{
            return(ans)
        }
    }))
    print(paste("There are", nrow(cdr3), "unique CDR3s"))
    
    return(cdr3)
}


##  read in and do preliminary data cleaning on sample info
readAndCleanSample <- function(sampleFile, cdr3){
    depths <- read.table(sampleFile, sep="\t", header=T, stringsAsFactors=F)
    
    depths$barcode <- apply(depths, 1, function(x){
        fp <- x[["fastq"]]
        fps <- strsplit(fp, "_")
        end <- fps[[1]][length(fps[[1]])]
        bar <- strsplit(end,"[.]")[[1]][1]
        return(bar)
    })
    
    ## get total TCR reads for each sample.
    depths$totTCR <- apply(depths, 1, function(x){
        return(sum(cdr3$abundance[cdr3$sample==x[["sample"]]]))
    })
    
    depths$totTCRa <- apply(depths, 1, function(x){
        return(sum(cdr3$abundance[cdr3$sample==x[["sample"]] & cdr3$chain=="alpha"]))
    })
    
    depths$totTCRb <- apply(depths, 1, function(x){
        return(sum(cdr3$abundance[cdr3$sample==x[["sample"]] & cdr3$chain=="beta"]))
    })
    
    depths$sampleFracTCR <- depths$totTCR/depths$seqDepth
    
    depths$label <- ""
    depths$label[grepl("_AbT",depths$sample)] <- "Aberrant"
    depths$label[grepl("_non-AbT",depths$sample)] <- "Non-Aberrant"
    depths$label[grepl("_Ctrl",depths$sample)] <- "Control"
    depths$label[grepl("PTCL182",depths$sample)] <- "Aberrant"  ## manually set as aberrant by flow after updated IHC analysis
    
    return(depths)
}

##  get fraction of TCRs that each clone is responsible for.
getCloneFraction <- function(cdr3, depths){
    
    cdr3$pcTotTCR <- apply(cdr3, 1, function(x){
        return(as.numeric(x[["abundance"]])/depths$totTCR[depths$sample==x[["sample"]]])
    })
    
    cdr3$pcTotTCRchain <- apply(cdr3, 1, function(x){
        if(is.na(x[["chain"]])){
            return(NA)
        }else if(x[["chain"]] == "alpha"){
            return(as.numeric(x[["abundance"]])/depths$totTCRa[depths$sample==x[["sample"]]])
        }else if(x[["chain"]] == "beta"){
            return(as.numeric(x[["abundance"]])/depths$totTCRb[depths$sample==x[["sample"]]])
        }else{
            ## Ambiguous chain
            return(as.numeric(x[["abundance"]])/depths$totTCR[depths$sample==x[["sample"]]])
        }
    })
    
    cdr3$pcTotReads <- apply(cdr3, 1, function(x){
        return(as.numeric(x[["abundance"]])/depths$seqDepth[depths$sample==x[["sample"]]])
    })
    
    return(cdr3)
}

## classify clones as dominant or not
getDominantClones <- function(cdr3){
    
    #thresh <- max(log10(datToUse$pcTotTCRchain*datToUse$pcTotReads))
    #threshAlpha <- max(log10(datToUse$pcTotTCRchain[datToUse$chain=="alpha"]*datToUse$pcTotReads[datToUse$chain=="alpha"]))
    #threshBeta <- max(log10(datToUse$pcTotTCRchain[datToUse$chain=="beta"]*datToUse$pcTotReads[datToUse$chain=="beta"]))
    
    ## metric = (pcTotReads x pcTotTCRchain)
    ##        = (abundance/seqDepth x abundance/totTCR[chain])
    ##        = (abundance)^2/(seqDepth x totTCR[chain])
        
    
    cdr3$dominantMetric <- (cdr3$pcTotReads*cdr3$pcTotTCRchain)
    
    ## Use control data to measure background rate, and classify CDR3s as dominant or not.
    datToUse <- subset(cdr3, grepl("Ctrl", sample))
    ctrlAlphaValues <- datToUse$dominantMetric[datToUse$chain=="alpha"]
    threshAlpha <- max(ctrlAlphaValues)
    ctrlBetaValues <- datToUse$dominantMetric[datToUse$chain=="beta"]
    threshBeta <- max(ctrlBetaValues)
    
    #print(paste("Threshold being used is ",thresh,sep=""))
    print(paste("Threshold for alpha is ",threshAlpha,sep=""))
    print(paste("Threshold for beta is ",threshBeta,sep=""))
    
    cdr3$dominant <- FALSE
    cdr3$dominant[cdr3$chain=="alpha"] <- cdr3$dominantMetric[cdr3$chain=="alpha"] > threshAlpha
    cdr3$dominant[cdr3$chain=="beta"] <- cdr3$dominantMetric[cdr3$chain=="beta"] > threshBeta
    #cdr3$dominant[cdr3$chain=="Ambiguous"] <- FALSE #implicit
    return(cdr3)
}

## calculate shannon entropy for each sample
calcEntropy <- function(cdr3, depths){
    #Calculate shannon entropy for each sample
        
    depths$shannon <- apply(depths, 1, function(x){
        return(entropy::entropy(cdr3$abundance[cdr3$sample==x[["sample"]]]))
    })
    depths$numClones <- apply(depths, 1, function(x){
        return(nrow(subset(cdr3, sample==x[["sample"]])))
    })
    depths$normShannon <- depths$shannon / log2(depths$numClones)
    
    return(depths)
}

## estimate tumor purity using relative abundance of dominant clonotype.
estimatePurity <- function(cdr3, depths){
    
    depths$estPurity <- apply(depths, 1, function(x){
        
        ## get all dominant clones for this sample
        doms <- subset(cdr3, sample==x[["sample"]] & dominant)
        
        ## take abundance of max beta if present, otherwise alpha.
        
        if("beta" %in% doms$chain){
            return(max(doms$pcTotTCRchain[doms$chain=="beta"]))
        }else{
            return(max(doms$pcTotTCRchain))
        }
        
    })
    
    depths$estPurity[is.na(depths$estPurity) | depths$estPurity==-Inf] <- 0
    
    depths$hasDominant <- apply(depths, 1, function(x){
        return(TRUE %in% cdr3$dominant[cdr3$sample==x[["sample"]]])
    })
    
    return(depths)
}


##############
#### MAIN ####
##############

setwd(WORKINGDIR)

#### Read in Data ####

shallowCDR3 <- readAndCleanCDR3(shallowCDR3File)
deepCDR3 <- readAndCleanCDR3(deepCDR3File)

shallowSample <- readAndCleanSample(shallowSampleFile, shallowCDR3)
deepSample <- readAndCleanSample(deepSampleFile, deepCDR3)

shallowCDR3 <- getCloneFraction(shallowCDR3, shallowSample)
deepCDR3 <- getCloneFraction(deepCDR3, deepSample)

shallowCDR3 <- getDominantClones(shallowCDR3)
deepCDR3 <- getDominantClones(deepCDR3)

shallowSample <- calcEntropy(shallowCDR3, shallowSample)
deepSample <- calcEntropy(deepCDR3, deepSample)

shallowSample <- estimatePurity(shallowCDR3, shallowSample)
deepSample <- estimatePurity(deepCDR3, deepSample)


## Merge cdr3 and sample to get all dat for plotting and have samples which have no yield.
shallowCDR3 <- merge(shallowCDR3, shallowSample, by="sample", all=T)
deepCDR3 <- merge(deepCDR3, deepSample, by="sample", all=T)


################################
#### Clonal Abundance plots ####
################################

##shallow
shallowCDR3$plotLabel <- paste(shallowCDR3$label,"::",shallowCDR3$sample,"\nEst. Purity: ",round(shallowCDR3$estPurity,4),sep="")
shallowCDR3$chain <- factor(shallowCDR3$chain, levels=c("alpha","beta","Ambiguous"))

(p <- ggplot(shallowCDR3, aes(x=aaSeq, y=pcTotTCRchain, shape=chain, size=log10(pcTotReads), color=dominant)) + geom_point() + facet_wrap(~plotLabel, scales="free_x") + scale_shape_manual(values=c(1,2,0), name="TCR chain") + theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), axis.line = element_line(color = 'black'), legend.key = element_blank(), panel.border = element_blank()) + ylab("Relative Abundance") + xlab("") + scale_color_manual(values=c("#888888","#f7941d"),name="Clonotype",breaks=c(TRUE,FALSE),labels=c("Dominant", "Background")) + scale_size_continuous(name="Read Abundance",breaks=c(-9,-8,-7,-6,-5,-4,-3),labels=c("1 in 1,000,000,000","1 in 100,000,000","1 in 10,000,000","1 in 1,000,000","1 in 100,000","1 in 10,000","1 in 1,000")))
ggsave("shallowAbundancePlot.pdf", p, width=24, height=15)



###############################
#### TCR C gene expression ####
###############################

## shallow

#AlphaBeta

shallowAB <- read.table(abTCR_TPM_shallow_file, header=T, sep="\t", stringsAsFactors=F)

## summing transcripts, then normalizing.
## ENST00000611116.1       ENSG00000277734.5       TRAC
## ENST00000616778.4       ENSG00000277734.5       TRAC

## ENST00000633705.1       ENSG00000211751.8       TRBC1
## ENST00000610416.2       ENSG00000211751.8       TRBC1
## ENST00000632136.1       ENSG00000281981.1       TRBC1
## ENST00000613903.2       ENSG00000281981.1       TRBC1
## ENST00000613594.4       ENSG00000281981.1       TRBC1

## ENST00000622053.4       ENSG00000276849.4       TRBC2
## ENST00000620987.4       ENSG00000276849.4       TRBC2
## ENST00000612153.4       ENSG00000276849.4       TRBC2
## ENST00000613720.4       ENSG00000276849.4       TRBC2
## ENST00000621068.2       ENSG00000276849.4       TRBC2
## ENST00000614992.4       ENSG00000276849.4       TRBC2
## ENST00000466254.1       ENSG00000211772.9       TRBC2

shallowAB$summed.trac <- apply(shallowAB[,c("ENST00000611116.1","ENST00000616778.4")], 1, sum)
shallowAB$summed.trbc1 <- apply(shallowAB[,c("ENST00000633705.1","ENST00000610416.2","ENST00000632136.1","ENST00000613903.2","ENST00000613594.4")], 1, sum)
shallowAB$summed.trbc2 <- apply(shallowAB[,c("ENST00000622053.4","ENST00000620987.4","ENST00000612153.4","ENST00000613720.4","ENST00000621068.2","ENST00000614992.4","ENST00000466254.1")], 1, sum)


shallowABnorm2 <- as.data.frame(cbind(sample=shallowAB$sample, as.data.frame(scale(shallowAB[,c("summed.trac","summed.trbc1","summed.trbc2")]))))
shallowABnorm2$summedNorm.ab <- apply(shallowABnorm2[,c(-1)], 1, sum)



#GammaDelta
shallowGD <- read.table(gdTCR_TPM_shallow_file, header=T, sep="\t", stringsAsFactors=F)

## summing transcripts, then normalizing.
## ENST00000443402.6       ENSG00000211689.7       TRGC1

## ENST00000436911.6       ENSG00000227191.7       TRGC2
## ENST00000610547.1       ENSG00000227191.7       TRGC2

## ENST00000390477.2       ENSG00000211829.7       TRDC

shallowGD$summed.trgc1 <- shallowGD$ENST00000443402.6
shallowGD$summed.trgc2 <- apply(shallowGD[,c("ENST00000436911.6","ENST00000610547.1")], 1, sum)
shallowGD$summed.trdc <- shallowGD$ENST00000390477.2

shallowGDnorm2 <- as.data.frame(cbind(sample=shallowGD$sample, as.data.frame(scale(shallowGD[,c("summed.trgc1","summed.trgc2","summed.trdc")]))))
shallowGDnorm2$summedNorm.gd <- apply(shallowGDnorm2[,c(-1)], 1, sum)


#merge

shallowExp <- merge(shallowABnorm2[,c("sample","summedNorm.ab")], shallowGDnorm2[,c("sample","summedNorm.gd")], by="sample", all=T)

## test the two which have no dominant clones in shallow and deep (PTCL068 and 069)

shallowExp$testOutlier <- 0
shallowExp$testOutlier[shallowExp$sample %in% c("PTCL068_AbT","PTCL029_AbT")] <- 1
shallowExp$testOutlier <- as.factor(shallowExp$testOutlier)

wilcox.test(summedNorm.ab ~ testOutlier, data=shallowExp)
wilcox.test(summedNorm.gd ~ testOutlier, data=shallowExp)

