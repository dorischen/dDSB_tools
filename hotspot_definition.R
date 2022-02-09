### hotspot_definition.R
### Definition of DSB hotspots and their top peak
### Input: dpp files
### Authors: Doris Chen, Franz Klein
### Version: 201221

## PACKAGES
library(data.table)
library(stringr)


## SUB-FUNCTIONS
process_dpp_d <- function(dppTable, chroms, depthCut=0, depthOption="excl", createIndex=FALSE)
{  setkey(dppTable, chrom, position)
   if(dppTable[position<0, .N]>0)
     warning("Negative positions !")
   # reduce to valid chromosomes
   dppTable <- dppTable[chrom %in% chroms,]
   
   # reduce to 3 columns
   if(dim(dppTable)[2]==3) 
   {  vDepthSum <- dppTable[,3,with=FALSE]
   } else 
   {  vDepthSum <- dppTable[,3,with=FALSE] + dppTable[,4,with=FALSE]
   }
   resultTable <- data.table::copy(dppTable[,1:2,with=FALSE])
   resultTable[,depth:=vDepthSum]
   
   if(depthOption=="incl")
   {  resultTable <- resultTable[depth>=depthCut,]
   } else
     resultTable <- resultTable[depth>depthCut,]
   setkey(resultTable, chrom, position)
   
   if(createIndex)
   {  resultTable[,index:=seq(1, nrow(resultTable))]
   }
   
   print(resultTable)
   print(paste(nrow(resultTable),"rows found."))
   
   return(resultTable)
}

# smooth one chromosome (w/o filling - just setting is.na = 0)
ksmooth_per_chrom <- function(inputTable, chromLength, bw, res)  
{  pointnumber <- round(chromLength/res)  
   smooth1 <- ksmooth(inputTable$position, inputTable$depth, "normal", bandwidth=bw, range.x=range(0,(chromLength-1)), n.points=pointnumber)
   smooth1$y[is.na(smooth1$y)] <- 0        
   return(smooth1)
}
 
# smooth one experiment
smooth_dpp_dt <- function(inputTable, bw, res, chroms, chromLengths)
{  dppSmoothed <- data.table()    
   for (chr in chroms)
   {  ks <- ksmooth_per_chrom(inputTable[chrom==chr,], chromLengths[chr], bw, res)
      dppSmoothed <- rbindlist(list(dppSmoothed, data.table(chrom=chr, position=round(ks$x), depth=ks$y)))
   }
   print(paste0("Max. Y:",max(dppSmoothed$depth)))
   return(dppSmoothed)
} 

save_table <- function(inputTable, workingDirOUT, outfile)
{  setwd(workingDirOUT)
   write.table(inputTable, file=outfile, append=FALSE, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
   print(paste(outfile,"saved in ",workingDirOUT))
}


## USER INPUT
# input file (dpp)
workingDirIND <- "input_data/dpp_files"
filePatternD <- "_5prime_dpp.txt"
splitString <- "_5prime"
workingDirOUT <- paste0(workingDirIND,"/kHotspots")

genome <- "ASMv1"
depthCut <- 6
maxGap <- 20
minLength <- 40  # Pan et al., hotspots
minDepthSum <- 175
minDepthSumPerNt <- 1

plusSmoothing <- FALSE
bw <- round(201/2.5)
res <- 1


## DEFAULTS
scriptSuffix <- "201221"
fileSuffix <- paste0("kHotspots",scriptSuffix,"-dCut",depthCut,"-maxG",maxGap,"-minSum",minDepthSum,"-pNt",minDepthSumPerNt)
if(plusSmoothing)
  fileSuffix <- paste0(fileSuffix,"-bw",bw,"res",res)

if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

if(genome=="R64")
{  chroms <- seq(1,16)
   chromNames <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "mitochondrion")
   chromLengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066)  # chrom. lengths 1 - 16  
} else
if(genome=="sk1-mvo1")
{  chroms <- seq(1,16)
   chromNames <- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16")
   chromLengths <- c(203893, 794508, 342718, 1490682, 602514, 284456, 1067526, 544538, 435585, 719294, 687260, 1008248, 908607, 812465, 1054033, 921188)
 } else
if(genome=="ASMv1")
{  chroms <- seq(1,16)
   chromNames <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "mitochondrion")
   chromLengths <- c(228861, 829469, 340914, 1486921, 589812, 299318, 1080440, 542723, 449612, 753937, 690901, 1054145, 923535, 791982, 1053869, 946846)
} else
if(genome=="spombe")
{  chroms <- seq(1,3)
   chromNames <- c("Chr108.2007", "Chr208.2007", "Chr308.2007") 
   chromLengths <- c(5579133, 4539804, 2452883)
} else
  warning(paste("Unknown genome",genome))



## FUNCTION
setwd(workingDirIND)
fileListD <- list.files(pattern=filePatternD) 
print(fileListD)

# go through dpp files
for(f in seq_along(fileListD))
{  inputfile <- fileListD[f]
   print(inputfile)
   
   fileNameBase <- unlist(strsplit(inputfile,splitString))[1]
   print(fileNameBase)
       
   # logfile
   setwd(workingDirOUT)      
   logfile <- paste0(fileNameBase,"_",fileSuffix,".log")
   sink(file=logfile, append=FALSE, type="output", split=TRUE)

   setwd(workingDirIND)
   dppTableSource0 <- fread(inputfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   dppTableSource <- process_dpp_d(dppTableSource0, chroms)  # needed for later calculations
   
   if(plusSmoothing)
   {  dppTable0 <- smooth_dpp_dt(dppTableSource, bw, res, chroms, chromLengths)
   } else
     dppTable0 <- data.table::copy(dppTableSource)
   
   print(paste("Depth cutoff:",depthCut))
   
   dppTable <- dppTable0[depth>depthCut,]
   
   resultTable <- data.table()
   for (chr in chroms) 
   {  chromTable <- dppTable[chrom==chr,]
      setkey(chromTable, position)
      
      # set left and right borders of hotspots
      firstPos <- chromTable[1,position]
      lastPos <- chromTable[nrow(chromTable),position]
      chromTable[,`:=`(distToPrev=position - c(firstPos,chromTable[-nrow(chromTable),position]), distToNext=c(chromTable[-1,position],lastPos) - position)]
      chromTable[,`:=`(leftBorder=ifelse(distToPrev>=maxGap, 1, 0), rightBorder=ifelse(distToNext>=maxGap, 1, 0))]
      chromTable[1, leftBorder:=1]
      chromTable[nrow(chromTable), rightBorder:=1]
               
      # set hotspot nr
      indexL <- which(chromTable$leftBorder==1)
      indexR <- which(chromTable$rightBorder==1)
      chromTable[,hotspotNr:=0]
      for(i in seq_along(indexL))
      {  chromTable[indexL[i]:indexR[i], hotspotNr:=i]
      }
      
      # get max. depth per hotspot
      chromTable[,maxDepth:=max(depth), by=hotspotNr]
       
      # get depth sums per hotspot
      hotspotTable0 <- unique(chromTable[hotspotNr>0,.(chrom, posL=min(position), posR=max(position), depth_sum=sum(depth)), by=hotspotNr])
              
      # add peaks per hotspot
      vPeakPos <- chromTable[depth==maxDepth, min(position), by=hotspotNr]$V1
      hotspotTable0[,peak_pos:=vPeakPos]
               
      # select hotspots above depth cutoff and length
      hotspotTable0[,length:=posR-posL+1]
      hotspotTable0[,depthSumPerNt:=depth_sum/length]
      hotspotTable <- hotspotTable0[depth_sum>=minDepthSum & length>=minLength & depthSumPerNt>=minDepthSumPerNt,]  
               
      # summarize
      resultTable <- rbindlist(list(resultTable, hotspotTable))
   }
   
   hotspotCount <- nrow(resultTable)
   print(paste(hotspotCount,"hotspots found in ", fileNameBase))
         
   # re-assign id
   setkey(resultTable, chrom, posL, posR)
   resultTable[,hotspotNr:=seq(1,nrow(resultTable))]
      
   # show stats
   print(summary(resultTable$length))
   print(summary(resultTable$depth_sum))
   
   # save table
   outfile <- paste0(fileNameBase,"_",fileSuffix,".txt")
   save_table(resultTable, workingDirOUT, outfile)
   
   # get DSBs in and outside of hotspots
   dppTableOri <- dppTableSource[,.(chrom, start=position, end=position, depth)]
   setkey(dppTableOri, chrom, start, end)
   totalSum <- sum(dppTableOri$depth)
   totalSites <- nrow(dppTableOri)
   
   overlapTableSource <- foverlaps(dppTableOri, resultTable, type="within")
   
   noOverlapTable <- overlapTableSource[is.na(hotspotNr),]
   print(paste0(sum(noOverlapTable$depth)," (",round(100*sum(noOverlapTable$depth)/totalSum,1),"%) DSBs at ",nrow(noOverlapTable)," (",round(100*nrow(noOverlapTable)/totalSites,1),"%) DSB sites found outside of hotspots."))
   overlapTable <- overlapTableSource[!is.na(hotspotNr),]
   print(paste0(sum(overlapTable$depth)," (",round(100*sum(overlapTable$depth)/totalSum,1),"%) DSBs at ",nrow(overlapTable)," (",round(100*nrow(overlapTable)/totalSites,1),"%) DSB sites found within hotspots."))
 
   sink(file=NULL)
}

 
 


