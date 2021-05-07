### Peak - valley predictions, optimized for (periodic) dDSB length distributions

## PACKAGES
library(stringr)
library(data.table)


## SUB-FUNCTIONS
add_alpha <- function(vColors, alpha=255) 
{  if(missing(vColors))
     stop("Please provide a vector of colours.")
 
   vColorsT <- as.character(apply(sapply(vColors, col2rgb), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255, alpha=alpha)))
   return(vColorsT)
}

fill_table <- function(distTable, start, end=0)
{  distTable$Length <- as.numeric(distTable$Length)
   if(end==0)
   {  end <- max(distTable$Length)
   }
    
   resultTable0 <- data.table(Length=seq(start,end))
   resultTable <- merge(resultTable0, distTable, by="Length", all.x=TRUE)
   resultTable[is.na(Count), Count:=0]
   print(resultTable)
   
   return(resultTable)
}

smooth_curve_plus_freq <- function(inputTable, bw, res)  # inputTable: Length - Count
{  valueCount <- nrow(inputTable)
   pointCount <- round(valueCount/res)  
   smoothedTable <- ksmooth(inputTable$Length, inputTable$Count, "normal", bandwidth=bw, range.x=range(1,valueCount), n.points=pointCount)
   smoothedTable$y[is.na(smoothedTable$y)] <- 0        
   
   resultTable0 <- data.table(Length=smoothedTable$x, CountS=smoothedTable$y)
   resultTable <- merge(inputTable, resultTable0, by="Length")
   resultTable[,FreqS:=CountS/sum(CountS)]
   print(resultTable)
   
   return(resultTable)
}

create_dist_plot <- function(plotOption, workingDirOUT, outfile, plotTitle, plotTable, countColumn, posColumn, minLength, maxLength, maxY0, tickBin, tickBinY, selColor, lineWidth, gridLinesY, xLabel, plusPeaks=FALSE, plusLegend, pvTable=NULL, period=0, periodSd=0, period2=02, periodSd2=0)   # plotTable .. Length, Freq
{  if(maxY0==0)
   {  maxY <- max(plotTable[,get(countColumn)])*1.1
   } else
     maxY <- maxY0
   
   gridColor <- add_alpha("gray40", alpha=0.5*255)
   axisLabelSize <- 2.8
   axisSize <- 2.2
    
   if(plotOption=="pdf")
   {  setwd(workingDirOUT)
      if(file.exists(outfile))
      {  file.remove(outfile)
         print(paste(outfile,"removed"))
      }
      pdf(file=outfile, bg="transparent", height=9, width=maxLength*15/375) 
   }

   par(xaxs="i", yaxs="i", mgp=c(4, 1.5, 0), mar=c(6,9,2,2), bg="transparent")  # sets position of axis title, tick labels and ticks
 
   plot(plotTable[,get(posColumn)], plotTable[,get(countColumn)], type="l", col="white", lwd=0.1, xlab="", ylab="", xlim=c(minLength, maxLength), ylim=c(0,maxY), xaxt="n", yaxt="n", main=plotTitle, font.main=1, cex.main=1, bty="n")
   axis(side=1, at=seq(minLength, maxLength, tickBin), lwd=2, labels=TRUE, tcl=-0.8, cex.axis=axisSize) 
   mtext(text="Dyad distance (bp)", side=1, line=4.2, cex=axisLabelSize, font=1)
   axis(side=2, at=seq(0, maxY, tickBinY), lwd=2, labels=TRUE, tcl=-0.8, cex.axis=axisSize, las=1) 
   mtext(text="Relative frequency", side=2, line=7, cex=axisLabelSize, font=1)
   lines(plotTable[,get(posColumn)], plotTable[,get(countColumn)], type="l", col=selColor, lwd=lineWidth)
    
   textX <- 0.6
   legendTextSize <- 2
   if(plusPeaks)
   {  pvLwd <- lineWidth-1
      pvTableSel <- pvTable[Length>=minPeak & Length<=peakThresh]
      peakValues <- pvTableSel[Type=="peak",Length]
      valValues <- pvTableSel[Type=="valley",Length]
      
      # valleys
      if(length(valValues)>0)
      {  lines(x=(valValues), y=pvTableSel[Type=="valley",get(countColumn)], type="h", col=valleyColor, lwd=3)
         points(valValues, pvTableSel[Type=="valley",get(countColumn)], col=valleyColor0, cex=legendTextSize, lwd=pvLwd)
         if(plusLegend)
           text(textX*maxLength, 0.65*maxY, paste0("Valley region: ",round(sum(pvTableSel[Type=="valley", .(sum=(get(countColumn)+prevFreq+nextFreq))]),3)*100, "%"), col=valleyColor0, cex=legendTextSize, offset=0, adj=0)
      }
      
      # peaks
      if(length(peakValues)>0)
      {  if(length(selPeakIndices)>0)
         {  text(peakValues[selPeakIndices], c(pvTableSel[Type=="peak",get(countColumn)]+maxY/25)[selPeakIndices], peakValues[selPeakIndices], col=peakColor0, cex=legendTextSize)
         } else
           text(peakValues, pvTableSel[Type=="peak",get(countColumn)]+maxY/33, peakValues, col=peakColor0, cex=legendTextSize)
      
         lines(x=(peakValues), y=pvTableSel[Type=="peak",get(countColumn)], type="h", col=peakColor, lwd=3)
         points(peakValues, pvTableSel[Type=="peak",get(countColumn)], col=peakColor0, cex=legendTextSize, lwd=pvLwd)
         
         if(plusLegend)
           text(textX*maxLength, 0.7*maxY, paste0("Peak region: ",round(sum(pvTableSel[Type=="peak", .(sum=(get(countColumn)+prevFreq+nextFreq))]),3)*100, "%"), col=peakColor0, cex=legendTextSize, offset=0, adj=0)
         
         if(plusLegend & !is.na(period) & period>0)
         {  text(textX*maxLength, 0.95*maxY, paste0("Mean peak distance (",minPeak,"-",peakThresh,"nt): ",round(period,3)," +/- ",round(periodSd,3)), col=selColor, cex=legendTextSize, offset=0, adj=0)
            text(textX*maxLength, 0.9*maxY, paste0("Mean peak distance (",minPeak,"-",peakThresh2,"nt): ",round(period2,3)," +/- ",round(periodSd2,3)), col=selColor, cex=legendTextSize, offset=0, adj=0)
            text(textX*maxLength, 0.85*maxY, paste0("Minimum peak length: ",min(peakValues),"nt"), col=selColor, cex=legendTextSize, offset=0, adj=0)
            text(textX*maxLength, 0.8*maxY, paste0("Maximum peak length: ",max(peakValues),"nt"), col=selColor, cex=legendTextSize, offset=0, adj=0)
            if(restrictPeriod)
            {  periodRange <- range(pvTableSel[!is.na(peakDist) & peakDist>0, peakDist])
               text(textX*maxLength, 0.75*maxY, paste0("Peak distance range: ",periodRange[1]," - ", periodRange[2]), col=selColor, cex=legendTextSize, offset=0, adj=0)
            }
         }
      }
      
   }
   
   if(plotOption!="")
   {  dev.off()
      print(paste(outfile,"saved in",workingDirOUT))
   }
}

call_peaks <- function(inputTable, countColumn, posColumn, thresh, minDist)
{  # collapse repeated values to allow determination of min/max
   uniqueTable <- inputTable[,head(.SD,1), by=.(get(countColumn), rleid(get(countColumn)))]  
   uniqueTable[,get:=NULL]
   lengthU <- nrow(uniqueTable)
  
   # add neigbouring frequencies; row index
   uniqueTable[,`:=`(prevFreq=c(0,uniqueTable[,get(countColumn)])[1:lengthU], nextFreq=c(uniqueTable[,get(countColumn)],0)[2:(lengthU+1)], Index=rleid)]  # add columns with Values shifted 1 row up or 1 row down
   
   # get local maxima and minima 
   maxTable <- uniqueTable[get(countColumn)>prevFreq & get(countColumn)>nextFreq,]   # maxima
   maxTable[,Type:="peak"]
   minTable <- uniqueTable[get(countColumn)<prevFreq & get(countColumn)<nextFreq,]   # minima
   minTable[,Type:="valley"]
 
   # create peak-valley table
   pvTable0 <- rbindlist(list(maxTable,minTable))
   setkeyv(pvTable0, posColumn)
   pvTable0[,`:=`(absDiff=c(abs(diff(pvTable0[,get(countColumn)])),0), pvIndex=seq(1,nrow(pvTable0)))]  # add height difference between peaks and valleys
   
   # remove peaks and valleys below threshold
   if(thresh>0)
   {  cutIndexP <- pvTable0[absDiff<thresh & Type=="peak", pvIndex]
      pvTable1 <- pvTable0[!pvIndex %in% c(cutIndexP,(cutIndexP+1)[cutIndexP+1<=nrow(pvTable0)]),] # remove peaks AND valleys below difference thresholds
      pvTable1[,`:=`(absDiff=c(abs(diff(pvTable1[,get(countColumn)])),0), pvIndex=seq(1,nrow(pvTable1)))]
      
      cutIndexV <- pvTable1[absDiff<thresh & Type=="valley",pvIndex]
      pvTable2 <- pvTable1[!pvIndex %in% c(cutIndexV,(cutIndexV+1)[cutIndexV+1<=nrow(pvTable1)]),] # remove valleys AND peaks below difference thresholds
   } else
     pvTable2 <- pvTable0
   
   # avoid overlapping peak valleys
   pvTable2[,`:=`(distToNext=c(diff(get(posColumn)),0), pvIndex=seq(1,nrow(pvTable2)))]  # add height difference between peaks and valleys
   
   if(minDist>0)
   {  # peaks
      cutIndexP <- pvTable2[distToNext<minDist & Type=="peak", pvIndex]
      if(length(cutIndexP)>0)
      {  pvTable <- pvTable2[!pvIndex %in% c(cutIndexP,(cutIndexP+1)[cutIndexP+1<=nrow(pvTable2)]),] # remove peaks AND valleys below difference thresholds
         pvTable[,`:=`(distToNext=c(diff(pvTable[,get(posColumn)]),0), pvIndex=seq(1,nrow(pvTable)))]
      } else
        pvTable <- pvTable2
      
      # valleys
      cutIndexV <- pvTable[distToNext<minDist & Type=="valley", pvIndex]
      if(length(cutIndexV)>0)
      {  resultTable <- pvTable[!pvIndex %in% c(cutIndexV,(cutIndexV-1)[cutIndexV-1>0]),] # remove peaks AND valleys below difference thresholds
         resultTable[,`:=`(distToNext=c(diff(resultTable[,get(posColumn)]),0), pvIndex=seq(1,nrow(resultTable)))]
      } else
       resultTable <- pvTable
   } else
     resultTable <- pvTable2
    
   resultTable[Type=="peak",peakDist:=c(0,diff(get(posColumn)))]
   print(resultTable)
   print(paste(nrow(resultTable[Type=="peak",]),"peaks found."))
   
   return(resultTable)
}

call_peaks_smoothed <- function(inputTable, thresh)
{  # collapse repeated values to allow determination of min/max
   uniqueTable <- inputTable[,head(.SD,1), by=.(Freq, rleid(Freq))]  
   lengthU <- nrow(uniqueTable)
  
   # add neigbouring frequencies; row index
   uniqueTable[,`:=`(Up=c(0,uniqueTable$Freq)[1:lengthU], Down=c(uniqueTable$Freq,0)[2:(lengthU+1)], prevFreq=c(0,uniqueTable$Freq)[1:lengthU], nextFreq=c(uniqueTable$Freq,0)[2:(lengthU+1)], Index=seq(1,nrow(uniqueTable)))]  # add columns with Values shifted 1 row up or 1 row down
   
   # get local maxima and minima from smoothed values
   maxTable <- uniqueTable[Freq>prevFreq & Freq>nextFreq,]   # maxima
   maxTable[,Type:="peak"]
   minTable <- uniqueTable[Freq<prevFreq & Freq<nextFreq,]   # minima
   minTable[,Type:="valley"]
 
   # create peak-valley table
   pvTable0 <- rbindlist(list(maxTable,minTable))
   pvTable0 <- pvTable0[order(Length),]
   pvTable0[,`:=`(absDiff=c(abs(diff(pvTable0$Freq)),0), pvIndex=seq(1,nrow(pvTable0)))]
   
   # remove peaks and valleys below threshold
   if(thresh>0)
   {  cutIndexP <- pvTable0[absDiff<thresh & Type=="peak",pvIndex]
      pvTable1 <- pvTable0[!pvIndex %in% c(cutIndexP,(cutIndexP+1)[cutIndexP+1<=nrow(pvTable0)]),]
      pvTable1[,`:=`(absDiff=c(abs(diff(pvTable1$Freq)),0), pvIndex=seq(1,nrow(pvTable1)))]
      
      cutIndexV <- pvTable1[absDiff<thresh & Type=="valley",pvIndex]
      pvTable <- pvTable1[!pvIndex %in% c(cutIndexV,(cutIndexV+1)[cutIndexV+1<=nrow(pvTable1)]),]
   } else
     pvTable <- pvTable0
   
   # avoid overlapping peak valleys
   pvTable[,pvIndex:=seq(1,nrow(pvTable))]  # re-assign index (in case of filtered peaks or valleys)
   pvTable[,pvDiff:=c(100,diff(Length))] # for valleys too close to left peak
   minDist <- 3
   if(nrow(pvTable[(Type=="valley") & (pvDiff<minDist),])>0)
   {  for(v in pvTable[(Type=="valley") & (pvDiff<minDist),pvIndex])
      {  if(v<nrow(uniqueTable))
         {  distToPeak <- pvTable[pvIndex==v,pvDiff]
            pvTable[pvIndex==v,Length:=(pvTable[pvIndex==v,Length]+(minDist-distToPeak))]  # move peak by difference minDist - distToPeak
            selIndex <- uniqueTable[Length==pvTable[pvIndex==v,Length],Index]
            if(length(selIndex)>0)
            {  pvTable[pvIndex==v, c("Length", "Count", "Freq", "CountS", "Freq", "Up", "Down", "prevFreq", "nextFreq", "Index"):=uniqueTable[Index==selIndex,]]
            } else
            {  pvTable[pvIndex==v, c("Length", "Count", "Freq", "CountS", "Freq"):=inputTable[Length==pvTable[pvIndex==v,Length],.(Length,Count,Freq,CountS,Freq)]]
            }
         }
      }
   }
   pvTable[,pvDiff:=c(diff(Length),100)] # for valleys too close to right peak
   if(nrow(pvTable[(Type=="valley") & (pvDiff<minDist),])>0)
   {  for(v in pvTable[(Type=="valley") & (pvDiff<minDist),pvIndex])
      {  if(v>1)
         {  distToVal <- pvTable[pvIndex==v,pvDiff]
            pvTable[pvIndex==v,Length:=(pvTable[pvIndex==v,Length]-(minDist-distToVal))]  # move peak by difference minDist - distToPeak
            selIndex <- uniqueTable[Length==pvTable[pvIndex==v,Length],Index]
            if(length(selIndex)>0)
            {  pvTable[pvIndex==v, c("Length", "Count", "Freq", "CountS", "Freq", "Up", "Down", "prevFreq", "nextFreq", "Index"):=uniqueTable[Index==selIndex,]]
            } else
              pvTable[pvIndex==v, c("Length", "Count", "Freq", "CountS", "Freq"):=inputTable[Length==pvTable[pvIndex==v,Length],.(Length,Count,Freq,CountS,Freq)]]
         }
      }
   }
   pvTable[,pvDiff:=c(100,diff(Length))]
   # remove valley if still too close
   pvTable <- pvTable[pvDiff>=minDist]
   
   pvTable[Type=="peak",peakDist:=c(0,diff(Length))]
   print(pvTable)
   
   return(pvTable)
}

save_table <- function(inputTable, workingDirOUT, outfile)
{  setwd(workingDirOUT)
   write.table(inputTable, file=outfile, append=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
   print(paste(outfile,"saved in ",workingDirOUT))
}


# move files to a (new) folder (adapted from http://stackoverflow.com/questions/10266963/moving-files-between-folders)
move_file <- function(from, to)
{  todir <- dirname(to)

   # check if directory is existing, if not then create
   if (!isTRUE(file.info(todir)$isdir)) 
     dir.create(todir, recursive=TRUE)
   
   file.rename(from=from,  to=to)
}



## USER INPUT
# length distribution tables
workingDirIND <- "input_data/tables"
filePatternD <- "_countsAll.txt"
sampleIdIndex <- 1

# plot options
maxX <- 350  # max. value at x-axis (inclusive); if 0, maximum derived from data 
maxY0 <- 0   # max. value at y-axis (inclusive); if 0, maximum derived from data  

# output options
workingDirOUT <- ""
fileNameBaseOv <- "dDSB_ASMv1" # for overview file

# other parameters
format <- "length"  # "length" .. input column names expected to be "Length" "Count"
minDist <- 3  # minimum distance between peaks or valleys
plusSmoothing <- TRUE
bw <- 3  # smoothing bandwidth
res <- 1   # smoothing resolution
pvThresh <- 1.7e-05   # peak - valley difference threshold
minPeak <- 20 # for peak distance calculation, affects max. periodic length
peakThresh <- 350  # max. peak length; for calculation of peak periodicity
peakThresh2 <- 110  # max. peak length; for calculation of peak periodicity
selPeakIndices <- c(1:9,seq(11,30,2))  # for peak labels
restrictPeriod <- FALSE  # restrict peaks to maximum periodic peak length
periodThresh1 <- 8  # for selection of max. periodic length  (otherwise wrongly called peaks could distort value)
periodThresh2 <- 13  # for selection of max. periodic length (otherwise uncalled peaks could distort value)
plotOption <- "pdf"
saveTable <- TRUE

colorTable <- TRUE  # if FALSE vColors needed
workingDirINC <- "input_data/samples"
inputfileC <- "sample_list.txt"
idColumnC <- "Sample.id"

plusLegend <- FALSE
tickBinY <- 0.004
peakColor0 <- "darkred"
valleyColor0 <- "navy"


## DEFAULT VALUES
scriptSuffix <- "200228"
fileSuffix <- paste0(scriptSuffix,"-maxX",maxX,"-minDist",minDist)
if(plusSmoothing)
  fileSuffix <- paste0(fileSuffix,"-bw",bw,"res",res)
if(pvThresh>0)
  fileSuffix <- paste0(fileSuffix,"-thresh",pvThresh)
fileSuffix <- paste0(fileSuffix,"-",minPeak,"-",peakThresh)
if(restrictPeriod)
  fileSuffix <- paste0(fileSuffix,"-periodicOnly")

if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

peakColor <- add_alpha(peakColor0, 0.3*255)  # make transparent
valleyColor <- add_alpha(valleyColor0, 0.3*255)   # make transparent



## SCRIPT
# get input file name(s)
setwd(workingDirIND)
fileListD <- list.files(pattern=filePatternD) 
print(fileListD)

if(colorTable)
{  setwd(workingDirINC)
   sampleTable <- fread(inputfileC, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   colNamesS <- str_replace_all(colnames(sampleTable), " ", ".")
   colnames(sampleTable) <- colNamesS
   sampleTable 
   
   vIds <- unlist(lapply(str_split(fileListD, "_"), "[[", sampleIdIndex))
   
   vColors <- sampleTable[match(vIds,get(idColumnC)), Color]
   print(vColors)
}
vColorsT <- add_alpha(vColors, alpha=0.8*255)

overviewTable <- data.table()
for(d in seq_along(fileListD))
{  inputfileD <- fileListD[d]
   currSample <- unlist(lapply(str_split(inputfileD, "_"), "[[", sampleIdIndex))
   fileNameBase <- currSample
   selColor <- vColorsT[d]

   # get input table 
   setwd(workingDirIND)
   inputTableD0 <- fread(inputfileD, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   
   inputTableD <- data.table::copy(inputTableD0)
   inputTableD <- fill_table(inputTableD, start=1)
   inputTableD[,Freq:=Count/sum(Count)]  # add frequency (before smoothing!!)
   
   if(plusSmoothing)
     inputTableD <- smooth_curve_plus_freq(inputTableD, bw, res)
  
   # create plot
   plotTable <- data.table::copy(inputTableD)
   
   # call peaks and calculate periodicity
   if(plusSmoothing)
   {  freqColumn <- "FreqS"
   } else
     freqColumn <- "Freq"
   pvTable <- call_peaks(plotTable, freqColumn, "Length", pvThresh, minDist)
   peakDistances <- pvTable[Length>=minPeak & Length<=peakThresh & !is.na(peakDist), peakDist][-1]
   periodicity <- mean(peakDistances)
   periodicitySD <- sd(peakDistances)
   peakDistances2 <- pvTable[Length>=minPeak & Length<=peakThresh2 & !is.na(peakDist), peakDist][-1]
   periodicity2 <- mean(peakDistances2)
   periodicitySD2 <- sd(peakDistances2)
   
   # restrict to periodic peaks only
   if(restrictPeriod)
   {  if(length(pvTable[(peakDist>0 & peakDist<periodThresh1) | peakDist>periodThresh2, which=TRUE])>0)  # if all periodic, no restriction necessary
      {  minNotPIndex <- min(pvTable[(peakDist>0 & peakDist<periodThresh1) | peakDist>periodThresh2, which=TRUE])
         pvTable <- pvTable[1:(minNotPIndex-1),]
         
         # re-calculate periodicity
         peakDistances <- pvTable[!is.na(peakDist), peakDist][-1]
         periodicity <- mean(peakDistances)
         periodicitySD <- sd(peakDistances)
      }
   }
   
   # select values below maxX (optional)
   if(maxX>0)
   {  plotTable <- plotTable[Length<=maxX,]
      pvTable <- pvTable[Length<=peakThresh,]
   } else
     maxX <- max(plotTable$Length)
   
   outfile0 <- paste0(fileNameBase,"_peaks",fileSuffix)
   create_dist_plot(plotOption, workingDirOUT, outfile=paste0(outfile0,".pdf"), plotTitle=outfile0, plotTable, countColumn=freqColumn, posColumn="Length", 0, maxX, maxY0, tickBin=10, tickBinY, selColor, lineWidth=2.8, gridLinesY=5, xLabel="Fragment length (nt)", plusPeaks=TRUE, plusLegend, pvTable, periodicity, periodicitySD, periodicity2, periodicitySD2)   # plotTable .. Length, Freq
  
   # save table 
   if(saveTable)
   {  outfile <- paste0(outfile0,"_table.txt")
      save_table(pvTable, workingDirOUT, outfile)
   }
   
   # for overview table
   pvTableSel <- pvTable[Length>=minPeak & Length<=peakThresh]
   periodRange <- range(pvTableSel[!is.na(peakDist) & peakDist>0, peakDist])
   vOverview <- list(currSample, periodicity, periodicitySD, periodicity2, periodicitySD2, periodRange[1], periodRange[2], min(pvTableSel[Type=="peak",Length]), max(pvTableSel[Type=="peak",Length]), max(pvTableSel[Type=="peak" & get(freqColumn)==max(get(freqColumn)),Length]), as.numeric(pvTableSel[Type=="peak",.(100*sum(get(freqColumn),prevFreq,nextFreq))]), as.numeric(pvTableSel[Type=="valley",.(100*sum(get(freqColumn),prevFreq,nextFreq))]))
   overviewTable <- rbindlist(list(overviewTable, vOverview))
}
names(overviewTable) <- c("Sample", "Peak distance avg", "Peak distance std", "Peak distance 2 avg", "Peak distance 2 std", "Min peak distance", "Max peak distance", "Min peak length", "Max peak length", "Max peak", "Perc peak region", "Perc valley region")
print(overviewTable)
if(saveTable)
{  outfile <- paste0(fileNameBaseOv,"_peaks",fileSuffix,"_overviewTable.txt")
   save_table(overviewTable, workingDirOUT, outfile)
}

# move files
fileListT <- list.files(path=workingDirOUT, pattern="_table.txt")
print(fileListT)
setwd(workingDirOUT)
for(file in fileListT)
{  move_file(from=file, to=paste0(workingDirOUT,"/tables/",file))
}












