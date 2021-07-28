### ddsb_HiC_plots.R
### Hi-C type representation of dDSB fragments
### Input: fd files
### Authors: Franz Klein, Doris Chen
### Version 200726

## PACKAGES
library(stringr)
library(data.table)


## SUB-FUNCTIONS
process_dpp_c <- function(dppTable, option="", cal=FALSE, calFactor=1)
{  setkey(dppTable, chrom, position)
   if(dppTable[position<0, .N]>0)
     warning("Negative positions !")

   if(option=="sum")
   {  # reduce to 3 columns
      if(dim(dppTable)[2]==3)
      {  vDepthSum <- dppTable[,3,with=FALSE]
      } else
      {  vDepthSum <- dppTable[,3,with=FALSE] + dppTable[,4,with=FALSE]
      }
      resultTable <- dppTable[,1:2,with=FALSE]
      resultTable[,depth:=vDepthSum]
      if(cal)
        resultTable[,depth:=depth*calFactor]
   } else
   {  resultTable <- dppTable
      if(cal)
        resultTable[,norm_depth_plus:=norm_depth_plus*calFactor, norm_depth_minus:=norm_depth_minus*calFactor]
   }

   return(resultTable)
}

process_fd_c <- function(fdTable, createIndex=FALSE, withCal=FALSE, calFactor=1)  # addition column "end", optionally addition of "index", optionally multiplying depth with calFactor
{  resultTable <- data.table::copy(fdTable)
   resultTable[,end := start + length - 1L]
   if(createIndex)
   {  resultTable[,index := seq(1, nrow(resultTable))]
      setkey(resultTable, index)
   }
   if(withCal)
     resultTable[,depth:=depth*calFactor]
   return(resultTable)
}

add_alpha <- function(vColors, alpha=255)  # alpha between 0 and 1!!
{  if(missing(vColors))
     stop("Please provide a vector of colours.")
   
   vColorsT <- as.character(apply(sapply(vColors, col2rgb), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue=255, alpha=alpha)))
   return(vColorsT)
}


## USER INPUT
# fd files (chrom - start - length - depth)
workingDirINF <- "input_data/fd_files"
filePatternF <- "_fd.txt"
workingDirOUT <- paste0(workingDirINF,"/HiCplots")

# plot
selRegion <- c("3", "214000", "220000") # chrom, start, end
plotHeight0 <- 0 # inch; if 0, derived from data
plotWidth0 <- 0  # inch; if 0, derived from data
maxY0 <- 1000  # max. Y value (dDSB fragment length) in plot; if '0' it will be derived from the data
symbolSize <- 0.6

# colors
rainbow <- TRUE    # if FALSE, then grey
fixedColCount <- 512  # number of colors in gradient; 0 for automatic count

# additional DSB track (vertical lines)
plusDSB <- TRUE  # FALSE or TRUE
workingDirIND <- "input_data/dpp_files"
inputfileD <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_ASMv1-n100_5prime_dpp.txt"
sampleNameD <- "Mohibullah (2017), wt t4"
sampleNameDF <- "keeneyMoh1" # for file name
optionD <- "sumStrands"  # "sepStrands" (Watson and Crick strands separate) or "sumStrands"
scaleFactorD <- 0.75  # scaling factor for signal heights (depth)
lineWidthD <- 0.75 
selColorD <- "#E67300" 
plotOrder <- "DSBtoBack"  # "DSBtoFront" | "DSBtoBack"


## DEFAULTS
if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

scriptSuffix <- "210726"
fileSuffix <- paste0("_HiC",scriptSuffix)
if(maxY0>0)
  fileSuffix <- paste0(fileSuffix,"-maxY",maxY0)
if(!rainbow)
  fileSuffix <- paste0(fileSuffix,"-grey")
if(plusDSB)
  fileSuffix <- paste0(fileSuffix,"-",sampleNameDF,"DSB.",scaleFactorD,".",lineWidthD,".",plotOrder)

selColorDT <- add_alpha(selColorD, alpha=0.5*255)  

if(plotWidth0==0)
  plotWidth <- 8 + (as.numeric(selRegion[3]) - as.numeric(selRegion[2])) * 7/6000



## FUNCTION
# load data
if(plusDSB)
{  # get reference files
   setwd(workingDirIND)
   dppTableSource <- fread(inputfileD, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   if(optionD=="sumStrands")
   {  dppTable <- process_dpp_c(dppTableSource, option="sum")  # chrom - position - depth
   } else
   {  dppTable <- dppTableSource
   }
   print(dppTable)
}

# fd files
fileListF <- list.files(path=workingDirINF, pattern=filePatternF)
print(fileListF)

oldpar <- par()  # save old plot settings

for(f in seq_along(fileListF))
{  inputfileF <- fileListF[f]
   fileNameBase <- str_replace(inputfileF, "_fd.txt", "")
   
   currChrom <- selRegion[1]
   currStart <- as.numeric(selRegion[2])
   currEnd <- as.numeric(selRegion[3])
   print(paste0("chrom", currChrom, ":", currStart, "-", currEnd))
   
   fileNameBaseP <- paste0(fileNameBase,fileSuffix,"_chr",currChrom,".",currStart,"-",currEnd)
   print(fileNameBaseP)

   # get fd table (chrom - start - length - depth - end)
   setwd(workingDirINF)
   fdTable0 <- fread(inputfileF, header=TRUE, sep="\t", stringsAsFactors=FALSE)
   fdTable <- process_fd_c(fdTable0)
   print(fdTable)
 
   # get fd in region
   fdTableSel <- fdTable[chrom==currChrom & start>=currStart & end<=currEnd,]
   fdTableSel <- fdTableSel[order(fdTableSel$depth),]  # for plotting order (strongest fragments last)

   # get dpp in region
   if(plusDSB)
   {  dppTableSel <- dppTable[chrom==currChrom & position>=currStart & position<=currEnd,]
      dppTableSel[,depth:=depth*scaleFactorD]
   }
   
   # set max. depth
   maxDepth <- ceiling(max(fdTableSel$depth))

   # create plot
   plotTitle <- paste0(fileNameBaseP,"\n(max. depth=",maxDepth,", symbol size=",symbolSize,")")
   print(plotTitle)
   
   # set colors
   if(rainbow)
   {  makeramp <- colorRampPalette(c(rgb(0,0,1,0.2),rgb(0,1,0,0.4), rgb(1,1,0,0.6),rgb(1,0.7,0,0.7),rgb(1,0.3,0,0.8), rgb(1,0.15,0,0.9), rgb(1,0,0,1)), alpha=TRUE, bias=4)
   } else
     makeramp <- colorRampPalette(c(rgb(0,0,0,0.2),rgb(0,0,0,0.4), rgb(0,0,0,0.6),rgb(0,0,0,0.8), rgb(0,0,0,1)), alpha=TRUE, bias=1)
   if(fixedColCount==0) 
   {  fixedColCount <- maxDepth
   }
   colorRange <- makeramp(fixedColCount)
   colCountLog <- ceiling(log2(fixedColCount))
   selColCounts <- 2^seq(0,colCountLog)  # for legend
 
   # plot
   factorX <- 10^floor(log10(max(fdTableSel$start)-min(fdTableSel$end)))
   #minX <- floor(currStart/factorX) * factorX
   #maxX <- ceiling(currEnd/factorX) * factorX
   minX <- currStart
   maxX <- currEnd
   factorY <- 10^floor(log10(max(fdTableSel$length)-min(fdTableSel$length))) * 2
   if(maxY0==0)
   {  maxY <- ceiling(max(fdTableSel$length)/factorY)*factorY
   } else
     maxY <- maxY0
   
   if(plotHeight0==0)
     plotHeight <- max(3, (plotWidth - 2.3) * maxY/(as.numeric(selRegion[3])-as.numeric(selRegion[2])))
   
   setwd(workingDirOUT)
   outfile <- paste0(fileNameBaseP,".pdf")
   if(file.exists(outfile))
   {  file.remove(outfile)
      print(paste(outfile,"removed"))
   }
   pdf(file=outfile, bg="transparent", height=plotHeight, width=plotWidth) 
   par(mar=c(3,6,3,2), cex=1, xaxs="i", yaxs="i", bty="n", mgp=c(0,0.8,0))  
   
   # preparation of plot
   plot(0, 0, xlim=c(minX,maxX), ylim=c(0,maxY), pch=".", main=plotTitle, font.main=1, cex.main=0.8, cex.lab=1.5, xaxt="n", yaxt="n", xlab="", ylab="", asp=0.5)
   mtext(text=paste("Position on chromosome",currChrom,"(kb)"), side=1, cex=1.5, line=1)
   axis(side=1, at=seq(minX, maxX, factorX), labels=seq(minX, maxX, factorX)/factorX, tcl=-0.5, cex.axis=1.1, lwd=2, pos=0) # thick ticks
   #axis(side=1, pos=yCross-0.02*yLength, at=seq(minX, maxX, factorX/10), cex.axis=1, lwd=0) # labels of thick ticks
   axis(side=1, at=seq(minX, maxX, factorX/10), labels=FALSE, lwd=1, pos=0) # thin ticks
   mtext(text="dDSB fragment\nlength (bp)", side=2, line=2.5, cex=1.5)
   axis(side=2, at=seq(0, maxY, factorY), labels=TRUE, cex.axis=1.1, lwd=2, pos=minX) 
   
   # DSBs (optional)
   if(plusDSB & plotOrder=="DSBtoBack")
   {  if(optionD=="sumStrands")
      {  lines(x=dppTableSel$position, y=dppTableSel$depth, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
      } else
      {  lines(x=dppTableSel$position, y=dppTableSel$norm_depth_plus, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
         lines(x=dppTableSel$position, y=-dppTableSel$norm_depth_minus, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
      }
      legend(x="topleft", inset=c(0,0), sampleNameD, cex=1, adj=c(0,0.2), pt.cex=1.5, pch="-", col=selColorD, text.col=selColorD, bty="n")  # legend at top left
   }
  
   # dDSBs
   points((fdTableSel$start+fdTableSel$length/2), fdTableSel$length, pch=18, col=colorRange[fdTableSel$depth], cex=symbolSize) 
   
  # DSBs (optional)
   if(plusDSB & plotOrder=="DSBtoFront")
   {  if(optionD=="sumStrands")
      {  lines(x=dppTableSel$position, y=dppTableSel$depth, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
      } else
      {  lines(x=dppTableSel$position, y=dppTableSel$norm_depth_plus, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
         lines(x=dppTableSel$position, y=-dppTableSel$norm_depth_minus, type="h", lwd=lineWidthD, col=selColorDT, asp=1)
      }
      legend(x="topleft", inset=c(0,0), sampleNameD, cex=1, adj=c(0,0.2), pt.cex=1.5, pch="-", col=selColorD, text.col=selColorD, bty="n")  # legend at top left
   }
   
   # legend
   legend(x="topright", horiz=TRUE, legend=as.character(selColCounts), col=colorRange[selColCounts], pch=15, bty="n", pt.cex=1.5, cex=1)

   # saving of pdf file
   dev.off()
   print(paste(outfile,"saved in",workingDirOUT))
}

par(oldpar)  # reset par settings







