### Simulation of DSBs per cell for "random" DSB distance distribution, with user-specified probabilities
### Authors: Doris Chen, Franz Klein
### Version: 210729


## PACKAGES
library(data.table)
library(stringr)


## SUB-FUNCTIONS
multiply_dpp <- function(dppTable, n, addIndex=FALSE)
{  rowCount <- nrow(dppTable)
   resultTable <- data.table(genome=rep(seq(1,n),each=rowCount), chrom=rep(dppTable$chrom,n), position=rep(dppTable$position,n), depth=rep(dppTable$depth,n))
   
   if(addIndex)
     resultTable[,index:=seq(1,nrow(resultTable))]
   
   print(resultTable)
   return(resultTable)
}

replace_na_in_table_b <- function(inputTable, colNames, replaceString, numeric=FALSE)
{  for(colName in colNames)
   {  if(numeric)
      {  inputTable[is.na(get(colName)), (colName):=as.numeric(replaceString)]
      } else
      {  inputTable[is.na(get(colName)), (colName):=replaceString]
      }
   }
 
   return(inputTable)
}

add_average_and_sd <- function(inputTable, countColumnBase, valueColumn, dsbCount)  # e.g. "Count" in case of "Count.1", "Count.2" etc. columns
{  vColumns <- grep(countColumnBase, names(inputTable), value=TRUE)
   resultTable <- replace_na_in_table_b(inputTable, vColumns, "0", numeric=TRUE)
   resultTable[,`:=`(avg=mean(as.numeric(.SD)), sd=sd(.SD)), by=get(valueColumn), .SDcols=vColumns]
   resultTable[,`:=`(avgFreq=avg/sum(avg), sdFreq=sd/sum(avg))]
    
   return(resultTable)
}

save_table <- function(inputTable, workingDirOUT, outfile, sep="\t", plusHeader=TRUE, plusQuote=FALSE)
{  setwd(workingDirOUT)
   write.table(inputTable, file=outfile, append=FALSE, sep=sep, dec=".", row.names=FALSE, col.names=plusHeader, quote=plusQuote)
   print(paste(outfile,"saved in ",workingDirOUT))
}

move_file <- function(from, to)
{  todir <- dirname(to)
   
   # check if directory is existing, if not then create
   if (!isTRUE(file.info(todir)$isdir)) 
     dir.create(todir, recursive=TRUE)

   file.rename(from=from,  to=to)
}


## USER INPUT
# reference DSB map (5prime dpp file)
#workingDirDSB <- "input_data/dpp_files"
workingDirDSB <- "Z:/scripts/GitHub/dDSB_tools/input_data/dpp_files"
inputfileDSB <- "keeneyMoh1_SRR3942949_Spo11-PrA_wt_t4_ASMv1-n100_5prime_dpp.txt"

# for output files
workingDirOUT <- paste0(workingDirDSB,"/ddsb_sim")
fileNameBase0 <- "keeneyMoh1_ASMv1"
saveTables <- TRUE

# further parameters
genome <- "ASMv1"
dsbCount <- 300  # per cell (per meiotic cell 4 haploid genomes)
addCountFactor <- 0.0065 # in case of minDist>0; for correction in simulation
simCount <- 1000  # number of iterations
maxDist <- 70 # incl.
minDist <- 30  # incl.
percdDSB <- 20  # estimated assumed dDSB percentage


## DEFAULTS
scriptSuffix <- "210729"
fileSuffix <- paste0("_dDSBsim",scriptSuffix,"-dsb",dsbCount,"-rep",simCount,"-dist",minDist,"-",maxDist)
fileNameBase <- paste0(fileNameBase0, fileSuffix)

if(!file.exists(workingDirOUT))
{  dir.create(workingDirOUT, recursive=TRUE)
   print(paste(workingDirOUT,"created"))
}

if(genome=="ASMv1")
{  chroms <- seq(1,16)
   chromLengths <- c(228861, 829469, 340914, 1486921, 589812, 299318, 1080440, 542723, 449612, 753937, 690901, 1054145, 923535, 791982, 1053869, 946846)
} 


## SCRIPT
setwd(workingDirOUT)
logfile <- paste0(fileNameBase,".log")
sink(file=logfile, append=FALSE, type="output", split=TRUE)

# load reference DSB table
setwd(workingDirDSB)
dsbTable0 <- fread(inputfileDSB, header=TRUE, sep="\t", stringsAsFactors=FALSE)
dsbTable <- dsbTable0[chrom %in% chroms, .(chrom, position, depth=norm_depth_plus+norm_depth_minus)]
dsbTable

# multiply and add probabilities
dsbTableMult <- multiply_dpp(dsbTable, n=4, addIndex=TRUE)
depthSum <- sum(dsbTableMult$depth)
dsbTableMult[,prob:=depth/depthSum]

# sample from positions
#ptm <- proc.time()

countTableAll <- data.table()
for(s in seq(1,simCount))
{  print(s)
 
   # sample from DSB map
   if(minDist>0)
   {  simIndex <- sample(dsbTableMult$index, dsbCount+addCountFactor*dsbCount, prob=dsbTableMult$prob)
   } else
     simIndex <- sample(dsbTableMult$index, dsbCount, prob=dsbTableMult$prob)
   simDsbTable <- dsbTableMult[index %in% simIndex,]
   
   # get distances
   simDsbTable[,`:=`(distToNext=c(diff(position),0), chromBorder=c(diff(chrom),1))]
   if(minDist>0)
   {  vDist <- simDsbTable[chromBorder==0 & distToNext>=minDist, distToNext]
   } else
     vDist <- simDsbTable[chromBorder==0, distToNext]
   
   # create countTable
   countTable0 <- data.table(Length=sort(vDist))
   countTable <- countTable0[,.(Count=.N), by=Length]
   
   # merge with overall table
   if(s==1)
   {  countTableAll <- data.table::copy(countTable)
      setnames(countTableAll, old="Count", new="Count.1")
   } else
   {  countTableAll <- merge(countTableAll, countTable, by="Length", all=TRUE)
      setnames(countTableAll, old="Count", new=paste0("Count.",s))
   }
}
#elapsedTimeS <- proc.time() - ptm
print(paste0(round(as.numeric(elapsedTimeS["elapsed"])/60,1),"min elapsed."))

setkey(countTableAll, Length)
countTableAll <- add_average_and_sd(countTableAll, "Count", "Length", dsbCount) 

# get random distribution
pnbinomDist <- pnbinom(countTableAll$Length, 1, prob=1/(4*sum(chromLengths)/dsbCount))  # cumulative distribution of (dsbCount) random lengths 
cumTableRnd <- data.table(Length=countTableAll$Length, Freq=pnbinomDist)

if(saveTables)
  save_table(cumTableRnd, workingDirOUT, outfile=paste0(fileNameBase,"_cumTableRnd.txt"))

negBinomCumFreq <- pnbinom(maxDist, 1, prob=1/(4*sum(chromLengths)/dsbCount))
print(paste0(round(100*negBinomCumFreq,3),"% of randomly (neg. binom.) distributed DSB distances are <=", maxDist,"nt"))

avgCountSim <- as.numeric(countTableAll[Length<=maxDist,.(sum(avg))])
avgCumFreq <- as.numeric(countTableAll[Length<=maxDist,.(sum(avgFreq))])
avgFreqSd <- as.numeric(countTableAll[Length<=maxDist,.(mean(sdFreq))])
print(paste0("On average ",round(avgCountSim,1)," simulated DSB distances per cell are <=", maxDist,"nt."))
print(paste0(round(100*avgCumFreq,3),"% of simulated DSB distances are >=",minDist," and <=", maxDist,"nt (avg. sd=",100*round(avgFreqSd,4),"%; rep=",simCount,", DSB count per cell=",dsbCount,", total avg. count=",round(sum(countTableAll$avg),1),")."))
print(paste0(round(avgCumFreq/negBinomCumFreq,1),"x more in simulated DSBs than expected in neg. binomial distribution."))
print(paste0(round(percdDSB/(100*avgCumFreq),1),"x more dDSBs observed than found in this simulation (with ",percdDSB,"% dDSBs per cell), or ",round(percdDSB/(100*negBinomCumFreq),1),"x more than expected randomly (neg. binom)."))

# close log file
sink(file=NULL)
print(paste0("Log saved in ",workingDirOUT,"/",logfile))

# move files to subfolders
workingDirOUTL <- paste0(workingDirOUT,"/log_files/")
fileListL <- list.files(path=workingDirOUT, pattern=paste0(fileNameBase,".log"))
print(fileListL)
setwd(workingDirOUT)
for(file in fileListL)
{  move_file(from=file, to=paste0(workingDirOUTL,file))
}

if(saveTables)
{  workingDirOUTT <- paste0(workingDirOUT,"/tables/")
   fileListT <- list.files(path=workingDirOUT, pattern=paste0(".txt"))
   print(fileListT)
   setwd(workingDirOUT)
   for(file in fileListT)
   {  move_file(from=file, to=paste0(workingDirOUTT,file))
   }
}


