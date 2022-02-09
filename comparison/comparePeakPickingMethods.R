## Part 2 - compare results of peak picking

rm(list=ls())
cat("\014")  
set.seed(2022)
setwd("D:/PeakBot_Data/PHM_comparison")
#setwd("~/Documents/ucloud/PeakBot/comparison")

figurePlots = list()



#############################
### implement some common functions
tic<-function(){ticStart__ptm<<-proc.time(); return(ticStart__ptm[3])}
toc<-function(unit="min"){return((proc.time()[3]-ticStart__ptm[3])/list(min=60, sec=1, hours=60*60, days=24*60*60)[unit][[1]])}
tocP<-function(unit="min"){print(sprintf("Time elapsed since last tic: %.1f %s", toc(), unit))}
source("./ext__exportAsFeatureML.R")

library(ggplot2)
library(grid)
library(gridExtra)
library(VennDiagram)
library(plotly)
library(cowplot)
library(MSnbase)
library(magrittr)
library(dplyr)
library(ggrepel)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")




checkForInt <- function(rts, ints, mzs, mzmin, mzmax, rtmin, rtmax, minInt){
  isPresent = FALSE
  maxInds = c(-1, -1, -1)
  for(rti in 1:length(rts)){
    if(rtmin <= rts[rti] && rts[rti] <= rtmax){
      a = which(mzmin <= mzs[[rti]] & mzs[[rti]] <= mzmax & minInt <= ints[[rti]])
      if(length(a) > 0){
        temp = a[which.max(ints[[rti]][a])]
        
        if(maxInds[1] == -1 || maxInds[2] < ints[[rti]][temp]){
          maxInds = c(rti, ints[[rti]][temp], temp)
        }
        isPresent = TRUE
      }
    }
  }
  apexRT = -1
  if(isPresent){
    apexRT = rts[maxInds[1]]
  }
  return(list(isPresent = isPresent, apexRT = apexRT))
}



## Comparison of detected chromatographic peaks with the different methods
files = c("05_EB3388_AOH_p_0","06_EB3389_AOH_p_10",
          "07_EB3390_AOH_p_20","08_EB3391_AOH_p_60", 
          "16_EB3392_AME_p_0","17_EB3393_AME_p_10",
          "18_EB3394_AME_p_20","19_EB3395_AME_p_60")
methods = c("PeakBot", "XCMS", "MSDial", "peakOnly")
commonIDs = list()
commonIDsFoundIn = data.frame(id=c(), avMZ=c(), avRT=c(), method=c())
maxPPM = 15
maxRT  = 0.15 * 60

venns = list()
plots = list()
allDat = list()
file = files[4]
for(file in files){
  message(sprintf("processing sample '%s'", file))
  
  res <- readMSData(paste0("./../PHM/", file, ".mzXML"), msLevel=1, centroided=FALSE)
  rts = rtime(res)
  ints <- intensity(res)
  mzs <- mz(res)
  
  resMethods = list()
  plotDat = data.frame(id=c(), rt=c(), mz=c(), area=c())
  for(method in methods){
    dat <- switch(method,
                  "PeakBot" = read.table(paste0("./PeakBot/", file, "_positivePeakBot.tsv"), sep="\t", quote="", header=TRUE), 
                  "XCMS"    = read.table(paste0("./XCMS/", file, ".mzXML_pos.tsv"), sep="\t", quote="", header=TRUE), 
                  "MSDial"  = read.table(paste0("./MSDial/", file, "_6e4Thres.txt"), sep="\t", quote="", header=TRUE), 
                  "peakOnly"  = read.table(paste0("./peakOnly/", file, ".csv"), sep=",", quote="", header=TRUE))
    
    if(method=="XCMS"){
      dat = cbind(Num=row.names(dat), dat)
      
      use = rep(FALSE, nrow(dat))
      for(rowi in 1:nrow(dat)){
        ## although XCMS reports min/max mz values, these are only valid in the centroid mode data and need to be converted
        u = checkForInt(rts, ints, mzs, dat[rowi,"mzmin"] * (1 - 5 / 1E6), dat[rowi,"mzmax"] * (1 + 5 / 1E6), dat[rowi,"rtmin"], dat[rowi,"rtmax"], 1E5)
        use[rowi] = u$isPresent
      }
      message(paste0("    Removing ", sum(!use), " features from XCMS evaluation with ", nrow(dat), " features due to a too low abundance"))
      dat = dat[use, ]
    }
    if(method=="MSDial"){
      colnames(dat)[colnames(dat)=="PeakID"] = "Num"
      dat[,"RT..min."] = dat[,"RT..min."] * 60
      colnames(dat)[colnames(dat)=="RT..min."] = "RT..sec."
      dat[,"RT.left.min."] = dat[,"RT.left.min."] * 60
      colnames(dat)[colnames(dat)=="RT.left.min."] = "RT.left.sec."
      dat[,"RT.right..min."] = dat[,"RT.right..min."] * 60
      colnames(dat)[colnames(dat)=="RT.right..min."] = "RT.right..sec."
      
      use = rep(FALSE, nrow(dat))
      for(rowi in 1:nrow(dat)){
        u = checkForInt(rts, ints, mzs, dat[rowi,"Precursor.m.z"] * (1 - 5 / 1E6), dat[rowi,"Precursor.m.z"] * (1 + 5 / 1E6), dat[rowi,"RT.left.sec."], dat[rowi,"RT.right..sec."], 1E5)
        use[rowi] = u$isPresent
      }
      message(paste0("    Removing ", sum(!use), " features from MSDial evaluation with ", nrow(dat), " features due to a too low abundance"))
      dat = dat[use, ]
    }
    if(method=="peakOnly"){
      colnames(dat)[1] = "Num"
      colnames(dat)[ncol(dat)] = "area"
      dat[,"rtmin"] = dat[,"rtmin"] * 60 * 60
      dat[,"rtmax"] = dat[,"rtmax"] * 60 * 60
      dat = cbind(dat, rt=-1)
      
      use = rep(FALSE, nrow(dat))
      for(rowi in 1:nrow(dat)){
        u = checkForInt(rts, ints, mzs, dat[rowi,"mz"] * (1 - 5 / 1E6), dat[rowi,"mz"] * (1 + 5 / 1E6), dat[rowi,"rtmin"], dat[rowi,"rtmax"], 1E5)
        use[rowi] = u$isPresent
        if(use[rowi]){
          dat[rowi, "rt"] = u$apexRT
        }
      }
      message(paste0("    Removing ", sum(!use), " features from peakOnly evaluation with ", nrow(dat), " features due to a too low abundance"))
      dat = dat[use, ]
    }
    
    dat = switch(method, 
                 "PeakBot" = dat[, c("Num", "RT", "MZ", "PeakArea")],
                 "XCMS"    = dat[, c("Num", "rt", "mz", "intb")],
                 "MSDial"  = dat[, c("Num", "RT..sec.", "Precursor.m.z", "Area")],
                 "peakOnly"  = dat[, c("Num", "rt", "mz", "area")])
    colnames(dat) = c("Num", "rt", "mz", "area")
    row.names(dat) = paste0("Feat", dat$Num)
    
    dat = dat[120 <= dat[,"rt"] & dat[,"rt"] <= 750,]
    
    dat = cbind(method=method, dat)
    for(otherMethod in methods){
      dat = cbind(dat, newCol = "")
      colnames(dat)[ncol(dat)] = paste0("FoundIn_", otherMethod)
    }
    for(otherMethod in methods){
      dat = cbind(dat, newCol = dat[,"area"])
      colnames(dat)[ncol(dat)] = paste0("Area_", otherMethod)
      
      dat = cbind(dat, newCol = 0)
      colnames(dat)[ncol(dat)] = paste0("ClosestRT_", otherMethod)
      dat = cbind(dat, newCol = 0)
      colnames(dat)[ncol(dat)] = paste0("ClosestMZ_", otherMethod)
    }
    resMethods[[method]] = dat
    plotDat = rbind(plotDat, dat)
  }
  
  message("  | .. peaks per data analysis method")
  print(lapply(resMethods, nrow))
  
  p <- ggplot(data=plotDat, mapping=aes(x=rt, y=mz, colour=method)) + geom_point(alpha=0.3) + theme_minimal()
  print(p)
  plots[[file]] = p
  #print(ggplot(data=plotDat, mapping=aes(x=rt, y=mz, colour=method)) + geom_point() + facet_wrap(~method) + theme_minimal())
  #ggplotly(ggplot(data=plotDat, mapping=aes(x=rt, y=mz, colour=method)) + geom_point(alpha=0.3) + theme_minimal())
  
  commonID = 1
  for(method in methods){
    commonIDs[[method]] = c()
  }
  
  found = list()
  for(method in methods){
    for(rowi in 1:nrow(resMethods[[method]])){
      mz = resMethods[[method]][rowi, "mz"]
      rt = resMethods[[method]][rowi, "rt"]
      area = resMethods[[method]][rowi, "area"]
      
      if(resMethods[[method]][rowi, paste0("FoundIn_", method)] != "X"){
        resMethods[[method]][rowi, paste0("FoundIn_", method)] = "X"
        commonIDs[[method]] = c(commonIDs[[method]], commonID)
        
        foundIn = c(method)
        
        mzs = c(mz)
        rts = c(rt)
        
        for(oMethod in methods){
          if(oMethod != method){
            js = which(abs(resMethods[[oMethod]][,"mz"] - mz)/mz*1E6 <= maxPPM & abs(resMethods[[oMethod]][,"rt"] - rt) <= maxRT)
            if(length(js) > 0){
              resMethods[[oMethod]][js, paste0("FoundIn_", oMethod)] = "X"
              resMethods[[oMethod]][js, paste0("FoundIn_", method)] = "X"
              resMethods[[method]][rowi, paste0("FoundIn_", oMethod)] = "X"
              
              resMethods[[oMethod]][js, paste0("Area_", method)] = area
              resMethods[[method]][rowi, paste0("Area_", oMethod)] = resMethods[[oMethod]][js[1], "area"]
              
              commonIDs[[oMethod]] = c(commonIDs[[oMethod]], commonID)
              foundIn = c(foundIn, oMethod)
              
              mzs = c(mzs, resMethods[[oMethod]][js,"mz"])
              rts = c(rts, resMethods[[oMethod]][js,"rt"])
            }
          }
        }
        
        toPut = paste0(sort(foundIn), collapse="_")
        if(!(toPut %in% names(found))){
          found[[toPut]] = resMethods[[method]][rowi,]
        }else{
          found[[toPut]] = rbind(found[[toPut]], resMethods[[method]][rowi,])
        }
      }
      commonIDsFoundIn = rbind(commonIDsFoundIn, data.frame(id=commonID, avMZ=mean(mzs), avRT=mean(rts), method=paste0(foundIn, collapse=", ")))
      commonID = commonID + 1
    }
  }
  if(file == "08_EB3391_AOH_p_60"){
    message(paste0("   Results of ", file, ":"))
    print(lapply(found, nrow))
    for(n in names(found)){
      uset = 1:nrow(found[[n]])
      if(length(uset) > 50){
        uset = sample(uset, 50)
      }
      exportAsFeatureXML(found[[n]][uset,], paste0("manual__", file, "___foundIn_", n, ".featureML"), cid = "Num", crt = "rt", cmz = "mz")
      exportAsFeatureXML(found[[n]], paste0(file, "___foundIn_", n, ".featureML"), cid = "Num", crt = "rt", cmz = "mz")
    }
    
    tab = resMethods[["PeakBot"]][resMethods[["PeakBot"]][,"FoundIn_PeakBot"]=="X" & resMethods[["PeakBot"]][,"FoundIn_XCMS"]=="X",]
    tab = tab[!is.na(tab[,"Area_PeakBot"]) & !is.na(tab[,"Area_XCMS"]),]
    
    p <- ggplot(mapping=aes(x = log2(Area_PeakBot), y = log2(Area_XCMS), colour="all")) + theme_minimal()
    p <- p + geom_point(data = tab, colour="slategrey", alpha=0.2)
    p <- p + stat_smooth(data = tab, method = "lm", col = "slategrey")
    fit <- lm(Area_XCMS ~ Area_PeakBot, data=tab)
    summary(fit)
    
    qs = quantile(tab$Area_PeakBot/tab$Area_XCMS, probs=c(0.05, 0.95))
    tabr = tab[qs[1] <= tab$Area_PeakBot/tab$Area_XCMS & tab$Area_PeakBot/tab$Area_XCMS <= qs[2],]
    
    p <- p + geom_point(data = tabr, mapping=aes(x = log2(Area_PeakBot), y = log2(Area_XCMS), colour="restricted"), colour="firebrick", alpha=0.2)
    p <- p + stat_smooth(data = tabr, method = "lm", col = "firebrick")
    fitr <- lm(Area_XCMS ~ Area_PeakBot, data=tabr)
    summary(fitr)
    
    p <- p + ggtitle("Comparison of PeakBot and XCMS peak areas") + labs(caption=sprintf("PHM dataset\ngrey: all data points, R2 of %.3f; red: lower/upper 5%% with highest deviations removed, R2 of %.3f", summary(fit)$r.squared, summary(fitr)$r.squared))
    p <- p + xlab("Peak area determined with PeakBot (log2)") + ylab("Peak area determined with XCMS (log2)")
    print(p)
    
    ggsave("./fig_lmPeakAreas.png", plot=p, width=6, height=4, dpi=300)
    figurePlots[["lm"]] = p
    
    message("    factors between PeakBot and XCMS determined peak areas:")
    quantile(tab$Area_PeakBot/tab$Area_XCMS, probs=c(0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
  }
  
  p = ggplot(data=commonIDsFoundIn, mapping=aes(x=avRT, y=avMZ, colour=method)) + geom_point() + facet_wrap(~method) + ggtitle("Comarison of methods", sub=sprintf("Sample: '%s'", file)) + xlab("Retention time (sec)") + ylab("m/Z") + theme(legend.position="bottom")
  #plots[[file]] = p
  print(p)
  
  temp = list()
  for(method in methods){
    temp[[method]] = commonIDs[[method]]
    if(!(method %in% names(allDat))){
      allDat[[method]] = c()
    }
    allDat[[method]] = c(allDat[[method]], paste0(file, "__", commonIDs[[method]]))
  }
  venn = venn.diagram(x = temp, filename=NULL, main=file, euler.d=FALSE, scaled=FALSE, print.mode=c("raw", "percent"), col=NA, fill=c("Firebrick", "Dodgerblue", "Olivedrab", "slategrey")[1:length(methods)], cex=0.7, sigdigs=2)
  venn = gList(venn, grid.text(paste0("PeakBot\ntotal: ", round(length(temp[["PeakBot"]]) / length(unique(c(temp[["PeakBot"]], temp[["XCMS"]], temp[["MSDial"]], temp[["peakOnly"]])))*100, 2), "%"), x=0.01, y=0.05, just="left", gp=gpar(col="slategrey", fontsize=12)))
  venns[[file]] = venn
  grid.newpage()
  grid.draw(venns[[file]])
}

temp = list()
for(file in files){
  temp[[file]] = plots[[file]] + theme(legend.position="None")
}

legend <- cowplot::get_legend(plots[[file]])
temp[["legend"]] = legend
p <- plot_grid(plotlist = temp, scale=0.95, labels="auto")
print(p)
save_plot("./fig_comparisonsFeatureMap.png", plot=p, base_width =10, base_height =10)


temp = list()
for(file in files){
  temp[[file]] = venns[[file]]
}

venn = venn.diagram(x = allDat, filename=NULL, main="average", euler.d=FALSE, scaled=FALSE, print.mode=c("percent"), col=NA, fill=c("Firebrick", "Dodgerblue", "Olivedrab", "slategrey")[1:length(methods)], cex=0.7, sigdigs=2)
venn = gList(venn, grid.text(paste0("PeakBot\ntotal: ", round(length(allDat[["PeakBot"]]) / length(unique(c(allDat[["PeakBot"]], allDat[["XCMS"]], allDat[["MSDial"]], allDat[["peakOnly"]])))*100, 2), "%"), x=0.01, y=0.05, just="left", gp=gpar(col="slategrey", fontsize=12)))
temp[["average"]] = venn

p <- plot_grid(plotlist = temp, scale=0.95, labels="auto")
print(p)
save_plot("./fig_comparisonsVenn.png", plot=p, base_width =10, base_height =10)
save_plot("./fig_comparisonsVennAverage.png", plot=temp[["average"]], base_width =6, base_height =6)
save_plot("./fig_comparisonsVenn08_EB3391_AOH_p_60.pdf", plot=temp[["08_EB3391_AOH_p_60"]], base_width =6, base_height =6)
  






rsd <- function(vals){
  if(sum(!is.na(vals))<=1){
    return(-1)
  }
  return(sd(vals)/mean(vals))
}

dat = read.delim("./PeakBot/WheatEar_Features.tsv")
dat = dat[8:16]

grps = rep(c("lvl1", "lvl1x2", "lvl1x4"), each=3)

df = data.frame(rowi = c(),
                grp = c(), 
                rsd = c())

for(grp in sort(unique(grps))){
  temp = apply(dat, 1, function(x){
    vals = c(x[grps==grp])
    vals = vals[vals>0]
    if(length(vals)<2){
      return(0)
    }else{
      return(rsd(vals))
    }
  })
  use = temp > 0 
  df = rbind(df, 
             data.frame(rowi = c(1:nrow(dat))[use], grp = rep(grp, times = sum(use)), rsd = temp[use]))
}


p <- ggplot(data = df, mapping = aes(x=rsd, group=grp)) + theme_minimal()
dfs = df %>% group_by(grp) %>% summarise(meanRSD = mean(rsd), medianRSD = median(rsd), q90RSD = quantile(rsd, probs=0.9))
p <- p + facet_wrap(~grp)
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = meanRSD), colour="grey")
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = medianRSD), colour="grey")
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = q90RSD), colour="grey")
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("avg: %.1f%%", 100*meanRSD), x = 0.22, y = 8000, hjust=0))
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("med: %.1f%%", 100*medianRSD), x = 0.22, y = 7200, hjust=0))
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("P90: %.1f%%", 100*q90RSD), x = 0.22, y = 6400, hjust=0))
p <- p + geom_histogram()
p <- p + ggtitle("Relative standard deviation of replicates") + labs(caption="WheatEar dataset\nlvl1 = 2\U003BCl, lvl1x2 = 4\U003BCl, lvl1x4 = 8 \U003BCl injection volumne\nreplicates with 0 peak area omitted; avg: average, med: median, P90: 90% percentile") + xlab("Relative standard deviation") + ylab("Number of features")
p <- p + xlim(c(0, 0.75))
print(p)
figurePlots[["rsd"]] = p

ggsave("./fig_rsds.png", plot = p, width = 6, height=4, dpi = 300)


p <- plot_grid(plotlist = list(figurePlots[["rsd"]], figurePlots[["lm"]]), ncol=1, scale=1, labels="auto")
print(p)
save_plot("./fig_quantification.png", plot=p, base_width =6, base_height =8, dpi = 300)


