
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggrepel)

rsd <- function(vals){
  if(sum(!is.na(vals))<=1){
    return(-1)
  }
  return(sd(vals)/mean(vals))
}

dat = read.delim("~/LCHRMS-GPU/peakbot_example/Data/WheatEar_positivedetectedFeatures.tsv")
dat = dat[8:16]

grps = rep(c("lvl1", "lvl1x2", "lvl1x4"), each=3)

df = data.frame(rowi = rep(-1, nrow(dat)*length(unique(grps))),
               grp = factor(0, levels = unique(grps)), 
               rsd = rep(-1, nrow(dat)*length(unique(grps))))

ug <- unique(grps)
cur = 1
for(rowi in 1:nrow(dat)){
  vals = dat[rowi,]
  for(grp in ug){
    df[cur, 1] = rowi
    df[cur, 2] = grp
    vs = t(vals[grps==grp])
    vs = vs[!is.na(vs) & vs>0]
    if(length(vs)>1){
      df[cur, 3] = rsd(vs)
      cur = cur + 1
    }
  }
}
df = df[1:cur,]

p <- ggplot(data = df, mapping = aes(x=rsd, group=grp)) + theme_minimal()
dfs = df %>% group_by(grp) %>% summarise(meanRSD = mean(rsd), medianRSD = median(rsd), q90RSD = quantile(rsd, probs=0.9))
p <- p + facet_wrap(~grp)
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = meanRSD), colour="grey")
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = medianRSD), colour="grey")
p <- p + geom_vline(data = dfs, mapping=aes(xintercept = q90RSD), colour="grey")
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("avg: %.1f%%", 100*meanRSD), x = 0.22, y = 8000, hjust=0))
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("med: %.1f%%", 100*medianRSD), x = 0.22, y = 7200, hjust=0))
p <- p + geom_label(data = dfs, mapping=aes(label=sprintf("Q90: %.1f%%", 100*q90RSD), x = 0.22, y = 6400, hjust=0))
p <- p + geom_histogram()
p <- p + ggtitle("Relative standard deviation of replicates") + xlab("Relative standard deviation") + ylab("Number of features")
p <- p + xlim(c(0, 0.75))
print(p)

ggsave("~/LCHRMS-GPU/rsds.png", plot = p, width = 6, height=4, dpi = 300)
