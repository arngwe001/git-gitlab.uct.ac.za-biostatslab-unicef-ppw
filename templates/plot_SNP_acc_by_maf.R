#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(plyr)
library(ggplot2)
library(reshape2)
# library(gtable)
# library(grid)
# library(gridExtra)

superpop.plot.colours <- c("#8A2BE2", "#CC0066", "#6495ED", "#66CDAA", "#A52A2A", "#CDAA7D", "#66CD00", "#7AC5CD", "#CD5B45",
                           "#CDC8B1", "#00CDCD", "#CD950C", "#8B7500", "#800000", "#808000", "#008000", "#800080", "#008080", "#000080")
linetypes <- c("dashed", "solid", "dotted", "dotdash", "longdash", "twodash")

data <- read.table("${SNP_acc_report}", h=T, sep='\t')
colnames(data) <- c("${group}", "MAF<1%", "MAF 1-5%", "MAF>5%")
mdata <- melt(data, id=c("${group}"))

p <- ggplot(mdata, aes(x = variable, y = value, fill = ${group})) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=superpop.plot.colours) +
  ylab("Concordance rate") + xlab("MAF bins") +
  scale_y_continuous(breaks=c(0.0, 0.25, 0.50, 0.75, 0.95, 1)) +
  theme_bw()
p
ggsave("${SNP_acc_by_maf_plot}", height=6, width=10, units='in', dpi=150)
