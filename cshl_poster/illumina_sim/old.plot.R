#!/usr/bin/Rscript

library(ggplot2)
library(grid)
library(gridExtra)

theme_set(theme_bw(base_size = 18))

# -----------------------------------------------------------------------------
png("random.png")

data = read.table("random.txt", header = T)
tp = subset(data, subset = data$Index == "TP")
fp = subset(data, subset = data$Index == "FP")

q1 = qplot(x = factor(Mapper), y = Value, data = tp)
q1 = q1 + geom_boxplot(aes(fill = factor(Mapper))) + ylim(0,6300) 
q1 = q1 + geom_jitter()
q1 = q1 + labs(x = "Pipeline", y = "Breakpoints correctly identified")
q1 = q1 + guides(fill=F)

q2 = qplot(x = factor(Mapper), y = Value, data = fp)
q2 = q2 + geom_boxplot(aes(fill = factor(Mapper))) + ylim(0,70)
q2 = q2 + geom_jitter()
q2 = q2 + labs(y = "Breakpoints incorrectly identified")
q2 = q2 + guides(fill=F)

grid.arrange(q1, q2, nrow = 1, ncol=2, top = "(a)")

dev.off()

# -----------------------------------------------------------------------------
theme_set(theme_bw(base_size = 18))

png("repeats.png")

data = read.table("repeats.txt", header = T)
tp = subset(data, subset = data$Index == "TP")
fp = subset(data, subset = data$Index == "FP")

q3 = qplot(x = factor(Mapper), y = Value, data = tp)
q3 = q3 + geom_boxplot(aes(fill = factor(Mapper))) + ylim(0,5700)
q3 = q3 + geom_jitter()
q3 = q3 + labs(x = "Pipeline", y = "TP breakpoint count")
q3 = q3 + guides(fill=F)

q4 = qplot(x = factor(Mapper), y = Value, data = fp)
q4 = q4 + geom_boxplot(aes(fill = factor(Mapper))) + ylim(0,120)
q4 = q4 + geom_jitter()
q4 = q4 + labs(x = "Pipeline", y = "FP breakpoint count")
q4 = q4 + guides(fill=F)

grid.arrange(q3, q4, nrow = 1, ncol=2, 
    top = "(b)")

dev.off()
