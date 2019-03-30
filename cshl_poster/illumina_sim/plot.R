#!/usr/bin/Rscript

library(ggplot2)
library(grid)
library(gridExtra)

# -----------------------------------------------------------------------------
theme_set(theme_bw(base_size = 16))
png("random.png", width = 600, height = 580, res = 100)

data = read.table("random.txt", header = T)
tp = subset(data, subset = data$Index == "TP")
fp = subset(data, subset = data$Index == "FP")

df = data.frame(
    Mapper = c("BWA", "SRTk"),
    means = c(mean(subset(tp, subset = tp$Mapper == "BWA")$Value),
              mean(subset(tp, subset = tp$Mapper == "SRTk")$Value)),
    sds = c(sd(subset(tp, subset = tp$Mapper == "BWA")$Value),
            sd(subset(tp, subset = tp$Mapper == "SRTk")$Value))
    )

q1 = qplot(x = factor(Mapper), y = means, data = df, fill = factor(Mapper))
q1 = q1 + geom_bar(data = df, stat = "identity", position = "dodge")
q1 = q1 + geom_errorbar(data = df, aes(ymin=means-sds, ymax=means+sds), width =
0.5)
q1 = q1 + geom_jitter()
q1 = q1 + labs(x = "\nTool used with LUMPY", y = "Breakpoint pairs correctly identified\n")
q1 = q1 + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12),
                plot.margin = unit(c(3,1,1,1), "mm"))
q1 = q1 + guides(fill=F)

df = data.frame(
    Mapper = c("BWA", "SRTk"),
    means = c(mean(subset(fp, subset = fp$Mapper == "BWA")$Value),
              mean(subset(fp, subset = fp$Mapper == "SRTk")$Value)),
    sds = c(sd(subset(fp, subset = fp$Mapper == "BWA")$Value),
            sd(subset(fp, subset = fp$Mapper == "SRTk")$Value))
    )


q2 = qplot(x = factor(Mapper), y = means, data = df, fill = factor(Mapper))
q2 = q2 + geom_bar(data = df, stat = "identity", position = "dodge")
q2 = q2 + geom_errorbar(data = df, aes(ymin=means-sds, ymax=means+sds), width =
0.5)
q2 = q2 + geom_jitter()
q2 = q2 + labs(x = "\nTool used with LUMPY", y = "Breakpoint pairs incorrectly identified\n")
q2 = q2 + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12),
                plot.margin = unit(c(3,1,1,1), "mm"))
q2 = q2 + guides(fill=F)

blank<-rectGrob(gp=gpar(col="white")) 

grid.arrange(q1, blank, q2, nrow = 1, ncol=3, widths=c(1, 0.25, 1), top = "(a)")

dev.off()

# -----------------------------------------------------------------------------
theme_set(theme_bw(base_size = 16))
png("repeats.png", width = 600, height = 580, res = 100)

data = read.table("repeats.txt", header = T)
tp = subset(data, subset = data$Index == "TP")
fp = subset(data, subset = data$Index == "FP")

df = data.frame(
    Mapper = c("BWA", "SRTk"),
    means = c(mean(subset(tp, subset = tp$Mapper == "BWA")$Value),
              mean(subset(tp, subset = tp$Mapper == "SRTk")$Value)),
    sds = c(sd(subset(tp, subset = tp$Mapper == "BWA")$Value),
            sd(subset(tp, subset = tp$Mapper == "SRTk")$Value))
    )

q1 = qplot(x = factor(Mapper), y = means, data = df, fill = factor(Mapper))
q1 = q1 + geom_bar(data = df, stat = "identity", position = "dodge")
q1 = q1 + geom_errorbar(data = df, aes(ymin=means-sds, ymax=means+sds), width =
0.5)
q1 = q1 + geom_jitter()
q1 = q1 + labs(x = "\nTool used with LUMPY", y = "Breakpoint pairs correctly identified\n")
q1 = q1 + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12),
                plot.margin = unit(c(3,1,1,1), "mm"))
q1 = q1 + guides(fill=F)

df = data.frame(
    Mapper = c("BWA", "SRTk"),
    means = c(mean(subset(fp, subset = fp$Mapper == "BWA")$Value),
              mean(subset(fp, subset = fp$Mapper == "SRTk")$Value)),
    sds = c(sd(subset(fp, subset = fp$Mapper == "BWA")$Value),
            sd(subset(fp, subset = fp$Mapper == "SRTk")$Value))
    )


q2 = qplot(x = factor(Mapper), y = means, data = df, fill = factor(Mapper))
q2 = q2 + geom_bar(data = df, stat = "identity", position = "dodge")
q2 = q2 + geom_errorbar(data = df, aes(ymin=means-sds, ymax=means+sds), width =
0.5)
q2 = q2 + geom_jitter()
q2 = q2 + labs(x = "\nTool used with LUMPY", y = "Breakpoint pairs incorrectly identified\n")
q2 = q2 + theme(axis.text=element_text(size=12),
                axis.title=element_text(size=12),
                plot.margin = unit(c(3,1,1,1), "mm"))
q2 = q2 + guides(fill=F)

blank<-rectGrob(gp=gpar(col="white")) 

grid.arrange(q1, blank, q2, nrow = 1, ncol=3, widths=c(1, 0.25, 1), top = "(b)")

dev.off()
