#!/usr/bin/Rscript

library(ggplot2)
library(reshape)

theme_set(theme_bw(base_size = 13))

pdf("tumor.pdf", width = 18)
data = read.table("input.txt", header = T)

data.m = melt(data)
names(data.m) = c("Gene1.Gene2", "variable", "Support")

q = qplot(x = factor(variable), y = Gene1.Gene2, data = data.m)
q = q + geom_tile(aes(fill = Support))
q = q + scale_fill_gradient(low="white", high="red")
q = q + labs(x = "SampleId", y = "Genes involved in fusion")
print(q)

dev.off()
