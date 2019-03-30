#!/usr/bin/env  Rscript

NotFancy <- function(l) {
 l <- format(l, scientific = FALSE)
 parse(text=l)
}

library(ggplot2)
theme_set(theme_bw(base_size = 18))

data = read.table("counts.txt", header = T)

pdf("expression.pdf", width = 20)
p = ggplot(data, aes(x = factor(Gene), y = Mean))
p = p + geom_bar(stat = "identity", position = "dodge", fill = "red")
p = p + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.3)
p = p + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p = p + scale_y_log10(labels = NotFancy)
p = p + labs(x = "Gene", y = "Normalized counts")
print(p)
dev.off()


