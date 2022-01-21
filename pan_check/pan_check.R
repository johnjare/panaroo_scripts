args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], header = TRUE)

jpeg("pan_check.jpeg", width = nrow(df)/4)
heatmap(t(as.matrix(df[,-1])), col = c("white", "black"), labCol = FALSE, xlab = "Genes")
legend(x="right", legend=c("Present", "Absent"),fill=c("black", "white"))
dev.off()
