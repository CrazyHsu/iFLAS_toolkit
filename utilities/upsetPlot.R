library(ggplot2)
library(UpSetR)

args <- commandArgs(TRUE)

four_tools_junc_intersect <- args[1]
junc_df <- read.table(four_tools_junc_intersect, sep="\t", header = T)
names(junc_df) <- c("junc_id", "iflas", "flair", "tama", "stringtie")
junc_df[junc_df >= 2] <- 1
myQuery <- list(list(query=intersects, params=list("iflas", "tama"), color="red", active=T),list(query=intersects, params=list("iflas", "tama", "stringtie"), color="green", active=T), list(query=intersects, params=list("iflas", "tama", "stringtie", "flair"), color="blue", active=T))
myQuery <- list(list(query=intersects, params=list("iflas", "tama", "stringtie"), color="#B22222", active=T), list(query=intersects, params=list("iflas"), color="#FA8072", active=T))

pdf(args[2])
upset(junc_df, sets = c("iflas", "flair", "tama", "stringtie"), order.by = c("degree", "freq"), queries=myQuery, sets.bar.color = c("maroon", "blue", "orange", "green"), text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1), mb.ratio=c(0.5,0.5))
dev.off()