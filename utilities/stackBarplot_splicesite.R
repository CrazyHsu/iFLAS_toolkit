library(ggplot2)
args <- commandArgs(TRUE)

four_tools_splicesite <- read.csv(args[1], sep = "\t", header = F)
names(four_tools_splicesite) <- c("tools", "splicesite", "counts")
four_tools_splicesite$tools <- factor(four_tools_splicesite$tools, levels=c("B73_ref", "iflas", "tama", "stringtie", "flair"))
four_tools_splicesite$splicesite <- factor(four_tools_splicesite$splicesite, levels=c("Other", "GC-AG", "GT-AG"))

pdf(args[2])
p <- ggplot(four_tools_splicesite, aes(x=tools, y=counts, fill=splicesite, label=counts)) +
  geom_bar(position="fill", stat="identity") +
  geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  theme_bw() +
  scale_y_continuous(expand = c(0.02, 0), labels = function(x) paste0(x*100)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines"))
print(p)
dev.off()