library(ggplot2)
library(scales)

args <- commandArgs(TRUE)
dd <- read.csv(args[1], sep="\t")
dd$Tool <- factor(dd$Tool, levels = c("iflas", "tama", "stringtie", "flair", "reference"))
dd$Category <- factor(dd$Category, levels = c("novel", "exon_alt_in_middle", "partial_truncated", "monoExon_truncated", "monoExon", "complete_identity", "covered"))

pdf(args[2], width=8, height=6)
p <- ggplot(dd, aes(x = Group, y = Value, fill = Category)) +
    geom_bar(stat = 'identity') +
    facet_grid(~ Tool) + theme_bw() +
    scale_y_continuous(expand = c(0.02, 0))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines")) +
    scale_fill_discrete(labels = c("Novel", "Exon Alt In Middle", "Partial Truncated", "MonoExon Truncated", "MonoExon", "Complete Identity", "Covered")) +
    scale_y_continuous(labels = comma_format())
print(p)
dev.off()