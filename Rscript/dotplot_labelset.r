# 加载所需的包
library(ggplot2)
library(ggrepel)

# 生成示例数据集
set.seed(123)
genes <- data.frame(
  Gene = paste0("Gene", 1:100),
  log2FoldChange = rnorm(100, mean = 0, sd = 1),
  pvalue = runif(100, min = 0, max = 0.05)
)

# 根据p-value和fold change进行分组
genes$Significant <- ifelse(genes$pvalue < 0.05 & abs(genes$log2FoldChange) >= 1,
                            ifelse(genes$log2FoldChange > 1, "Up", "Down"), "Stable")

# 绘制火山图
ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  # geom_text_repel(data = subset(genes, pvalue < 0.05 & abs(log2FoldChange) >= 1),
  #                 aes(label = Gene), size = 5, box.padding = unit(0.35, "lines"),
  #                 point.padding = unit(0.3, "lines"),max.overlaps = 25) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 1, lty = 4, col = "black", lwd = 0.8) +
  labs(x = "log2 (fold change)", y = "-log10 (p-value)") +
  theme(legend.position = "bottom")+
  geom_text_repel(data = genes[1:25,],
                  color = "red",
                  aes(label = Gene), size = 5, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),max.overlaps = 25)+
  geom_text_repel(data = genes[25:30,],
                  color = "green",
                  aes(label = Gene), size = 5, box.padding = unit(0.35, "lines"),show.legend = T,
                  point.padding = unit(0.3, "lines"),max.overlaps = 25)+
  theme(legend.position = "right")
