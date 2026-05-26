# 微生物组β多样性PCoA分析（复现Nature Plants论文）
set.seed(123) # 固定随机种子，100%可复现
rm(list=ls())

# 加载包（你已经装好了，直接用）
library(ggplot2)
library(vegan)

# --------------------------
# 模拟论文数据：健康vs患病油菜根际微生物
# 完全不需要外部csv文件，数据直接在代码里生成
# --------------------------
# 100个样本：50个健康，50个患病
n_samples <- 100
group <- data.frame(
  sample = paste0("sample_", 1:n_samples),
  group = rep(c("健康油菜", "患病油菜"), each = 50)
)

# 模拟200个细菌OTU的丰度表
otu <- matrix(rpois(n_samples * 200, lambda = 8), nrow = n_samples, ncol = 200)
colnames(otu) <- paste0("OTU_", 1:200)
rownames(otu) <- group$sample

# 关键：让Sphingopyxis（OTU_1）在健康油菜中丰度高5倍（符合论文结论）
otu[1:50, 1] <- otu[1:50, 1] * 5

# --------------------------
# 标准微生物组分析流程
# --------------------------
# 1. 数据标准化（Hellinger转化，微生物组必做）
otu_hellinger <- decostand(otu, method = "hellinger")

# 2. 计算Bray-Curtis距离
bray_dist <- vegdist(otu_hellinger, method = "bray")

# 3. PCoA分析
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

# 4. 整理绘图数据
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$group <- group$group

# 5. 计算主成分解释度
eig_percent <- pcoa$eig / sum(pcoa$eig) * 100

# --------------------------
# 绘制论文同款PCoA图
# --------------------------
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) + # 95%置信椭圆
  labs(title = "健康vs患病油菜根际微生物群落差异",
       x = paste0("PC1 (", round(eig_percent[1], 1), "%)"),
       y = paste0("PC2 (", round(eig_percent[2], 1), "%)"),
       color = "植株状态") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# 显示图
print(p)

# 保存图到你的电脑
ggsave("pcoa_result.png", p, width = 8, height = 6, dpi = 300)

# --------------------------
# 统计检验：PERMANOVA（论文必须有）
# --------------------------
permanova_result <- adonis2(bray_dist ~ group, data = group, permutations = 999)
cat("\n===== PERMANOVA统计检验结果 =====\n")
print(permanova_result)