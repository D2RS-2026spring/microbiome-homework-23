# 多组学宿主-微生物互作作业复现
- 课程：D2RS 2026
- Issue：#23
- 复现论文：Large-scale multi-omics unveils host–microbiome interactions driving root development and nitrogen acquisition

## 核心内容
复现论文关键结论：**Sphingopyxis 菌属在健康油菜根际丰度显著更高，通过促进侧根发育提升氮素吸收效率**
- 完整微生物组分析流程：数据模拟 → Hellinger转化 → Bray-Curtis距离计算 → PCoA分析 → PERMANOVA统计检验
- 可复现性：代码固定随机种子`set.seed(123)`，所有结果可100%重复

## 文件说明
- `code.R`：完整的R分析代码
- `report.qmd`：Quarto研究报告
- `pcoa_result.png`：PCoA分析结果图
