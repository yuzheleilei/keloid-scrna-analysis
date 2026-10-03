# keloid-scrna-analysis
Single-cell RNA-seq analysis of keloid scars comparing **pruritic (AK, itchy)** vs **non-pruritic (IAK, non-itchy)** lesions, with matched mature-stage (MS) and normal skin (NS) controls. 
| 软件 | 版本 | 用途 |
| --- | --- | --- |
| Cell Ranger | v6.1.0 | 10X BCL→FASTQ 比对、UMI 定量 |
| R | v4.0.3 | 全部分析运行环境 |
| Seurat | v3.1.1 | QC、归一化、降维聚类、FindMarkers（Wilcoxon）、AddModuleScore |
| DoubletFinder | v2.0.3 | 双细胞剔除 |
| Harmony |v1.0| 批次整合（前 30 PCs，batch=数据集来源） |
| monocle3 | v1.0.0 | 拟时序 |
| CytoTRACE2 | v1.1.0 | 发育潜能打分 |
| CellChat | v1.6.0 | 细胞通讯 |
