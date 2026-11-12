# FACHB-416 Paper1 Round16：S0 候选 DTU 定位预测整合报告

生成时间：2026-07-05 11:52:01

## 输入与边界

- 本轮在用户明确同意后，将 S0 小批量蛋白序列提交到 DTU 在线服务；交付文件不包含原始蛋白序列。
- 使用工具：DeepTMHMM 1.0、SignalP 6.0、TargetP 2.0（Plant 模式）和 DeepLoc 2.1。
- 本报告只把定位结果作为生信预测证据；不能写成已经实验验证的分泌、跨膜、细胞壁或细胞器定位。
- 功能层证据沿用前序 dbCAN5/Pfam/InterPro/motif 的 S0 证据表，定位层用于提升候选优先级和机制叙事。

## DTU 任务链接

- DeepTMHMM：https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDE000304C20E9C081DF
- SignalP6：https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDED00304C381A6DDA90
- TargetP2：https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDF100304C3F1BDE4213
- DeepLoc2：https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDFF00304C486001E21F

## 核心结果

- SignalP6 预测 SP：5/27。
- TargetP2 预测 SP：5/27；mTP：3；cTP：1。
- DeepTMHMM 检出 signal：6/27；检出至少 1 个 TMR：9/27。
- DeepLoc2 含 Extracellular 定位：7/27；含 Golgi apparatus 定位：10/27。

## 支持 Paper1 主线的候选层级

### A1/B1：分泌或胞外强定位线索

g1615.t1(GT49), g5532.t1(GT49), g651.t1(GT49), g12848.t1(GT32), g12483.t1(-)

解释：这些候选同时获得 SignalP6、DeepTMHMM、TargetP2 或 DeepLoc2 中多条分泌/胞外证据。特别是 GT49 候选 g651.t1、g1615.t1、g5532.t1 以及 GT32 弱 dbCAN 候选 g12848.t1，可作为“扩张/保留的表面糖链或胞外基质相关模块”的优先候选；S0C 的 g12483.t1 则更适合写作细胞壁糖蛋白成熟或胞外基质蛋白加工相关候选。

### A2/B2：膜系统或高尔基强定位线索

g11399.t1(GT32), g9209.t1(GT49), g9850.t1(GT32), g11586.t1(-), g12701.t1(-), g12816.t1(-), g13195.t1(-), g7484.t1(-)

解释：这些候选含跨膜段或 DeepLoc2 高尔基/膜系统定位，适合支撑“糖链装配、糖前体转运和分泌通路辅助模块被保留或强化”的叙事。其中 GT32/GT49 候选 g11399.t1、g9209.t1、g9850.t1 可作为膜系统糖基转移酶候选；S0B 多跨膜候选 g7484.t1、g11586.t1、g12701.t1、g12816.t1、g13195.t1 更适合作为前体转运或膜系统辅助候选。

### C 层：细胞器定位或表面机制需降级

g12482.t1, g2036.t1, g4696.t2, g7252.t1, g12484.t2, g1702.t1

解释：这部分可能参与碳氮代谢、能量或细胞器状态，而不是直接定位到细胞表面。可保留为间接机制或备用候选，不建议放入 EPS/细胞壁主证据链。

## 对 Paper1 假说的增益

本轮结果把原来“S0 家族有糖基化/细胞壁相关功能”的证据推进到“其中一部分候选还具备合理的分泌、胞外、跨膜或高尔基定位预测”。这对第一篇的主线有两个作用：

1. 支撑长期实验室培养后，少数扩张或保留模块并非随机，而是集中在生长维持、膜系统糖基化、胞外/表面结构更新和形态响应相关过程。
2. 将“表面糖链/EPS/细胞壁糖蛋白”从泛泛功能注释，升级为可以按候选等级排序的机制模型：分泌型 GT/壁蛋白候选负责胞外糖链或基质输出，膜系统 GT 与多跨膜候选负责糖链装配和前体转运，细胞器定位候选作为能量或代谢背景调节层。

## 建议写法

- 可以写：S0 候选中存在一组 dbCAN 支持的 GT32/GT49/GT47/GT10 糖基转移酶候选，其中若干成员同时获得 SignalP6、DeepTMHMM、TargetP2 和 DeepLoc2 的分泌/胞外或膜系统定位预测，提示其可能参与表面糖链、EPS 或细胞壁相关糖基化过程。
- 不建议写：这些基因已经被证明定位于细胞壁或直接合成 EPS。
- 主图建议：把 Round16 表中的“功能证据等级、定位等级、S0 家族变化方向、DeepLoc 定位、TMR 数量”合成一个多轨热图或 circos 内圈，用 A1/A2/B1/B2 候选做醒目标注。

## 输出文件

- 整合明细表：`FACHB416_paper1_round16_S0_DTU_localization_integrated_zh_20260705.tsv`
- 工具分类计数：`FACHB416_paper1_round16_S0_DTU_localization_tool_counts_zh_20260705.tsv`
- R 作图源数据：`FACHB416_paper1_round16_S0_DTU_localization_figure_source_data_20260705.tsv`
- 本报告：`FACHB416_paper1_round16_S0_DTU_localization_summary_zh_20260705.md`
- 运行日志：`FACHB416_paper1_round16_S0_DTU_localization_run_log_20260705.md`
