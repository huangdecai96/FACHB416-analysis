# Paper1 Round11 dbCAN5 结果交接稿

## 可写入方法

为系统评估 S0 表面糖链/细胞壁糖基化相关候选基因的 CAZyme 证据，并建立全蛋白组背景，本研究在服务器 `/PROJECT_ROOT/Gene/paper1_round11_annotation_20260704` 下建立独立 dbCAN5 环境，使用 dbCAN 5.2.9 和 DIAMOND 2.2.2。数据库目录为 `db/dbcan3`，包含 `CAZy.dmnd`、`dbCAN.hmm`、`dbCAN-sub.hmm` 及底物映射文件。注释采用 `run_dbcan CAZyme_annotation --mode protein --methods diamond,hmm,dbCANsub --threads 64`，先对 27 个 S0 蛋白进行试跑，确认完整三方法可用后，再对 13,767 个最长异构体蛋白进行全蛋白组 CAZyme 背景重建。

## 可写入结果

S0 的 27 个候选蛋白中，13 个获得至少 1 个 dbCAN5 证据；其中 12 个获得 2 个以上工具支持，作为高可信 CAZyme 候选。高可信 S0 CAZyme 全部属于 GT 类糖基转移酶，主要包括 GT32、GT49、GT47 和 GT10。另有 `g12848.t1` 仅获得 DIAMOND-only GT32 证据，应作为弱候选放入补表，不应与 HMM/dbCAN-sub 支持的主证据等同。

全蛋白组背景中，750/13,767 个蛋白获得至少 1 个 CAZyme 证据，302/13,767 个蛋白获得 2 个以上工具支持。以高可信 CAZyme 为口径，S0 的 CAZyme 比例为 12/27，而全蛋白组背景为 302/13,767；Fisher 单侧检验显示 S0 显著富集 CAZyme 候选（OR=37.34，P=1.29e-13）。以高可信 GT 类为口径，S0 为 12/27，全蛋白组背景为 120/13,767，富集更强（OR=101.33，P=1.71e-18）。

## 建议叙事

这一结果支持 Paper1 的核心假说：长期实验室培养株并非均匀保留所有环境适应相关基因，而是在大规模基因家族缩减背景下，保留或扩张了一组与生长维持、细胞表面糖基化、形态/细胞壁结构调节相关的候选家族。S0 中 GT32/GT49/GT47/GT10 的集中出现，提示这些基因可作为“实验室营养适应型”中形态响应和表面结构调节能力被维持或增强的分子证据。表述边界建议写作“表面糖链/细胞壁糖基化/EPS 相关候选”，不要写成“已证实 EPS 合成基因”。

## 图表建议

1. 主图或补图 Panel A：全蛋白组 CAZyme 大类背景柱图，显示 GT、GH、CBM、CE、AA、PL 的全部候选和高可信候选。
2. Panel B：S0 候选证据矩阵，行是 27 个 S0 基因，列是 dbCAN-HMM、dbCAN-sub、DIAMOND、证据等级、S0 类别、正交家族变化方向。
3. Panel C：S0 与全蛋白组背景的富集图，展示任意高可信 CAZyme 和高可信 GT 的比例、OR 和 P 值。
4. Panel D：S0 GT 家族组成图，突出 GT32 和 GT49 是主构成，GT47/GT10 是较少但高可信的补充证据。

## 对应交付文件

- `FACHB416_paper1_round11_dbcan5_summary_report_zh_20260704.md`
- `FACHB416_paper1_round11_dbcan5_s0_candidate_detail_zh_20260704.tsv`
- `FACHB416_paper1_round11_dbcan5_full_cazyme_background_zh_20260704.tsv`
- `FACHB416_paper1_round11_dbcan5_family_counts_zh_20260704.tsv`
- `FACHB416_paper1_round11_dbcan5_class_counts_zh_20260704.tsv`
