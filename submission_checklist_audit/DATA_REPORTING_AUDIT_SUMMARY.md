# Data Reporting Checklist 审计概览

## 1. 当前准备程度

本次审计基于仓库真实存在的文件完成，只做了扫描、统计、核对和整理，没有修改正文、主图、补充材料或核心数据。当前仓库已经能为 Data Reporting Checklist 提供相当一部分证据，尤其是研究设计、数据来源边界、外部 sensory panel 聚合数据、M2OR 匹配结果、benchmark 指标表和基础代码可用性说明。

按 evidence matrix 统计，`READY` 项 3 个，`NEEDS AUTHOR CONFIRMATION` 项 11 个，`MISSING` 项 1 个。整体状态不是“空缺很多”，而是“已有较强仓库证据，但仍缺若干投稿系统级文字声明和少量关键归档证据”。

## 2. 已经可以直接支持 checklist 的内容

- 研究范围可由 `README.md` 和 `docs/workflow_overview.md` 直接支持：仓库明确覆盖 benchmark、多种分子表征、structure-disjoint validation、external sensory-panel validation、receptor-consistency analysis，并明确说明这不是 receptor mechanism proof。
- 第三方数据来源与限制已有直接证据：`data/raw_manifest/source_inventory.csv`、`data/raw_manifest/third_party_source_restrictions.md`、`docs/data_availability.md`。
- 外部 sensory panel 的聚合层面说明较完整：`docs/sensory_panel_metadata.md`、`data/external_panel/aggregated_panel_ratings.csv`、`panel_derived_topk_labels.csv`、`pairwise_panel_jaccard_similarity.csv`。
- receptor annotation reporting 证据较完整：`data/receptor_annotations/` 下的 matched、vectors、pairwise concordance、concordant/divergent、exclusion log 均存在。
- benchmark 结果指标表存在：`results/benchmark_metrics/` 下 micro-AUC、macro-AUC、labelwise delta AUC 已找到。

## 3. GitHub 中虽有证据，但投稿系统仍需要作者补文字说明的内容

- 数据集规模目前存在“请求值”与“仓库实际值”不一致问题。仓库实际反复写的是 `4,495` benchmark molecules、`87` external molecules，而不是你在任务中列出的 `4,403` 和 `90`。这类内容必须由作者确认最终应在 checklist 中写哪一组数字。
- `docs/workflow_overview.md` 说明了 training-only scaling、iterative multi-label stratification、Tanimoto > 0.8 的 structure-disjoint 逻辑，但仓库没有完整公开 `iterative_stratified_splits/` 目录，也没有 `feature_scaler_parameters.csv`，因此投稿系统里仍需要作者补解释。
- `configs/model_hyperparameters.yaml` 和 `configs/graph_model_settings.yaml` 证明模型家族、随机种子和部分训练设定存在，但大量超参数仍是 `TODO`，不能直接当作最终超参数表述。
- `.zenodo.json` 已存在，说明仓库为归档做过准备，但当前未看到 DOI 或 release tag 证据，因此 checklist 里的 data/code availability 仍需作者补全。

## 4. 当前缺失且必须由作者确认或补充的内容

- 最优先缺失项是 sensory panel 相关的 ethics approval / waiver / informed consent 明确声明。仓库里只有匿名化和聚合发布说明，没有找到明确伦理或同意文本。
- Zenodo DOI 或正式 release tag 目前未见证据。
- 迭代分层 split 文件、feature scaler 参数文件、以及更强的 leakage-prevention 归档证据不完整。
- `supplementary_tables_source_data.xlsx` 未找到；`figures/` 目录和 `make_fig1.py` 到 `make_fig4.py` 也未找到，因此图和 source data 的可追溯性还需要人工补说明。
- QM 相关公开目录 `qm_calculations/` 不存在，所以 Gaussian / Multiwfn 原始或中间输出类证据不能从当前公共仓库直接支持，只能说明“当前 public snapshot 未包含”。

## 5. 是否建议现在上传 Data Reporting Checklist

不建议直接按当前状态原样上传最终版 checklist，但建议基于本次生成的草稿快速补齐缺失项后提交。原因不是仓库质量差，而是投稿系统通常对伦理、归档 DOI、第三方限制、split leakage 证据和最终超参数描述更敏感，这些地方目前还需要作者级确认。

## 6. 最终判断

当前判断：**需要补充后提交**。

补充说明：

- 轻量测试状态：TESTS NOT RUN: dependency unavailable
- leakage 检查摘要：cas overlap=0; smiles overlap=0
- 第三方限制说明文件：已找到
- ethics / informed consent 明确声明：未找到
- Zenodo DOI：未看到
