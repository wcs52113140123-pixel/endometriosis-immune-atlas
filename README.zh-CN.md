# 单细胞研究中的统计单位：子宫内膜异位症免疫微环境

[English](README.md)

对八个公开的子宫内膜异位症单细胞、空间与 bulk 数据集所做的病人层面再分析的代码与派生结果，
并包含对同一统计问题的跨疾病校准研究。

作者：X. Wang, K. Liang, D. Yang, B. Xiong, J. Liao, Z. Wang —
遵义医科大学附属医院妇产科，贵州遵义。

---

## 要问的问题

子宫内膜异位症的病例—对照单细胞研究，几乎都把一组内所有供者的细胞合并，再比较合并后的细胞层面分布。
但同一供者的细胞并非彼此独立的观测。本工作要问的是：把**病人**而不是**细胞**作为统计单位时，结论会有什么变化。

三个结果：

1. **Treg 富集在合并细胞时可复现，按病人统计则消失。** 病人层面各队列 *P* = 0.20、0.63、0.46；
   全部 81 个样本、以队列为分层时，对数尺度 β = 0.064（95% CI −0.56 ~ 0.69，*P* = 0.84）。
   每样本的 M2、CD8 耗竭、NK 细胞毒性评分同样为阴性。
2. **有一项发现通过了病人层面检验。** 一个九基因的 M2／免疫抑制程序，在病灶巨噬细胞分割段中高于
   同一位女性的在位内膜，5 位配对病人全部升高（配对 *P* = 0.031），且仅见于巨噬细胞分割段
   （基质 *P* = 0.41，上皮 *P* = 0.97）。另外两个检验相邻对比的队列不支持该发现，
   因此该结论限于腹膜病灶的巨噬细胞分割段。
3. **这种校准失效是普遍的。** 在 124 种疾病、523 个公开病例—对照情形（6,111 次细胞类型比较）中，
   合并细胞检验对真实零假设的拒绝率中位数为 0.84，而供者层面检验为 0.04；
   合并检验显著的结果中有 69% 在供者层面不成立。

---

## 目录结构

```
.
├── scripts/                     分析与作图代码
├── results/
│   ├── figures/                 已渲染的图
│   └── tables/                  上述每个数字对应的源表
├── data/
│   ├── README.md                数据清单：登记号、用途、获取方式
│   ├── interim/                 派生的中间结果
│   └── bulk/                    bulk 表达与生存矩阵
├── UNIT_OF_REPLICATION.md       可复用诊断工具的说明
├── environment.yml              conda 环境（python 3.11）
└── requirements.txt             pip 等价依赖
```

## 如何运行

```bash
conda env create -f environment.yml && conda activate endo
python scripts/00_download_data.py           # 获取公开矩阵到 data/raw
python scripts/01_load_qc_integrate.py       # 质控、Harmony 整合、谱系注释
python scripts/02_epithelial_infercnv.py     # 上皮子集、推断拷贝数
python scripts/05_statistics.py --part all   # 病人层面检验、功效、方法比较
python scripts/build_fig2_print.py           # 出图，一张图一个脚本
```

`scripts/run_all.py --list` 会列出已脚本化的步骤。`scripts/figure_layout.py` 保存共用的图件版式流程，
并记录每张图所用的参数。

## 那个诊断工具

`scripts/unit_of_replication.py` 可把同一套检验用到任何数据集上。给它一张供者 × 细胞类型的表和一个分组标签，
它会同时跑合并细胞检验与供者层面检验，通过在供者层面置换分组标签来校准，并报告哪些"合并显著"的结果能够存活：

```python
from unit_of_replication import diagnose, risk_note
res = diagnose(df, donor="donor_id", group="group", cell_type="cell_type", case="case")
print(risk_note(res))
```

`UNIT_OF_REPLICATION.md` 说明了各个参数与小样本组所用的确切 Mann–Whitney 零分布。

## 数据

所有原始数据均为公开数据。`data/README.md` 列出每个登记号在分析中的用途与获取命令；
`scripts/00_download_data.py` 会从原始数据库获取。

## 许可与引用

代码以 MIT 许可发布（见 `LICENSE`）。第三方数据集仍受其原始提交者设定的条款约束。
若使用本代码，请按 `CITATION.cff` 引用。

## 字体

图件使用 Arial、11/10/9 pt 三级字号。Arial 需单独授权，此处不包含；
请把 `ARIAL.TTF`、`ARIALBD.TTF`、`ARIALI.TTF`、`ARIALBI.TTF` 放入 `assets/fonts/` 以完全复现排版。
