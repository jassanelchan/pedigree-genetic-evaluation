[English](project-spec.md) | [简体中文](project-spec.zh-CN.md)

# 项目配置与运行

JSON 配置记录表型来源、模型和软件选择。数据路径可以是绝对路径，也可以相对于配置文件所在目录；程序路径应使用绝对路径。Windows 路径建议使用正斜杠。命令行在全新的 R 会话中运行，不依赖已有 R 工作区。

## 项目字段

| 字段 | 含义 |
| --- | --- |
| `schema_version` | 整数 `1` |
| `engine` | 必须明确指定 `blupf90`、`dmu` 或 `compare` |
| `output_dir` | 输出根目录，每次新运行使用独立子目录 |
| `software` | 覆盖本机配置中的 `renumf90`、`blupf90`、`dmu_dir` |
| `timeout` | 单个进程的最长秒数，默认 7200 |
| `threads` | 请求的线程数，默认 2，不超过本机容量 |
| `pedigree` | 共用系谱来源，可被单项分析的设置覆盖 |
| `analyses` | 非空分析列表，各项名称唯一、分别准备模型 |

`dmu_dir` 指向包含 `dmu1.exe`、`dmuai.exe` 和 `dmu4.exe` 的目录。R5.6 还需要发行包自带的共享运行库，仓库不分发这些程序。下载前请阅读[后端说明](dmu.zh-CN.md)。依赖的 R 包为 `jsonlite`、`digest`、`processx`、`openxlsx`、`Matrix`，以及读取 Excel 所需的 `readxl`。

## 数据表来源

```json
{
  "path": "phenotypes.xlsx",
  "sheet": "Sheet1",
  "columns": {"id": "个体号", "sire": "父亲", "dam": "母亲"},
  "missing_values": [-999],
  "date_origin": "1899-12-30"
}
```

- 表型来源只要求 `columns.id`。可选的父号、母号字段用于补充不在系谱文件中的动物的真实亲本信息，补充前会检查是否矛盾。
- 系谱来源要求 `columns.id`、`columns.sire` 和 `columns.dam`。
- 可选的 `pedigree.supplements` 是额外数据表来源列表，必须带有经核实的个体/父号/母号映射，例如生长导出表中的真实父母字段。这些来源补充给所有使用该系谱的分析；已知亲本矛盾时停止。共用群体的补充来源需明确声明，避免某模型仍把有已知亲本的祖先当成创始动物。
- 列可以用准确列名或从 1 起算的整数位置指定，Excel 优先使用准确列名。
- 多工作表文件需要指定工作表。所有字段按文本读取，保留字母 ID 和前导零；若 Excel 已经丢失前导零，不能在没有其他来源的情况下恢复。
- 文本文件支持 `header`（默认 true）、`sep`（CSV 为逗号，TSV 为制表符，其余默认空白）和 `encoding`（默认 UTF-8-BOM）。无表头的原生输入文件使用 `header:false`。
- 空值、NA、N/A、NULL 视为缺失。数值性状还识别 NAN 及明确给出的 `missing_values`。亲本 ID 中的 `0`、`-`、`.` 表示未知。
- 未指定 `date_origin` 时，数值 XLSX 日期使用经核实的工作簿日期系统。CSV/XLS 中的数值日期需明确给出已核实的起始日期。文本日期默认接受 YYYY-MM-DD 或 YYYY/MM/DD，也可用 R 日期格式指令指定 `date_format`。

## 单项分析字段

| 字段 | 含义 |
| --- | --- |
| `name`、`species` | 必填的分析名称与已知物种 |
| `phenotypes` | 表型数据表来源 |
| `model` | `animal`、`litter` 或 `repeatability` |
| `traits` | 有顺序的一个或两个性状定义 |
| `fixed` | 固定分类效应或协变量列表，程序另加显式均值 |
| `shared` | 共同窝模型使用已核实窝号的 `column`，或 `dam` 与 `birth_date` |
| `event_date`、`birth_date` | 来源中的列名/列号，只有被效应或记录键使用时才必须提供 |
| `record_key` | 唯一标识一条记录的来源列，重复记录模型必填 |
| `qc` | `method:"sd"`、`sd:3`，可选分组来源列 `group_by`；或明确指定 `method:"none"` |
| `inbreeding` | 布尔值，默认 true |
| `predict_ids` | 可选的无表型预测对象原始 ID，需在提供的系谱中有表示 |
| `allow_observed_founders` | 是否明确允许将系谱中完全没有表示的观测/目标 ID 加为互不相关的创始动物，默认 false |
| `trend_basis` | 默认 `event_year`，也可明确使用 `birth_year` |
| `initial` | 可选的 additive/shared/residual 方阵嵌套数组，性状顺序与定义一致 |
| `zero_covariances` | 明确的约束 `{component,trait1,trait2}`；性状索引从 1 起算，component 为 additive/shared/residual |

每个性状需提供 `name`、`column`、`unit`，可选 `label` 与 `distribution`。支持 `gaussian` 和 `count_gaussian`；后者先验证非负整数，再拟合高斯模型，不是 Poisson 模型。性状和单位必须有用户或数据依据，流程不会自动编造体重校正值。

固定效应使用 `name`、`column` 和 `type`，其中 type 默认为 `class`，连续协变量为 `cov`。由日期生成分类效应时，改用 `from_date:"event"` 或 `"birth"`，并指定 `grouping:"year"`、`"month"` 或 `"quarter"`。两个性状使用相同固定效应。每个性状分别统计固定效应水平的观测覆盖，包括某性状完全没有观测的水平。秩检查采用参照编码的可估计设计矩阵；出现空单元格本身不构成删除整年的理由。

记录键应包含动物和事件。键中使用的日期在查重前统一格式，日期写法不同不能让同一分娩变成两条记录。重复产仔数记录归属于母猪，不使用其仔猪的个体 ID。母猪永久环境效应与母猪自己的出生窝效应不同。

质控在声明的分组内对有限观测值执行一遍。筛除前计算阈值，标准差使用 n−1 分母。超界性状改为缺失，只要任一性状仍有效就保留该记录。完全无表型的行不要求固定效应完整，但所提供的亲本信息仍保留。筛除后重新编码分类和共同窝水平。有效记录缺少必要固定效应或窝号时，报告具体错误，不自动填补。

默认的自由协方差初值对应相关系数 0.1，组分比例为加性 0.3、共享效应 0.1（如有）、残差 0.6/0.7。这些只是初值，不是最终估计。零协方差属于模型约束，必须单独声明。v1 仅支持非对角零约束，固定方差不属于通用接口。

## 运行命令

```powershell
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'
& '<Rscript-path>' --vanilla '<skill-path>\scripts\run_eval.R' --spec '<project-path>\project.json'
```

使用本机实际的 Rscript 路径。`--engine` 和 `--output` 可覆盖配置中的对应字段。`--prepare-only` 只准备文件，不运行求解器。`--resume <prior-run-directory>` 仅在指纹一致时读取经验证的缓存。每次调用都把脚本和当前配置复制到新的调用目录。数据/模型指纹不包含仅用于趋势图的年份分组，因此更换作图年份口径不会重新执行 REML。

## 输出与失败状态

`delivery/evaluation_summary.xlsx` 汇总模型、描述性统计/CV、质控、固定效应水平覆盖、输入审计、方差组分/遗传参数、年度均值及斜率。各软件子目录包含原始 ID 的 EBV CSV、详细工作簿和图片。对照目录包含方差组分、EBV/PEV，以及原生输出支持时的拟合值比较。

`analysis_<n>/` 保留准备后的模型上下文、ID 映射、秩检查、准确模型定义和来源哈希。后端目录保留生成的参数卡/数据、REML 尝试、日志、进程状态，以及经验证的拟合/EBV 缓存。`code_records/` 保留完整可读的源代码快照。`run_manifest.json` 记录软件、R 和包的来源信息，状态为 `prepared`、`complete` 或 `partial`。退出码 2 表示需要检查错误和部分结果，不能据此宣布评估全部成功。

`input_sources/` 保留全部原始输入和亲本补充表的 SHA256 校验副本，相同输入按哈希去重。`input_manifest.csv` 记录原路径与快照路径。`reproduce_project.json` 使用相对快照路径，使整个运行目录转移后不依赖原始桌面目录。使用保存的 `code_records/invocation_.../run_eval.R` 和该配置文件复现，并在目标机器配置有效程序路径与 R 包。流程不覆盖源文件。原始快照含有群体私有数据，应保留在项目复现目录中，不能放入可共享的技能包。
