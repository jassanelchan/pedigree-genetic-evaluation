[English](project-spec.md) | [简体中文](project-spec.zh-CN.md)

# 项目配置与运行

JSON 配置记录表型来源、模型和软件选择。数据路径可以是绝对路径，也可以相对于配置文件所在目录；程序路径应使用绝对路径。Windows 路径建议使用正斜杠。命令行在全新的 R 会话中运行，不依赖已有 R 工作区。

## 项目字段

| 字段 | 含义 |
| --- | --- |
| `schema_version` | 整数 `1` 或 `2`，原有版本 1 项目仍然有效 |
| `engine` | 必须明确指定 `blupf90`、`dmu` 或 `compare` |
| `output_dir` | 输出根目录，每次新运行使用独立子目录 |
| `software` | 覆盖本机配置中的 `renumf90`、`blupf90`、`dmu_dir` |
| `timeout` | 每个 F90 进程，或 DMU1 与后续模块合计的最长秒数，默认 7200 |
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
| `model` | `animal`、`litter`、`repeatability` 或 `general` |
| `traits` | 有顺序的非空性状列表，不再固定限制为两个性状 |
| `fixed` | 固定分类效应或协变量，可指定适用的性状名称；程序另加显式均值 |
| `shared` | 仅用于共同窝快捷模型，使用已核实窝号的 `column`，或 `dam` 与 `birth_date` |
| `random` | 仅用于通用模型，声明独立命名的 `litter`、`pe` 或 `iid` 项 |
| `event_date`、`birth_date` | 来源中的列名/列号，只有被效应或记录键使用时才必须提供 |
| `record_key` | 唯一标识一条记录的来源列，重复记录模型必填 |
| `qc` | `method:"sd"`、`sd:3`，可选分组来源列 `group_by`；或明确指定 `method:"none"` |
| `inbreeding` | 布尔值，默认 true |
| `predict_ids` | 可选的无表型预测对象原始 ID，需在提供的系谱中有表示 |
| `allow_observed_founders` | 是否明确允许将系谱中完全没有表示的观测/目标 ID 加为互不相关的创始动物，默认 false |
| `trend_basis` | 默认 `event_year`，也可明确使用 `birth_year` |
| `initial` | 可选的完整性状×性状嵌套数组，以 additive、随机项名称和 residual 为键 |
| `zero_covariances` | 明确的约束 `{component,trait1,trait2}`，指定组分键及从 1 起算的全局性状索引 |

每个性状需提供 `name`、`column`、`unit`，可选 `label`、`distribution`、`record_type` 和 `record_filter`。支持 `gaussian` 和 `count_gaussian`；后者先验证非负整数，再拟合高斯模型，不是 Poisson 模型。性状和单位必须有用户或数据依据，流程不会自动编造体重校正值。

固定效应使用 `name`、`column` 和 `type`，其中 type 默认为 `class`，连续协变量为 `cov`。由日期生成分类效应时，改用 `from_date:"event"` 或 `"birth"`，并指定 `grouping:"year"`、`"month"` 或 `"quarter"`。可选的 `traits:["age120","bf120"]` 将该效应限定到这些性状，不填则适用于全部性状。例如，性别可只用于生长性状，胎次只用于产仔数。每个适用性状分别检查固定效应覆盖与秩。只有行内存在适用性状的有效观测时，才要求该固定效应有值；不适用行在内部使用中性占位值，不代表填补了观测。秩检查采用参照编码的可估计设计矩阵；出现空单元格本身不构成删除整年的理由。

## 快捷模型与命名随机项

所有模型都包含作用于全部性状的一个系谱加性遗传项，以及一个残差协方差矩阵。原有快捷模型均可用于一个或多个性状：

- `animal`：不添加其他随机项。
- `litter`：增加名为 `shared`、按单位矩阵建模的共同窝项，由 `shared` 定义分组。
- `repeatability`：增加名为 `shared` 的动物永久环境项，性状默认使用重复记录。

需要在不同性状子集上组合多个独立随机项时，使用 `general`，不再填写 `shared`。例如：

```json
"random": [
  {"name": "birth_litter", "role": "litter", "column": "litter", "traits": ["age120", "bf120"]},
  {"name": "sow_pe", "role": "pe", "traits": ["TNB"]}
]
```

加性项由程序自动加入，不要再写入 `random`。每个额外随机项需要唯一名称，符合 `[A-Za-z][A-Za-z0-9_]*`；`additive`、`residual` 和 `shared` 为保留名称。`traits` 可省略，默认适用于全部性状。role 决定分组和解释：

| role | 分组 | 含义 |
| --- | --- | --- |
| `litter` | 已核实的分组列 `column`，或 `dam` 与分析中的 `birth_date` | 共同出生窝效应 |
| `pe` | 动物 ID，不允许另设 `column` 或 `dam` | 同一动物多次观测共享的永久环境效应 |
| `iid` | 必填分组列 `column` | 其他按单位矩阵建模的分组随机效应 |

每个随机项可在自身适用的性状间估计协方差，不同命名随机项为相互独立的块。共同窝和一般 `iid` 方差不能因为属于环境效应就计入重复力分子。每个单位矩阵随机项在每个适用性状上至少需要两个有观测的组，且至少一组有重复；同一性状上具有完全相同分组划分的两个随机项会被拒绝。增加随机项仍需科学依据和可识别的设计。

当前范围为高斯系谱模型，不支持基因组关系、母体遗传效应、未知亲本遗传组、异质残差、阈值模型和随机回归。性状与协方差参数越多，对数据和计算量的要求也越高，并非模型越大越好。

## 单次记录与重复记录联合分析

`record_type` 默认 `single`，也可设为 `repeated`；重复力快捷模型默认 `repeated`。单次记录性状每头动物最多有一个有效观测。只要有重复记录性状，就必须提供含动物 ID 与真实事件标识或日期的 `record_key`。不要把同一条生长表型复制到每次分娩行来制造重复记录。

通用模型混合单次与重复记录性状时，每个性状都必须声明所属事件。例如：

```json
{"name": "age120", "column": "age120", "unit": "days", "record_type": "single", "record_filter": {"column": "event_type", "values": ["growth"]}}
```

重复 TNB 性状使用对应的 `farrowing` 筛选。筛选定义哪些行属于该性状，也包括属于该性状但没有测量值的行。筛选范围外出现非缺失性状值时会报错，不会静默丢弃。某性状不属于该行时保持缺失，不适用的固定效应和随机分组字段可留空。

[多性状示例](example/multitrait-project.json)将两个单次记录生长性状与重复 TNB 联合分析。生长和分娩分别作为真实事件行；性别用于生长，胎次用于 TNB，共同出生窝用于生长，母猪 PE 用于 TNB。动物和系谱通过遗传协方差连接这些性状。没有任何行同时观测生长性状与 TNB，因此在 `zero_covariances` 中明确将它们的残差协方差固定为 0。

记录键应包含动物和事件。键中使用的日期在查重前统一格式，日期写法不同不能让同一分娩变成两条记录。重复产仔数记录归属于母猪，不使用其仔猪的个体 ID。母猪永久环境效应与母猪自己的出生窝效应不同。

质控在声明的分组内对有限观测值执行一遍。筛除前计算阈值，标准差使用 n−1 分母。超界性状改为缺失，只要任一性状仍有效就保留该记录。完全无表型的行不要求固定效应完整，但所提供的亲本信息仍保留。筛除后重新编码分类和随机分组。适用且有观测的行缺少必要固定效应或随机分组时，报告具体错误，不自动填补。

## 协方差矩阵与约束

`initial` 中每个组分都使用完整、有顺序的性状列表，即使某随机项只适用于部分性状。键为 `additive`、`residual` 及声明的随机项名称，快捷模型保留 `shared` 名称。例如性状顺序为 `[age120,bf120,TNB]` 时，仅作用于前两个性状的共同窝项使用 3×3 矩阵，最后一行和一列为 0。适用性状对应的主子矩阵必须正定，不适用的元素是结构性零，不参与估计。

用户约束中的 `trait1`、`trait2` 必须为整数，始终对应从 1 起算的全局性状顺序，不是后端的局部编号。例如：

```json
"zero_covariances": [
  {"component": "residual", "trait1": 1, "trait2": 3},
  {"component": "residual", "trait1": 2, "trait2": 3}
]
```

保留的真实记录中没有任何一行同时观测某对性状时，它们的残差协方差不可估，必须明确声明有依据的零约束。单位矩阵随机项中，即使跨记录，某对适用性状也没有任何共同的有观测分组时，同样必须声明零协方差。流程不会自动替用户决定这些假设。若动物与系谱有足够连接，遗传协方差仍可能估计，但取决于数据支持；残差协方差固定为 0 不等于遗传协方差固定为 0。

默认初值将每个性状样本方差的 0.3 分给加性遗传。存在额外随机项时，这些项平分总比例 0.1，残差为 0.6；不存在时，残差为 0.7。某块适用 k 个性状时，自由相关初值为 `0.1/max(1,k−1)`。这些是数值初值，不是最终估计。自由协方差初值必须非零；固定零协方差需单独声明。只支持非对角零约束，固定适用性状的方差不属于通用接口。

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
