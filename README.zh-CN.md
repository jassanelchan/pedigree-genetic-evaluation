# 系谱遗传评估

[English](README.md) | **简体中文**

基于 R、BLUPF90 和 DMU 的畜禽遗传评估流程。从表型和系谱文件开始，完成方差组分估计、育种值预测和年度趋势汇总。

可选择 BLUPF90、DMU，或让两个软件使用相同数据和模型进行对照。每次运行保存输入文件、分析代码、参数卡、编号映射和求解器日志。

## 开始使用

有三种使用方式：导入 skill、命令行、R / RStudio。三种方式调用同一套 R 脚本，完整评估都需要完成下方的[软件配置](#软件配置)。

### 1. 导入 skill

安装 Node.js 和 npm 后，在终端复制运行：

```sh
npx skills add jassanelchan/pedigree-genetic-evaluation --skill pedigree-genetic-evaluation
```

根据提示选择应用和安装位置。安装器会添加技能说明及配套文件，R 和遗传评估软件需要另行安装。

然后提供表型、系谱文件和分析要求，例如：

> 使用 $pedigree-genetic-evaluation，用 BLUPF90 分析达 120 kg 日龄和背膘厚。采用双性状模型，性别和结测年月作为固定效应，加性遗传和共同窝效应作为随机效应；进行一次三倍标准差质控，保留代码、参数卡和日志。

[技能说明](skills/pedigree-genetic-evaluation/SKILL.zh-CN.md)列出了支持的模型和检查内容。导入命令使用 [skills 安装器](https://github.com/vercel-labs/skills)。

### 2. 命令行运行

下载仓库并进入目录：

```sh
git clone https://github.com/jassanelchan/pedigree-genetic-evaluation.git
cd pedigree-genetic-evaluation
```

安装 R 并将 `Rscript` 加入 `PATH`，首次使用时安装依赖包：

```sh
Rscript --vanilla -e "install.packages(c('jsonlite', 'digest', 'processx', 'openxlsx', 'readxl', 'Matrix'), repos='https://cloud.r-project.org')"
```

先检查仓库自带的模拟示例。下面这条命令只准备输入，不启动任何求解器：

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/project.json --engine blupf90 --prepare-only
```

配置好 BLUPF90 后，运行完整评估：

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/project.json --engine blupf90
```

将 `blupf90` 换成 `dmu` 或 `compare`，即可使用 DMU 或两套软件。分析自己的数据时，复制示例 JSON，修改列映射和模型，再将文件路径传给 `--spec`。具体字段见[项目配置说明](skills/pedigree-genetic-evaluation/references/project-spec.zh-CN.md)。

### 3. 在 R / RStudio 中运行

把工作目录设为下载后的仓库。如果尚未安装依赖包，先在 R 控制台运行：

```r
install.packages(
  c("jsonlite", "digest", "processx", "openxlsx", "readxl", "Matrix"),
  repos = "https://cloud.r-project.org"
)
```

然后运行示例：

```r
script <- normalizePath(
  "skills/pedigree-genetic-evaluation/scripts/run_eval.R"
)

project <- normalizePath(
  "skills/pedigree-genetic-evaluation/references/example/project.json"
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

prepare_only <- TRUE

status <- system2(
  rscript,
  c(
    "--vanilla", shQuote(script),
    "--spec", shQuote(project),
    "--engine", "blupf90",
    if (prepare_only) "--prepare-only"
  )
)

stopifnot(status == 0L)
```

配置好求解器后，将 `prepare_only <- TRUE` 改成 `prepare_only <- FALSE` 即可运行完整评估。这段代码会启动独立的 R 进程，不要用 `source()` 加载命令行入口脚本。

## 软件配置

已验证的环境为 **Windows、R 4.5.3**。其他操作系统和程序版本尚未验证。

| 程序 | 已验证版本 | 使用场景 |
| --- | --- | --- |
| RENUMF90 | 1.169 | BLUPF90 或双软件对照 |
| BLUPF90+ | 2.76 | BLUPF90 或双软件对照 |
| DMU | Version 6，release 5.6 | DMU 或双软件对照 |

请另行安装所需程序：[BLUPF90](https://nce.ads.uga.edu/wiki/) 和 [DMU](https://dmu.ghpc.au.dk/dmu/)。仓库不包含可执行文件。DMU 适配器使用包含独立 `dmu1`、`dmuai`、`dmu4` 程序的版本。

在下载后的仓库中，用 PowerShell 复制本机配置模板：

```powershell
Copy-Item skills/pedigree-genetic-evaluation/host_profile.example.json skills/pedigree-genetic-evaluation/host_profile.json
```

用文本编辑器打开 `host_profile.json`，把示例路径改成实际安装路径，保留所选软件需要的条目：

```json
{
  "software": {
    "renumf90": "C:/software/BLUPF90/renumf90.exe",
    "blupf90": "C:/software/BLUPF90/blupf90+.exe",
    "dmu_dir": "C:/software/DMU/bin"
  }
}
```

如果通过安装器导入 skill，就在安装后的技能目录内完成同样的复制和修改。`host_profile.json` 已加入 Git 忽略规则。项目 JSON 中的 `software` 设置优先于本机配置。RENUMF90 和 BLUPF90+ 也可通过 `PATH` 查找。

如果 Windows 会话无法正确读取中文路径或列名，在启动 R 前先在 PowerShell 中设置：

```powershell
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'
```

## 支持的模型

| 模型 | 性状数 | 随机效应 |
| --- | --- | --- |
| 动物模型（`animal`） | 1 个或多个 | 加性遗传效应、残差 |
| 共同窝模型（`litter`） | 1 个或多个 | 加性遗传效应、共同窝效应、残差 |
| 重复力模型（`repeatability`） | 1 个或多个 | 加性遗传效应、永久环境效应、残差 |
| 通用模型（`general`） | 1 个或多个 | 加性遗传效应、多个命名的共同窝、永久环境或其他独立分组效应、残差 |

前三种是常见设计的快捷配置，原有项目文件仍可使用。通用模型还可以在一次分析中合并一次测定性状和重复记录性状。配置接受 `schema_version` 为 `1` 或 `2`。性状数量仍需通过数据、可辨识性和软件适配检查；三性状和四性状设计已在两套已验证版本的软件中实际运行。

固定效应可以作用于全部性状，也可以通过 `fixed[].traits` 指定部分性状，支持分类效应和连续协变量。加性遗传效应始终作用于全部性状；其他随机项各自指定名称、角色（`litter`、`pe` 或 `iid`）及适用性状。同一随机项可以估计适用性状之间的协方差，不同命名随机项之间相互独立。

共同窝效应需要窝号，或母亲编号与仔猪出生日期；永久环境效应按动物编号分组。重复记录需要由个体和真实事件组成的记录键。混合一次与重复记录时，每个性状都要明确指定所属事件：生长测定保留在原来的测定行，每次分娩各占一行，不能把生长表型复制到每一胎。没有成对观测的性状，其残差协方差需要明确约束为零，程序不会自动替用户决定。

目前支持高斯系谱模型。计数性状会检查是否为非负整数，再使用高斯模型拟合。暂不支持基因组、母体遗传、异质残差、阈值、随机回归和未知亲本遗传组模型，也不支持不同命名随机项之间的协方差。

### 生长与繁殖联合示例

[多性状示例](skills/pedigree-genetic-evaluation/references/example/multitrait-project.json)将一次测定的达 120 kg 日龄、背膘厚与多胎总产仔数合并分析。性别用于两个生长性状，胎次和分娩季节用于总产仔数；`birth_litter` 作用于生长性状，`sow_pe` 作用于总产仔数。三个性状通过加性遗传协方差联系；生长与总产仔数来自不同事件，对应残差协方差在示例中明确约束为零。

修改前先复制配置。在仓库根目录打开 PowerShell，运行：

```powershell
Copy-Item skills/pedigree-genetic-evaluation/references/example/multitrait-project.json skills/pedigree-genetic-evaluation/references/example/multitrait-local.json
```

先检查复制后的配置，不启动求解器：

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/multitrait-local.json --engine blupf90 --prepare-only
```

软件配置完成后，去掉 `--prepare-only` 即可运行完整评估。使用 `--engine dmu` 或 `--engine compare`，可改用 DMU 或进行双软件对照。数据路径以复制后的 JSON 所在目录为基准；如果将配置移到别处，需要同步修改路径或使用绝对路径。原来的[两套分析示例](skills/pedigree-genetic-evaluation/references/example/project.json)仍可用于分别分析生长与繁殖性状。

## 输入与结果

项目 JSON 指定文件、列名、性状、单位、模型和质控规则。支持 Excel、CSV、TSV 和分隔符文本。[仓库示例](skills/pedigree-genetic-evaluation/references/example/project.json)全部使用模拟数据，来源见[示例说明](skills/pedigree-genetic-evaluation/references/example/provenance.zh-CN.md)。

默认按每个性状独立进行一次均值 ±3 SD 质控。一个性状缺失或未通过质控，不会连带删除同一个体的另一个有效性状。父母信息可从其他文件补充；关系冲突和重复事件会报告出来，供核实。

| 输出 | 内容 |
| --- | --- |
| `delivery/evaluation_summary.xlsx` | 质控、描述性统计与变异系数、方差组分、遗传力、适用时的重复力或遗传相关、标准误及年度汇总 |
| `delivery/analysis_*/EBV.csv` | 原始个体号、EBV、SEP，以及条件允许时的可靠性 |
| `delivery/analysis_*/*.png` | 质控直方图和年度趋势图 |
| 运行目录 | 输入副本、复现配置、代码、参数卡、编号映射、求解器日志、版本及运行清单 |

每个命名随机项的方差组分分别输出。某个性状的表型方差包括所有作用于它的方差项；只有该性状包含永久环境效应时才输出重复力，分子为加性遗传方差与永久环境方差之和，共同窝及其他独立分组效应不计入分子。

遗传参数的标准误使用完整采样协方差矩阵与 Delta 方法。年度汇总默认按测定或分娩年份进行，同一动物的每个性状在同一年只计一次 EBV。混合记录设计只使用属于该性状的事件行，因此母猪的分娩年份不会被当作生长性状的测定年份。结果含义和数值核验见[结果说明](skills/pedigree-genetic-evaluation/references/results.zh-CN.md)。

输入路径和相对形式的 `--output` 都以**项目 JSON 所在目录**为基准。自带示例因此会写入 `skills/pedigree-genetic-evaluation/references/example/example_runs/`，每次运行各建一个子目录。

继续已有运行时，添加 `--resume <运行目录>`。数据、模型和程序的指纹一致时，会复用已核验结果。只修改年度汇总或图像，不需要重新估计方差组分。

## 详细文档

- [项目配置](skills/pedigree-genetic-evaluation/references/project-spec.zh-CN.md)
- [BLUPF90 适配](skills/pedigree-genetic-evaluation/references/blupf90.zh-CN.md)
- [DMU 适配](skills/pedigree-genetic-evaluation/references/dmu.zh-CN.md)
- [结果与标准误](skills/pedigree-genetic-evaluation/references/results.zh-CN.md)
- [技能说明](skills/pedigree-genetic-evaluation/SKILL.zh-CN.md)

每份说明都有英文和简体中文版本，可通过页首链接切换。

## 测试

数据处理和统计测试不需要安装求解器：

```sh
Rscript --vanilla tests/run_tests.R
```

测试和仓库示例均使用模拟数据。离线测试覆盖数据准备与统计计算，包括性状专属效应、一次与重复记录混合，以及多个随机项。原生 REML、BLUP 实测覆盖原有模型、三性状与四性状设计、部分性状缺失、近交，以及上述已验证版本的跨软件对照。

运行原生软件回归检查前，在本机配置中填好两套软件路径，再运行：

```sh
Rscript --vanilla tests/run_native_tests.R --profile skills/pedigree-genetic-evaluation/host_profile.json
```

该命令会对两份自带项目配置运行双软件对照，并检查已完成运行在续跑时不会重新启动求解器。添加 `--output <绝对目录>` 可指定测试结果保存位置。
