# 系谱遗传评估

基于 R、BLUPF90 和 DMU 的畜禽系谱遗传评估流程，可作为 Codex skill 使用，也可以从命令行运行。

从表型和系谱文件开始，完成质控、方差组分估计、育种值预测和年度趋势汇总。每次运行会保存输入副本、代码、参数卡和终端日志，方便复核和重跑。

## 支持的模型

| 模型 | 性状数 | 随机效应 |
| --- | --- | --- |
| 动物模型 | 1 或 2 | 加性遗传效应、残差 |
| 共同窝模型 | 1 或 2 | 加性遗传效应、共同窝效应、残差 |
| 重复力模型 | 1 | 加性遗传效应、永久环境效应、残差 |

两个软件使用相同的数据和模型定义。可选择 `blupf90`、`dmu`，或 `compare` 对照运行。

当前版本不包含基因组模型、三性状模型、阈值模型、随机回归和母体遗传效应。计数性状使用高斯模型，输入会先检查是否为非负整数。

## 环境

已验证环境：Windows、R 4.5.3、RENUMF90 1.169、BLUPF90+ 2.76、DMU version 6 release 5.6。

需要以下 R 包：

```r
install.packages(c(
  "jsonlite",
  "digest",
  "processx",
  "openxlsx",
  "readxl",
  "Matrix"
))
```

BLUPF90 和 DMU 程序需要自行安装，本仓库不包含可执行文件。其他程序版本需先核对输入和输出格式。

- [BLUPF90](https://nce.ads.uga.edu/wiki/)
- [DMU](https://dmu.ghpc.au.dk/dmu/)

第三方程序的使用条件以各自官方说明为准。

## 使用

将 `skills/pedigree-genetic-evaluation` 安装到 Codex 的技能目录，然后提供表型、系谱和分析要求：

> 使用 $pedigree-genetic-evaluation，分析这些表型和系谱文件。用 BLUPF90，按一次三倍标准差质控，输出遗传参数、育种值和年度趋势，保留代码与日志。

也可以直接运行 R 脚本。匿名示例在 `references/example` 中：

```powershell
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'

Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R `
  --spec skills/pedigree-genetic-evaluation/references/example/project.json `
  --engine blupf90 `
  --output ./evaluation_runs
```

运行前，在项目 JSON 的 `software` 字段填写软件的绝对路径。也可以复制技能目录中的 `host_profile.example.json` 为 `host_profile.json`，填写本机配置；该文件已加入忽略规则。

首次配置可加 `--prepare-only`，只检查数据和模型，不调用求解器。续跑可使用 `--resume <运行目录>`，相同数据、模型和程序会复用已经核验的结果。

## 数据处理和结果

默认对每个性状独立做一次均值 ±3 SD 质控。某个性状缺失或不通过质控，不会连带删除另一个有效性状。共同窝需要明确的窝号，或母亲编号与仔猪出生日期；重复记录需要个体与事件组成的记录键。

父母信息可以由多个文件补充。重复耳号、冲突的父母关系、重复分娩记录和无法识别的日期会报出具体问题，不会自动猜测修补。

结果包含：

- 质控前后描述性统计和变异系数；
- 方差组分、表型方差、狭义遗传力、重复力、遗传相关及标准误；
- 原始耳号对应的 EBV、SEP，以及条件允许时的可靠性；
- 年度表型均值、平均育种值、趋势拟合和图像；
- 原始输入副本、完整代码、参数卡、编号映射和运行日志。

遗传参数的标准误使用完整采样协方差和 Delta 方法。年度 EBV 汇总时，同一动物在同一年只计一次，避免多胎母猪被重复加权。

F90 启用近交校正时，使用经核验的近交系数和 Henderson 稀疏更新生成 A-inverse，减少原生整数系谱编码的舍入影响。近交系数本身仍有有限输出精度，具体方法和精度信息会随结果保存。

## 测试

不需要安装两套评估软件即可运行公共数据处理和统计测试：

```powershell
Rscript --vanilla tests/run_tests.R
```

测试数据全部由脚本生成，仓库中的示例也均为模拟数据。软件适配器另经实际 REML、BLUP、部分性状缺失、近交、共同窝、永久环境和跨软件对照验证。

## 文件说明

`skills/pedigree-genetic-evaluation/SKILL.md` 是技能入口。参数字段、软件适配和结果定义分别见：

- [项目配置](skills/pedigree-genetic-evaluation/references/project-spec.md)
- [BLUPF90 适配](skills/pedigree-genetic-evaluation/references/blupf90.md)
- [DMU 适配](skills/pedigree-genetic-evaluation/references/dmu.md)
- [结果与标准误](skills/pedigree-genetic-evaluation/references/results.md)
