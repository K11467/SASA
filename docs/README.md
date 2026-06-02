# 文档导航

这个目录用于集中放置项目说明、模块说明和整体思路，避免根目录出现大量零散文档。

建议阅读顺序：

1. [pipeline_overview.md](pipeline_overview.md)
   先看项目整体目标、为什么要做 `SASA / ΔSASA / ESM / EGNN` 这条链路。
2. [module_1_sasa.md](module_1_sasa.md)
   了解第一部分的 `SASA` 计算逻辑与示例输出。
3. [module_2_delta_sasa.md](module_2_delta_sasa.md)
   了解第二部分的 `ΔSASA` 标签构造、批量数据集和训练主表。
4. [module_3.md](module_3.md)
   查看 ESM-2 650M、EGNN、cross-chain EGNN 与外部 benchmark 的当前结果。
5. [新增创新.md](新增创新.md)
   查看已经落地的升级和删除的历史产物。
6. [project_spec.md](project_spec.md)
   保留仓库原始说明文件，便于回看最初交接内容。

如果只想快速上手运行代码，优先看根目录的 [README.md](../README.md)。
