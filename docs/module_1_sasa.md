# SASA 工具升级说明

本模块基于 Shrake-Rupley 打点法实现 SASA 计算。

输入：
- `data/raw/examples/2iww_H.pdb` 或其他 PDB 文件
- `data/raw/examples/Dot.txt` 球面采样点文件
- 溶剂探针半径，默认 `1.4`

输出：
- `data/processed/examples/atom_sasa.csv`：原子级 SASA
- `data/processed/examples/residue_sasa.csv`：残基级 SASA
- `data/processed/examples/chain_sasa.csv`：链级 SASA

当前核心代码在 `src/sasa_project/sasa.py` 中，主要包括：
- `Atom`
- `parse_pdb()`
- `load_dots()`
- `calculate_sasa()`
- `aggregate_residue_sasa()`
- `aggregate_chain_sasa()`
- `filter_atoms_by_chain()`
- `write_atom_sasa_csv()`
- `write_residue_sasa_csv()`
- `write_chain_sasa_csv()`

运行方式：

```bash
PYTHONPATH=src python3 scripts/run_sasa_example.py
```

校验项：
- `total_sasa == sum(atom.sasa for atom in atoms)`
- `total_sasa == sum(residue_sasa.values())`
