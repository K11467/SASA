import sys
from pathlib import Path
import freesasa
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# 确保脚本能找到 src 目录下的 sasa_project 包
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from sasa_project.sasa import parse_pdb, load_dots, calculate_sasa, aggregate_residue_sasa
from sasa_project.paths import EXAMPLE_DATA_DIR

def get_freesasa_results(pdb_filepath):
    """调用 FreeSASA 计算 Total SASA 和 Residue-level SASA (兼容最新版API)"""
    structure = freesasa.Structure(str(pdb_filepath))
    result = freesasa.calc(structure)
    
    total_sasa = result.totalArea()
    
    fs_residues = {}
    # 最新版 API：直接获取所有残基面积的嵌套字典
    areas = result.residueAreas()
    
    for chain_id, residues in areas.items():
        for res_num_str, area_obj in residues.items():
            try:
                # 将 FreeSASA 输出的字符串编号转换为整数，以便与自研代码对齐
                res_num = int(res_num_str)
                fs_residues[(chain_id, res_num)] = area_obj.total
            except ValueError:
                # 忽略一些奇怪的、无法转为数字的插入码残基
                continue
                
    return total_sasa, fs_residues

def plot_correlation(my_areas, fs_areas, correlation_score, output_path):
    """绘制自研与 FreeSASA 的残基级相关性散点图"""
    plt.figure(figsize=(8, 6), dpi=150)
    plt.scatter(fs_areas, my_areas, alpha=0.6, edgecolors='w', s=50, label='Residues')
    
    # 绘制完美拟合参考线 (y=x)
    max_val = max(max(my_areas), max(fs_areas))
    plt.plot([0, max_val], [0, max_val], 'r--', lw=2, label='Ideal y=x')
    
    # 图表美化
    plt.title(f"Custom SASA vs FreeSASA (Pearson R: {correlation_score:.4f})", fontsize=14, fontweight='bold')
    plt.xlabel(r"FreeSASA Area ($\AA^2$)", fontsize=12)
    plt.ylabel(r"Custom SASA Area ($\AA^2$)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"-> 散点对比图已保存至: {output_path}")

def main():
    # 1. 准备路径
    pdb_file = EXAMPLE_DATA_DIR / "2iww_H.pdb"
    dot_file = EXAMPLE_DATA_DIR / "Dot.txt"
    plot_output = EXAMPLE_DATA_DIR / "sasa_comparison_plot.png"
    solvent_radius = 1.4
    
    print(f"=== 开始验证 PDB: {pdb_file.name} ===")
    
    # 2. 运行自研 SASA 算法
    print("[1/3] 运行自研 Shrake-Rupley 算法...")
    atoms = parse_pdb(pdb_file)
    dots = load_dots(dot_file)
    my_total_sasa = calculate_sasa(atoms, dots, solvent_radius)
    my_residue_dict = aggregate_residue_sasa(atoms)
    
    # 3. 运行 FreeSASA 基线算法
    print("[2/3] 运行 FreeSASA 基线算法...")
    fs_total_sasa, fs_residue_dict = get_freesasa_results(pdb_file)
    
    # 4. 对齐残基并提取数据
    print("[3/3] 对齐残基数据并计算评价指标...")
    my_aligned_areas = []
    fs_aligned_areas = []
    
    for key, my_area in my_residue_dict.items():
        # key 格式: (chain_id, residue_id, insertion_code, residue_name)
        chain_id, residue_id = key[0], key[1]
        
        if (chain_id, residue_id) in fs_residue_dict:
            my_aligned_areas.append(my_area)
            fs_aligned_areas.append(fs_residue_dict[(chain_id, residue_id)])
            
    # 5. 计算评价指标
    rel_error = abs(my_total_sasa - fs_total_sasa) / fs_total_sasa * 100
    correlation, _ = pearsonr(my_aligned_areas, fs_aligned_areas)
    
    # 6. 打印报告
    print("\n" + "="*40)
    print("           验证实验结果报告           ")
    print("="*40)
    print(f"自研 Total SASA   : {my_total_sasa:.2f} Å²")
    print(f"FreeSASA Total SASA: {fs_total_sasa:.2f} Å²")
    print(f"Total SASA 相对误差: {rel_error:.2f}%")
    print("-" * 40)
    print(f"有效对齐残基数量   : {len(my_aligned_areas)} 个")
    print(f"残基级 皮尔逊相关系数: {correlation:.4f}  <-- (越接近 1 越好)")
    print("="*40 + "\n")
    
    # 7. 绘制散点回归图
    plot_correlation(my_aligned_areas, fs_aligned_areas, correlation, plot_output)

if __name__ == "__main__":
    main()