import sys
import time
from pathlib import Path
import freesasa
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# 确保能导入 src 目录下的包
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from sasa_project.sasa import parse_pdb, load_dots, calculate_sasa, aggregate_residue_sasa
from sasa_project.paths import PROJECT_ROOT, EXAMPLE_DATA_DIR

def get_freesasa_results(pdb_filepath):
    """调用 FreeSASA 获取总面积和残基级面积"""
    structure = freesasa.Structure(str(pdb_filepath))
    result = freesasa.calc(structure)
    total_sasa = result.totalArea()
    
    fs_residues = {}
    areas = result.residueAreas()
    for chain_id, residues in areas.items():
        for res_num_str, area_obj in residues.items():
            try:
                fs_residues[(chain_id, int(res_num_str))] = area_obj.total
            except ValueError:
                continue
    return total_sasa, fs_residues

def plot_batch_results(complex_names, correlations, errors, output_path):
    """绘制极具学术质感的批量验证双子图"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), dpi=150)
    
    x = np.arange(len(complex_names))
    
    # --- 图 1: 皮尔逊相关系数 (证明残基级趋势一致性) ---
    bars1 = ax1.bar(x, correlations, color='#2563eb', alpha=0.8, width=0.6)
    # 动态设置 Y 轴下限，让近乎完美的 0.99 更加明显
    y_min = min(correlations) - 0.01 if min(correlations) > 0.9 else 0.8
    ax1.set_ylim(y_min, 1.005)
    ax1.set_xticks(x)
    ax1.set_xticklabels(complex_names, rotation=45, ha='right', fontsize=9)
    ax1.set_ylabel('Pearson R', fontsize=11, fontweight='bold')
    ax1.set_title('Residue-Level SASA Correlation across Complexes', fontsize=13, fontweight='bold')
    ax1.grid(axis='y', linestyle='--', alpha=0.5)
    
    # 为柱子加上具体数值
    for bar in bars1:
        yval = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, yval + 0.002, f'{yval:.4f}', ha='center', va='bottom', fontsize=8)

    # --- 图 2: 总面积相对误差 (证明总量计算精度) ---
    bars2 = ax2.bar(x, errors, color='#dc2626', alpha=0.8, width=0.6)
    ax2.set_ylim(0, max(errors) * 1.25)
    ax2.set_xticks(x)
    ax2.set_xticklabels(complex_names, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Relative Error (%)', fontsize=11, fontweight='bold')
    ax2.set_title('Total SASA Relative Error across Complexes', fontsize=13, fontweight='bold')
    ax2.grid(axis='y', linestyle='--', alpha=0.5)
    
    # 为柱子加上具体数值
    for bar in bars2:
        yval = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, yval + (max(errors)*0.05), f'{yval:.2f}%', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path)
    print(f"\n[+] 批量对比可视化图表已保存至: {output_path}")

def main():
    # 参数设置
    pdb_dir = PROJECT_ROOT / "data" / "raw" / "pdb_complexes"
    dot_file = EXAMPLE_DATA_DIR / "Dot.txt"
    plot_output = EXAMPLE_DATA_DIR / "batch_sasa_validation.png"
    solvent_radius = 1.4
    max_test_files = 10  # 选取前10个复合物进行测试（如果全跑可能需要几分钟）

    print("=" * 60)
    print("🚀 启动自动化批量验证管线 (Custom vs FreeSASA)")
    print("=" * 60)

    pdb_files = list(pdb_dir.glob("*.pdb"))[:max_test_files]
    dots = load_dots(dot_file)

    if not pdb_files:
        print(f"❌ 未在 {pdb_dir} 找到 PDB 文件！")
        return

    # 数据收集列表
    complex_names = []
    my_totals = []
    fs_totals = []
    errors = []
    correlations = []

    for i, pdb_file in enumerate(pdb_files, 1):
        complex_name = pdb_file.stem
        print(f"[{i}/{len(pdb_files)}] 正在评测: {complex_name} ...", end="", flush=True)
        
        start_time = time.time()
        
        # 自研算法计算
        atoms = parse_pdb(pdb_file)
        my_total = calculate_sasa(atoms, dots, solvent_radius)
        my_res_dict = aggregate_residue_sasa(atoms)
        
        # FreeSASA计算
        fs_total, fs_res_dict = get_freesasa_results(pdb_file)
        
        # 对齐残基求相关性
        my_aligned = []
        fs_aligned = []
        for key, my_area in my_res_dict.items():
            chain_id, res_id = key[0], key[1]
            if (chain_id, res_id) in fs_res_dict:
                my_aligned.append(my_area)
                fs_aligned.append(fs_res_dict[(chain_id, res_id)])
                
        # 计算指标
        rel_error = abs(my_total - fs_total) / fs_total * 100
        corr, _ = pearsonr(my_aligned, fs_aligned)
        
        # 存入列表
        complex_names.append(complex_name)
        my_totals.append(my_total)
        fs_totals.append(fs_total)
        errors.append(rel_error)
        correlations.append(corr)
        
        elapsed = time.time() - start_time
        print(f" 完成! (耗时 {elapsed:.2f}s, Error: {rel_error:.2f}%, Corr: {corr:.4f})")

    # --- 输出能在报告里直接用的 Markdown 表格 ---
    print("\n\n" + "=" * 60)
    print("📊 可直接复制进报告的 Markdown 验证表格：")
    print("=" * 60)
    print("| PDB ID | 自研总面积 (Å²) | FreeSASA (Å²) | 相对误差 | 残基级 Pearson R |")
    print("|--------|-----------------|---------------|----------|------------------|")
    for i in range(len(complex_names)):
        print(f"| {complex_names[i]:<6} | {my_totals[i]:<15.2f} | {fs_totals[i]:<13.2f} | {errors[i]:>6.2f}% | {correlations[i]:>16.4f} |")
    
    # 绘制结果统计图
    plot_batch_results(complex_names, correlations, errors, plot_output)

if __name__ == "__main__":
    main()