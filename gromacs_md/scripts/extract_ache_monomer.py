#!/usr/bin/env python3
"""
自动从 AChE-Aβ 复合物 PDB (如 alllhrc_complex.pdb) 中提取纯 AChE 单体蛋白结构，
去除结合肽 (如 Aβ 7肽) 及结晶水，生成独立的单体对照组结构 input/ache.pdb。

用法:
    python3 extract_ache_monomer.py [输入复合物PDB, 默认 ../input/alllhrc_complex.pdb] [输出单体PDB, 默认 ../input/ache.pdb]
"""
import sys
from pathlib import Path

def extract_ache_pdb(src_pdb: Path, dst_pdb: Path) -> bool:
    if not src_pdb.exists():
        print(f"!!! 错误: 找不到参考复合物 PDB: {src_pdb}", file=sys.stderr)
        return False

    print(f">> 正在从复合物 {src_pdb} 提取纯 AChE 单体蛋白结构 ...")

    # 优先尝试使用 MDAnalysis
    try:
        import MDAnalysis as mda
        u = mda.Universe(str(src_pdb))
        ache = u.select_atoms("chainID A")
        if len(ache) == 0:
            ache = u.select_atoms("segid A")
        if len(ache) == 0:
            ache = u.select_atoms("protein and resid 1-530")
        if len(ache) > 0:
            dst_pdb.parent.mkdir(parents=True, exist_ok=True)
            ache.write(str(dst_pdb))
            print(f">> [MDAnalysis 提取成功] 已提取 AChE 单体 ({len(ache)} 个原子) 保存至: {dst_pdb}")
            return True
    except Exception as e:
        print(f">> [提示] MDAnalysis 提取未完成 ({e})，切换为标准 Python 纯解析 PDB 语法 ...")

    # 标准纯 Python 逐行解析 PDB 语法 (免额外依赖)
    out_lines = []
    atom_count = 0
    with open(src_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("CRYST1"):
                out_lines.append(line)
            elif line.startswith(("ATOM  ", "HETATM")):
                resname = line[17:20].strip()
                chain_id = line[21:22]
                try:
                    resseq = int(line[22:26].strip())
                except ValueError:
                    resseq = 9999
                # 过滤常见水分子与离子
                if resname in ("SOL", "HOH", "WAT", "NA", "CL", "Na+", "Cl-"):
                    continue
                # 选择链 A 或前 530 个蛋白质残基
                if chain_id == "A" or resseq <= 530:
                    out_lines.append(line)
                    atom_count += 1
            elif line.startswith("END"):
                out_lines.append(line)

    if atom_count > 0:
        dst_pdb.parent.mkdir(parents=True, exist_ok=True)
        with open(dst_pdb, "w", encoding="utf-8") as fo:
            fo.writelines(out_lines)
            if not out_lines[-1].startswith("END"):
                fo.write("END\n")
        print(f">> [标准 Python 提取成功] 已提取 AChE 单体 ({atom_count} 个原子) 保存至: {dst_pdb}")
        return True

    print("!!! 错误: 未能在复合物文件中成功选择到有效 AChE 蛋白质原子！", file=sys.stderr)
    return False


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../input/alllhrc_complex.pdb")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("../input/ache.pdb")
    success = extract_ache_pdb(src, dst)
    sys.exit(0 if success else 1)
