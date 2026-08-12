#!/usr/bin/env python3
"""
将复合物 PDB (如 alllhrc_complex.pdb) 自动拆分为单独的受体蛋白 (pro.pdb) 和多肽配体 (tai.pdb)，
严格配合用户的“分体生拓扑 + gro合并 + topol.top追加”的 12 步动力学模拟标准工作流。

用法:
    python3 split_complex_pdb.py <输入复合物PDB> <输出蛋白PDB: pro.pdb> <输出配体PDB: tai.pdb>
"""
import sys
from pathlib import Path

def split_complex(src_pdb: Path, pro_pdb: Path, tai_pdb: Path) -> bool:
    if not src_pdb.exists():
        print(f"!!! 错误: 找不到被拆分的复合物文件: {src_pdb}", file=sys.stderr)
        return False

    print(f">> 正在将复合物 {src_pdb} 自动拆解为单独的受体蛋白 pro.pdb 和肽段 tai.pdb ...")

    pro_lines = []
    tai_lines = []
    pro_atoms = 0
    tai_atoms = 0

    with open(src_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("CRYST1"):
                pro_lines.append(line)
                tai_lines.append(line)
            elif line.startswith(("ATOM  ", "HETATM")):
                resname = line[17:20].strip()
                chain_id = line[21:22]
                try:
                    resseq = int(line[22:26].strip())
                except ValueError:
                    resseq = 9999
                # 排除水与离子
                if resname in ("SOL", "HOH", "WAT", "NA", "CL", "Na+", "Cl-"):
                    continue
                # 分配至蛋白 (1-530) 或 多肽 (531-537 / 链B)
                if chain_id == "A" or resseq <= 530:
                    pro_lines.append(line)
                    pro_atoms += 1
                elif chain_id == "B" or resseq > 530:
                    tai_lines.append(line)
                    tai_atoms += 1
            elif line.startswith("END"):
                pro_lines.append(line)
                tai_lines.append(line)

    pro_pdb.parent.mkdir(parents=True, exist_ok=True)
    tai_pdb.parent.mkdir(parents=True, exist_ok=True)

    with open(pro_pdb, "w", encoding="utf-8") as fo:
        fo.writelines(pro_lines)
        if not pro_lines[-1].startswith("END"):
            fo.write("END\n")

    with open(tai_pdb, "w", encoding="utf-8") as fo:
        fo.writelines(tai_lines)
        if not tai_lines[-1].startswith("END"):
            fo.write("END\n")

    print(f">> [拆分完成] 受体蛋白 ({pro_atoms} 原子) -> {pro_pdb}")
    print(f">> [拆分完成] 多肽配体 ({tai_atoms} 原子) -> {tai_pdb}")
    return pro_atoms > 0 and tai_atoms > 0


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../input/alllhrc_complex.pdb")
    pro = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("pro.pdb")
    tai = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("tai.pdb")
    success = split_complex(src, pro, tai)
    sys.exit(0 if success else 1)
