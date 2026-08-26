#!/usr/bin/env python3
"""
从 AChE-肽复合物 PDB (如 alllhrc_complex.pdb) 提取纯 AChE 单体 (仅链 A),
去除肽链 (链 B)、结晶水与离子, 生成独立的单体对照组结构 input/ache.pdb。

v2.7.3 修复: 旧版纯 Python 后备逻辑 `chain_id == "A" or resseq <= 530` 会把
链 B 中编号 <= 530 的残基 (肽 1-7 号) 也保留下来, 导致生成的 ache.pdb 有
两条链 (不是单体对照)。现在的规则:
  - 有链 ID 的 PDB: 只保留链 A (肽在链 B, 编号无关);
  - 完全没有链 ID 的 PDB: 视为单链, 保留全部蛋白原子;
  - 输出后强制校验必须是单链, 否则报错返回失败;
  - 输出带 TER / END, 链间不会粘连。

用法:
    python3 extract_ache_monomer.py [输入复合物PDB, 默认 ../input/alllhrc_complex.pdb] [输出单体PDB, 默认 ../input/ache.pdb]
"""
import sys
from pathlib import Path

WATER_ION_NAMES = {
    "SOL", "HOH", "WAT", "H2O", "NA", "CL", "NA+", "CL-", "K", "K+",
    "MG", "MG2", "ZN", "CA", "CA2", "CL+", "NA1", "CL1",
}


def _is_atom_line(line: str) -> bool:
    return line.startswith(("ATOM  ", "HETATM"))


def _resname(line: str) -> str:
    return line[17:20].strip().upper()


def _chain_id(line: str) -> str:
    return line[21:22].strip()


def _count_chains_and_residues(path: Path):
    """返回 (链ID集合, 残基编号集合), 只统计蛋白/配体, 排除水与离子。"""
    chains = set()
    resids = set()
    if not path.exists():
        return chains, resids
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not _is_atom_line(line):
                continue
            if _resname(line) in WATER_ION_NAMES:
                continue
            chains.add(_chain_id(line))
            resids.add(line[22:26].strip())
    return chains, resids


def extract_ache_pdb(src_pdb: Path, dst_pdb: Path) -> bool:
    if not src_pdb.exists():
        print(f"!!! 错误: 找不到参考复合物 PDB: {src_pdb}", file=sys.stderr)
        return False

    print(f">> 正在从复合物 {src_pdb} 提取纯 AChE 单体 (仅链 A) ...")

    # ---------- 路径 1: MDAnalysis ----------
    try:
        import MDAnalysis as mda
        u = mda.Universe(str(src_pdb))
        ache = u.select_atoms("chainID A")
        if len(ache) == 0:
            ache = u.select_atoms("segid A")
        if len(ache) == 0:
            # 完全没有链 ID: 视为单链, 取全部蛋白 (输出校验会兜底)
            ache = u.select_atoms("protein")
        if len(ache) > 0:
            dst_pdb.parent.mkdir(parents=True, exist_ok=True)
            ache.write(str(dst_pdb))
            chains, resids = _count_chains_and_residues(dst_pdb)
            if len(chains) == 1:
                print(f">> [MDAnalysis 提取成功] 已提取 AChE 单体 ({len(ache)} 原子, "
                      f"{len(resids)} 残基, 链 {sorted(chains)}) -> {dst_pdb}")
                return True
            print(f"!! [MDAnalysis 结果校验失败] 输出仍有 {len(chains)} 条链, 改用纯 Python 提取 ...")
    except Exception as e:
        print(f">> [提示] MDAnalysis 提取未完成 ({e})，切换为标准 Python 纯解析 PDB 语法 ...")

    # ---------- 路径 2: 标准纯 Python 逐行解析 (免额外依赖) ----------
    with open(src_pdb, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    has_chain_ids = any(_chain_id(line) for line in lines if _is_atom_line(line))

    def keep(line: str) -> bool:
        if _resname(line) in WATER_ION_NAMES:
            return False
        if not has_chain_ids:
            return True                    # 单链 (无链 ID): 全部保留
        return _chain_id(line) == "A"      # 有链 ID: 只保留链 A (肽在链 B)

    out_lines = []
    atom_count = 0
    for line in lines:
        if line.startswith("CRYST1"):
            out_lines.append(line)
        elif _is_atom_line(line) and keep(line):
            out_lines.append(line)
            atom_count += 1

    if atom_count == 0:
        print("!!! 错误: 未能在复合物文件中成功选择到有效的 AChE 链 A 原子！", file=sys.stderr)
        return False

    # 保证 TER / END 存在 (防止链粘连)
    while out_lines and out_lines[-1].strip() in ("TER", "END"):
        out_lines.pop()
    out_lines.append("TER\n")
    out_lines.append("END\n")

    dst_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(dst_pdb, "w", encoding="utf-8") as fo:
        fo.writelines(out_lines)

    # ---------- 输出校验: 必须是单链 ----------
    chains, resids = _count_chains_and_residues(dst_pdb)
    if len(chains) != 1:
        print(f"!!! 错误: 提取结果仍有 {len(chains)} 条链 ({sorted(chains)}) —— 不是单体对照！"
              f"请检查 {src_pdb} 的链 ID 是否正确 (A=AChE, B=肽)。", file=sys.stderr)
        return False

    print(f">> [纯 Python 提取成功] 已提取 AChE 单体 ({atom_count} 原子, "
          f"{len(resids)} 残基, 链 {sorted(chains)}) -> {dst_pdb}")
    return True


if __name__ == "__main__":
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../input/alllhrc_complex.pdb")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("../input/ache.pdb")
    success = extract_ache_pdb(src, dst)
    sys.exit(0 if success else 1)
