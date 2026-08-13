#!/usr/bin/env python3
"""
Backward-compatible wrapper around compute_peptide_ss.parse_ss_dat.

用法:
    python3 dssp_bins.py ss_pep.dat output.dat [window_ns]
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_peptide_ss import (  # noqa: E402
    auto_window_ns,
    parse_ss_dat,
    seqs_to_fracs,
    write_outputs,
)


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("ss_pep.dat")
    window_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    times, seqs = parse_ss_dat(src, nres=7)
    if not seqs:
        print(">> 提示: 二级结构序列数据为空，跳过统计。")
        return
    if window_ns <= 0:
        window_ns = auto_window_ns(float(times[0]), float(times[-1]))
    fracs = seqs_to_fracs(seqs)
    write_outputs(src.parent if src.parent.as_posix() != "" else Path("."),
                  times, fracs, seqs, window_ns,
                  source=f"dssp_bins.py:{src.name}")


if __name__ == "__main__":
    main()
