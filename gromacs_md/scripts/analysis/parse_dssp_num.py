#!/usr/bin/env python3
"""
Backward-compatible wrapper. Real parsing lives in compute_peptide_ss.py.

用法:
    python3 parse_dssp_num.py [ss_pep_num.xvg] [ss_pep_bins.dat] [窗口ns] [残基数]
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_peptide_ss import (  # noqa: E402
    auto_window_ns,
    parse_ss_num_xvg,
    write_outputs,
)


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("ss_pep_num.xvg")
    # dst is ignored: write_outputs always uses ss_pep_bins.dat next to src
    window_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    nres = int(float(sys.argv[4])) if len(sys.argv) > 4 else 7

    if not src.exists():
        print(f"!!! 提示: 未找到 {src}，无法执行分箱统计。")
        return

    times, fracs = parse_ss_num_xvg(src, nres=nres)
    if times is None or fracs is None:
        print(f"!!! 提示: {src} 中没有有效数据点。")
        return
    if window_ns <= 0:
        window_ns = auto_window_ns(float(times[0]), float(times[-1]))
    write_outputs(src.parent if src.parent.as_posix() != "" else Path("."),
                  times, fracs, seqs=[], window_ns=window_ns,
                  source=f"parse_dssp_num.py:{src.name}")


if __name__ == "__main__":
    main()
