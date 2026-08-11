#!/usr/bin/env bash
# ============================================================
# 图4: 肽(Aβ) 二级结构倾向随时间的演变
# 论文 3.2 节 (helix, beta-turn, bend 倾向, 每 50 ns 统计)
# gmx do_dssp 基于 DSSP 判定二级结构, -ssdump 输出逐帧 SS 字符串
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

for ana_dir in "/mnt/f/anaconda" "/mnt/f/Anaconda" "/f/anaconda" "/f/Anaconda" "/mnt/c/anaconda3" "/c/anaconda3" "$HOME/anaconda3" "$HOME/Anaconda3"; do
    if [ -d "$ana_dir" ]; then
        export PATH="$ana_dir:$ana_dir/Scripts:$ana_dir/Library/bin:$PATH"
        break
    fi
done

PY_FOUND=""
if [ -n "${PYTHON:-}" ]; then
    PY_FOUND="${PYTHON}"
else
    for cand in \
        "/mnt/f/anaconda/python.exe" \
        "/mnt/f/Anaconda/python.exe" \
        "/mnt/c/anaconda3/python.exe" \
        "/f/anaconda/python.exe" \
        "/f/Anaconda/python.exe" \
        "/c/anaconda3/python.exe" \
        "python.exe" \
        "python" \
        "python3"; do
        if command -v "$cand" >/dev/null 2>&1 || [ -x "$cand" ]; then
            if "$cand" -c "import numpy" >/dev/null 2>&1; then
                PY_FOUND="$cand"
                break
            fi
        fi
    done
fi
PY="${PY_FOUND:-python}"

if ! grep -q "\[ *Peptide *\]" index.ndx 2>/dev/null; then
    echo ">> [单体蛋白模式] 未检测到 Peptide 组，跳过小肽二级结构计算。"
    exit 0
fi

# 计算肽的二级结构随时间的分布
if gmx dssp -h 2>&1 | grep -q -- "-sel "; then
    # GROMACS >= 2023 新版内置 dssp (支持 -sel Peptide -o .dat -num .xvg)
    gmx dssp -s md.tpr -f md.xtc -n index.ndx -sel "Peptide" -o ss_pep.dat -num ss_pep_num.xvg
    SS_FILE="ss_pep.dat"
else
    # GROMACS <= 2022 传统版 (do_dssp / dssp -sc / -ssdump)
    if gmx help dssp >/dev/null 2>&1; then
        DSSP_CMD="dssp"
    else
        DSSP_CMD="do_dssp"
    fi
    if gmx help ${DSSP_CMD} 2>&1 | grep -q -- "-ssdump"; then
        SC_FLAG="-ssdump"
    else
        SC_FLAG="-sc"
    fi
    gmx ${DSSP_CMD} -s md.tpr -f md.xtc -n index.ndx \
                ${SC_FLAG} ss_pep.sc -o ss_pep.xpm << EOF
Peptide
EOF
    SS_FILE="ss_pep.sc"
fi

# 将二级结构输出分箱统计每 50 ns 的 helix / turn / bend 倾向 (测试模式下自适应窗口)
TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    WIN_NS=0.05
else
    WIN_NS=50
fi

if command -v wslpath >/dev/null 2>&1 && [[ "${PY}" == *".exe"* || "${PY}" == *"/mnt/"* || "${PY}" == *":"* ]]; then
    SCRIPT_PATH=$(wslpath -w "$(dirname "$0")/dssp_bins.py" 2>/dev/null || echo "$(dirname "$0")/dssp_bins.py")
else
    SCRIPT_PATH="$(dirname "$0")/dssp_bins.py"
fi
"${PY}" "${SCRIPT_PATH}" "${SS_FILE}" ss_pep_bins.dat ${WIN_NS}

echo "二级结构输出: ${SS_FILE}, ss_pep_bins.dat"
echo "ss_pep_bins.dat 列: 时间窗口(ns)  helix倾向  turn倾向  bend倾向"
