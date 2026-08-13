#!/usr/bin/env bash
# ============================================================
# 图4: 肽二级结构倾向随时间的演变
# v2.4: GROMACS 2025 正确选区 + 逐帧 DSSP + 1 ns 分箱
#       gmx dssp 失败时自动回退到 Kabsch-Sander / Ramachandran
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

TPR_FILE="md_0_1.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="md.tpr"
fi

TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md_0_1.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md.xtc"
    fi
fi

run_gmx_dssp() {
    local sel="$1"
    echo ">> [gmx dssp] -sel ${sel}  -hmode dssp -clear  (只分析 Peptide, 避开 4ey6 三段受体)"
    gmx dssp -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx \
        -sel "${sel}" -o ss_pep.dat -num ss_pep_num.xvg \
        -hmode dssp -clear
}

# gmx dssp 失败不得中断整个分析: 后面的 compute_peptide_ss.py 会回退
set +e
DSSP_OK=0
if gmx dssp -h 2>&1 | grep -q -- "-sel "; then
    for SEL in 'group "Peptide"' 'group Peptide' 'Peptide' 'resid 531 to 537'; do
        if run_gmx_dssp "${SEL}"; then
            NROW=$(grep -v '^[#@;]' ss_pep_num.xvg 2>/dev/null | grep -c '[0-9]' || true)
            if [ "${NROW}" -ge 10 ]; then
                echo ">> [gmx dssp] 成功  sel=${SEL}  帧数=${NROW}"
                DSSP_OK=1
                break
            fi
            echo ">> [gmx dssp] 输出帧数过少 (${NROW})，换一种 -sel 再试 ..."
        else
            echo ">> [gmx dssp] 命令失败，换一种 -sel 再试 ..."
        fi
    done
else
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
    gmx ${DSSP_CMD} -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx \
                ${SC_FLAG} ss_pep.sc -o ss_pep.xpm << EOF
Peptide
EOF
    if [ -f ss_pep.sc ]; then
        DSSP_OK=1
    fi
fi
set -e

TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    WIN_NS=0.05
else
    # 100 ns 轨迹用 1 ns 窗 (~100 点); 绝不能再用 50 ns (只会得到 2 个点)
    WIN_NS=1
fi

to_py_path() {
    local p="$1"
    if command -v wslpath >/dev/null 2>&1 && [[ "${PY}" == *".exe"* || "${PY}" == *"/mnt/"* || "${PY}" == *":"* ]]; then
        wslpath -w "$p" 2>/dev/null || echo "$p"
    else
        echo "$p"
    fi
}

SS_SCRIPT="$(dirname "$0")/compute_peptide_ss.py"
echo ">> [DSSP rebuild] 逐帧解析 / 必要时 Kabsch-Sander 回退 ..."
"${PY}" "$(to_py_path "${SS_SCRIPT}")" --dir . --nres 7 --window-ns "${WIN_NS}" || {
    echo "!!! compute_peptide_ss.py 返回非零。若已有 ss_pep_frac.xvg 仍可继续绘图。"
}

echo "二级结构输出: ss_pep_frac.xvg (逐帧), ss_pep_bins.dat (${WIN_NS} ns 窗), ss_pep_perres.dat, ss_pep_summary.txt"
