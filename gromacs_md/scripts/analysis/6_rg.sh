#!/usr/bin/env bash
# ============================================================
# Radius of gyration (complex / AChE). Inserted before H-bond analysis.
#
# Robustness (v2.7.3):
#   - group fallbacks Protein -> System -> Backbone -> AChE
#   - output is only accepted when it really contains >=2 numeric data rows
#     (GROMACS can leave an empty/header-only xvg behind on failure)
#   - stale empty files are removed before each attempt so they can never
#     "poison" later runs that only check file existence
# ============================================================
set -u

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
elif command -v gmx >/dev/null 2>&1; then
    gmx() { gmx "$@"; }
else
    echo "!! 找不到 gmx / gmx.exe (PATH 里没有 GROMACS)"
    exit 1
fi

TPR_FILE="md.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="md_0_1.tpr"
fi
TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md_0_1.xtc"
    fi
fi

if [ ! -f "${TPR_FILE}" ] || [ ! -f "${TRAJ_FILE}" ]; then
    echo "!! 缺少 ${TPR_FILE} 或 ${TRAJ_FILE}，无法计算 Rg"
    exit 1
fi

xvg_has_data() {
    [ -f "$1" ] || return 1
    local n
    n=$(grep -vE '^[#@;]' "$1" | awk 'NF>=2 && $1 ~ /^[+-]?[0-9.]+$/ && $2 ~ /^[+-]?[0-9.]+$/ {c++} END{print c+0}')
    [ "${n:-0}" -ge 2 ]
}

try_gyrate() {
    local group="$1" out="$2"
    rm -f "$out"
    echo ">> [Rg] gmx gyrate (group ${group}) -> ${out} ..."
    if [ -f index.ndx ]; then
        gmx gyrate -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -o "$out" << EOF || return 1
${group}
EOF
    else
        gmx gyrate -s "${TPR_FILE}" -f "${TRAJ_FILE}" -o "$out" << EOF || return 1
${group}
EOF
    fi
    xvg_has_data "$out"
}

# ---- complex Rg (whole protein) ----
OK=0
for g in Protein System Backbone AChE; do
    if grep -q "\[ *${g} *\]" index.ndx 2>/dev/null || [ ! -f index.ndx ]; then
        if try_gyrate "$g" gyrate_complex.xvg; then
            OK=1
            break
        fi
    fi
done
if [ "${OK}" = "1" ]; then
    echo ">> [Rg] gyrate_complex.xvg OK"
else
    echo "!! [Rg] gyrate_complex.xvg 未能生成 (gmx gyrate 全部分组失败)"
fi

# ---- AChE-only Rg (if an AChE group exists) ----
if grep -q "\[ *AChE *\]" index.ndx 2>/dev/null; then
    OK=0
    for g in AChE Protein; do
        if try_gyrate "$g" gyrate_ache.xvg; then
            OK=1
            break
        fi
    done
    [ "${OK}" = "1" ] && echo ">> [Rg] gyrate_ache.xvg OK"
fi

echo "Rg output: gyrate_complex.xvg (and gyrate_ache.xvg if AChE group exists)"
