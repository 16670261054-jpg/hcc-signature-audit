#!/usr/bin/env bash
# M9 正式运行链（预注册 N20）：等官方等价性验证结束 → 只有 WT_identical 且 DR_identical 都为 TRUE 才开跑 2124 个基因
# 同一时间只开一个 R：等待期间不启动任何 R
ROOT="$(cd "${HCC_AUDIT_ROOT:-.}" && pwd)"; export HCC_AUDIT_ROOT="$ROOT"
export HCC_AUDIT_DATA="${HCC_AUDIT_DATA:-$ROOT/data/raw}"
L="$ROOT/analysis"
O="$ROOT/N11_out/hiend_m9_hvg"
CHK="$O/equivalence_check_ALAD.tsv"
while ! grep -q "^exit " "$L/hiend_m9_knk_official_check.log" 2>/dev/null; do sleep 300; done
if ! grep -q "^exit 0" "$L/hiend_m9_knk_official_check.log"; then echo "official check failed -> formal run NOT started"; exit 2; fi
[ -f "$CHK" ] || { echo "no equivalence table -> formal run NOT started"; exit 3; }
cat "$CHK"
WI=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR==2{print $h["WT_identical"]}' "$CHK")
DI=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR==2{print $h["DR_identical"]}' "$CHK")
echo "WT_identical=$WI DR_identical=$DI"
if [ "$WI" != "TRUE" ] || [ "$DI" != "TRUE" ]; then echo "NOT equivalent -> formal run NOT started (per N20)"; exit 4; fi
cd /tmp && KNK_CORES=4 KNK_KODIR=ko_formal KNK_KOFILE="$O/ko_list_formal.tsv" \
  "${RSCRIPT:-Rscript}" "$L/hiend_m9_knk_hvg.R" > "$L/hiend_m9_knk_formal.log" 2>&1
S=$?; echo "exit $S" >> "$L/hiend_m9_knk_formal.log"; exit $S
