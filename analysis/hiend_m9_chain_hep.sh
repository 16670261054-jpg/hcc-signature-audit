#!/usr/bin/env bash
# M9 方案 A 运行链（2026-10-17，用户选 A）：同一时间只开一个 R，逐步串行
#  0 等官方等价性验证（全组织网络，ALAD）结束
#  1 全组织测速两个 KO 的判别性检查（留档）
#  2 肿瘤肝细胞单独建网（BUILD_ONLY）
#  3 试点 20 个基因 KO → 判别性检查；不 PASS 就停（用户要求），不跑正式
#  4 PASS 才跑正式（网络内全部签名基因）→ 正式结果再做一次判别性检查
ROOT="$(cd "${HCC_AUDIT_ROOT:-.}" && pwd)"; export HCC_AUDIT_ROOT="$ROOT"
export HCC_AUDIT_DATA="${HCC_AUDIT_DATA:-$ROOT/data/raw}"
L="$ROOT/analysis"
N="$ROOT/N11_out"
RS="${RSCRIPT:-Rscript}"
PY="${PYTHON:-python}"
run() { local log=$1; shift; "$@" > "$L/$log" 2>&1; local s=$?; echo "exit $s" >> "$L/$log"; echo "[$(date '+%F %T')] $log exit $s"; return $s; }
cd /tmp
while ! grep -q "^exit " "$L/hiend_m9_knk_official_check.log" 2>/dev/null; do sleep 300; done
echo "[$(date '+%F %T')] official check finished: $(grep '^exit ' $L/hiend_m9_knk_official_check.log)"
KNK_OUT=hiend_m9_hvg KNK_KODIR=ko_speedtest run hiend_m9_discrim_mixed.log "$RS" "$L/hiend_m9_discrim_check.R"
KNK_OUT=hiend_m9_hep KNK_CELLTYPE=Hepatocyte KNK_BUILD_ONLY=1 KNK_CORES=8 run hiend_m9_hep_build.log "$RS" "$L/hiend_m9_knk_hvg.R" || exit 11
run hiend_m9_hep_kolist_pilot.log "$PY" -I "$L/hiend_m9_make_kolist.py" hiend_m9_hep pilot || exit 12
KNK_OUT=hiend_m9_hep KNK_CELLTYPE=Hepatocyte KNK_KODIR=ko_pilot KNK_KOFILE="$N/hiend_m9_hep/ko_list_pilot.tsv" KNK_CORES=4 \
  run hiend_m9_hep_pilot.log "$RS" "$L/hiend_m9_knk_hvg.R" || exit 13
KNK_OUT=hiend_m9_hep KNK_KODIR=ko_pilot run hiend_m9_discrim_hep_pilot.log "$RS" "$L/hiend_m9_discrim_check.R" || exit 14
P=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR==2{print $h["PASS"]}' "$N/hiend_m9_hep/discrim_check_ko_pilot.tsv")
echo "[$(date '+%F %T')] pilot PASS=$P"
if [ "$P" != "TRUE" ]; then echo "pilot NOT discriminative -> STOP (formal not started)"; exit 20; fi
run hiend_m9_hep_kolist_formal.log "$PY" -I "$L/hiend_m9_make_kolist.py" hiend_m9_hep formal || exit 21
KNK_OUT=hiend_m9_hep KNK_CELLTYPE=Hepatocyte KNK_KODIR=ko_formal KNK_KOFILE="$N/hiend_m9_hep/ko_list_formal.tsv" KNK_CORES=4 \
  run hiend_m9_hep_formal.log "$RS" "$L/hiend_m9_knk_hvg.R" || exit 22
KNK_OUT=hiend_m9_hep KNK_KODIR=ko_formal run hiend_m9_discrim_hep_formal.log "$RS" "$L/hiend_m9_discrim_check.R"
echo "[$(date '+%F %T')] ALL DONE"
