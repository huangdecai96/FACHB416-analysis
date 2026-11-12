from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(r"/PROJECT_ROOT/2026-07-03\new-chat-2")
DTU_DIR = ROOT / "work" / "round15_dtu_results_local"
OUT_DIR = ROOT / "outputs"

BASE_TABLE = OUT_DIR / "FACHB416_paper1_round11_s0_dbcan_interpro_motif_integrated_zh_20260705.tsv"

OUT_DETAIL = OUT_DIR / "FACHB416_paper1_round16_S0_DTU_localization_integrated_zh_20260705.tsv"
OUT_TOOL_SUMMARY = OUT_DIR / "FACHB416_paper1_round16_S0_DTU_localization_tool_counts_zh_20260705.tsv"
OUT_FIGURE_SOURCE = OUT_DIR / "FACHB416_paper1_round16_S0_DTU_localization_figure_source_data_20260705.tsv"
OUT_REPORT = OUT_DIR / "FACHB416_paper1_round16_S0_DTU_localization_summary_zh_20260705.md"
OUT_LOG = OUT_DIR / "FACHB416_paper1_round16_S0_DTU_localization_run_log_20260705.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_tsv_utf8(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def norm_gene(protein_id: str) -> str:
    protein_id = (protein_id or "").strip()
    if not protein_id:
        return ""
    return re.sub(r"\.t\d+$", "", protein_id).rstrip(".")


def fmt_prob(value: str | float | None) -> str:
    if value in (None, ""):
        return ""
    try:
        return f"{float(value):.3f}"
    except Exception:
        return str(value)


def parse_signalp(path: Path) -> dict[str, dict[str, str]]:
    data: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            pid = parts[0]
            data[pid] = {
                "SignalP6预测": parts[1],
                "SignalP6_OTHER概率": fmt_prob(parts[2]),
                "SignalP6_SP概率": fmt_prob(parts[3]),
                "SignalP6切割位点": parts[4] if len(parts) > 4 else "",
            }
    return data


def parse_targetp(path: Path) -> dict[str, dict[str, str]]:
    data: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            pid = parts[0]
            data[pid] = {
                "TargetP2预测": parts[1],
                "TargetP2_OTHER概率": fmt_prob(parts[2]),
                "TargetP2_SP概率": fmt_prob(parts[3]),
                "TargetP2_mTP概率": fmt_prob(parts[4]),
                "TargetP2_cTP概率": fmt_prob(parts[5]),
                "TargetP2_lTP概率": fmt_prob(parts[6]),
                "TargetP2切割位点": parts[7] if len(parts) > 7 else "",
            }
    return data


def parse_deeploc(path: Path) -> dict[str, dict[str, str]]:
    numeric_cols = [
        "Cytoplasm",
        "Nucleus",
        "Extracellular",
        "Cell membrane",
        "Mitochondrion",
        "Plastid",
        "Endoplasmic reticulum",
        "Lysosome/Vacuole",
        "Golgi apparatus",
        "Peroxisome",
        "Peripheral",
        "Transmembrane",
        "Lipid anchor",
        "Soluble",
    ]
    data: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            pid = row["Protein_ID"]
            entry = {
                "DeepLoc2定位": row.get("Localizations", ""),
                "DeepLoc2信号": row.get("Signals", ""),
                "DeepLoc2膜类型": row.get("Membrane types", ""),
            }
            for col in numeric_cols:
                entry[f"DeepLoc2_{col}概率"] = fmt_prob(row.get(col, ""))
            data[pid] = entry
    return data


def parse_deeptmhmm_gff(path: Path) -> dict[str, dict[str, object]]:
    data: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "DeepTMHMM_TMR数量": "",
            "DeepTMHMM信号肽区间": "",
            "DeepTMHMM_TMR区间": [],
        }
    )
    current_pid = ""
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("##") or line == "//":
                continue
            if line.startswith("#"):
                m_len = re.match(r"#\s+(\S+)\s+Length:\s+(\d+)", line)
                if m_len:
                    current_pid = m_len.group(1)
                    data[current_pid]["蛋白长度"] = m_len.group(2)
                m_tmr = re.match(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)", line)
                if m_tmr:
                    current_pid = m_tmr.group(1)
                    data[current_pid]["DeepTMHMM_TMR数量"] = m_tmr.group(2)
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            pid, feature, start, end = parts[:4]
            if feature == "signal":
                data[pid]["DeepTMHMM信号肽区间"] = f"{start}-{end}"
            elif feature == "TMhelix":
                data[pid]["DeepTMHMM_TMR区间"].append(f"{start}-{end}")
            current_pid = pid
    for pid, entry in data.items():
        entry["DeepTMHMM_TMR区间"] = ";".join(entry["DeepTMHMM_TMR区间"])
    return dict(data)


def parse_deeptmhmm_header_class(path: Path) -> dict[str, str]:
    classes: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith(">"):
                continue
            m = re.match(r">(\S+)\s+\|\s+(.+)$", line)
            if m:
                classes[m.group(1)] = m.group(2).strip()
    return classes


def choose_protein_id(gene: str, base_pid: str, tool_ids_by_gene: dict[str, list[str]]) -> str:
    clean_base = (base_pid or "").strip()
    if re.search(r"\.t\d+$", clean_base):
        return clean_base
    ids = tool_ids_by_gene.get(gene, [])
    if ids:
        return sorted(ids)[0]
    return clean_base


def contains(text: str, token: str) -> bool:
    return token.lower() in (text or "").lower()


def evidence_call(row: dict[str, object]) -> tuple[str, str, str]:
    dbcan = str(row.get("dbCAN证据等级", ""))
    s0_class = str(row.get("S0类别", ""))
    cazy = str(row.get("CAZyme家族", ""))
    signalp_sp = row.get("SignalP6预测") == "SP"
    targetp_pred = str(row.get("TargetP2预测", ""))
    targetp_sp = targetp_pred == "SP"
    targetp_organelle = targetp_pred in {"mTP", "cTP", "lTP"}
    tm_signal = bool(row.get("DeepTMHMM信号肽区间"))
    try:
        tm_count = int(row.get("DeepTMHMM_TMR数量") or 0)
    except ValueError:
        tm_count = 0
    dl_loc = str(row.get("DeepLoc2定位", ""))
    dl_sig = str(row.get("DeepLoc2信号", ""))
    dl_mem = str(row.get("DeepLoc2膜类型", ""))
    extracellular = contains(dl_loc, "Extracellular")
    golgi = contains(dl_loc, "Golgi")
    cell_membrane = contains(dl_loc, "Cell membrane")
    plastid = contains(dl_loc, "Plastid")
    mitochondrion = contains(dl_loc, "Mitochondrion")
    dl_signal = contains(dl_sig, "Signal peptide")
    dl_tm = contains(dl_sig, "Transmembrane") or contains(dl_mem, "Transmembrane")
    is_cazy_gt = cazy.startswith("GT") and "Tier0" not in dbcan and cazy != "-"

    support_notes = []
    if signalp_sp:
        support_notes.append("SignalP6=SP")
    if targetp_sp:
        support_notes.append("TargetP2=SP")
    if tm_signal:
        support_notes.append("DeepTMHMM=signal")
    if tm_count:
        support_notes.append(f"DeepTMHMM={tm_count}TMR")
    if extracellular:
        support_notes.append("DeepLoc2=extracellular")
    if golgi:
        support_notes.append("DeepLoc2=Golgi")
    if plastid:
        support_notes.append("DeepLoc2=plastid")
    if mitochondrion:
        support_notes.append("DeepLoc2=mitochondrion")

    if is_cazy_gt and signalp_sp and tm_signal and extracellular and (targetp_sp or dl_signal):
        tier = "A1_强分泌/胞外GT候选"
        role = "糖基转移酶候选同时具有多工具信号肽和胞外定位支持，可优先写作表面糖链/EPS或细胞壁糖基化相关的分泌型候选。"
    elif is_cazy_gt and tm_count >= 1 and (golgi or cell_membrane or dl_tm):
        tier = "A2_强膜系统/高尔基GT候选"
        role = "糖基转移酶候选同时具有跨膜或高尔基/膜系统定位线索，更适合写作分泌通路、膜系统糖基化或糖链装配相关候选。"
    elif "S0C" in s0_class and signalp_sp and tm_signal and extracellular and (targetp_sp or dl_signal):
        tier = "B1_分泌型壁糖蛋白成熟候选"
        role = "S0C 壁糖蛋白成熟相关候选具有强信号肽/胞外定位线索，可作为细胞壁或胞外基质蛋白成熟模块的关键候选。"
    elif "S0B" in s0_class and tm_count >= 2 and (golgi or cell_membrane or dl_tm):
        tier = "B2_多跨膜膜系统/前体转运候选"
        role = "S0B 候选表现为多跨膜并偏高尔基/膜系统定位，支持糖前体转运或膜系统装配辅助模块。"
    elif extracellular and (dl_signal or signalp_sp or tm_signal):
        tier = "B3_胞外/分泌通路中等候选"
        role = "有胞外或分泌通路线索，但功能证据或多工具一致性不足，适合放在候选补充层。"
    elif targetp_organelle or plastid or mitochondrion:
        tier = "C_细胞器定位需谨慎"
        role = "定位结果偏质体或线粒体，若用于表面/EPS机制写作需要降级，除非后续表达或同源证据能说明间接作用。"
    elif tm_count >= 1 or golgi or cell_membrane or dl_tm:
        tier = "C+_膜系统单线索候选"
        role = "有膜系统或跨膜线索，但功能证据不足或与 GT/EPS 主链关联较弱，可作为辅助机制候选。"
    else:
        tier = "D_定位支持弱/偏胞内"
        role = "当前定位工具不支持强分泌、跨膜或细胞器靶向，不宜作为表面/EPS主候选。"

    boundary = "所有定位均为 DTU 在线工具预测，属于生信证据；可支撑候选优先级和机制假说，不能写成已实验验证的分泌、跨膜或细胞壁定位。"
    if not support_notes:
        support_notes.append("无强定位信号")
    return tier, "；".join(support_notes), role + boundary


def load_job_meta() -> dict[str, str]:
    return {
        "DeepTMHMM": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDE000304C20E9C081DF",
        "SignalP6": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDED00304C381A6DDA90",
        "TargetP2": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDF100304C3F1BDE4213",
        "DeepLoc2": "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid=6A49CDFF00304C486001E21F",
    }


def s0_layer_en(s0_class: str) -> str:
    if "S0A" in s0_class:
        return "Direct GT"
    if "S0B" in s0_class:
        return "Precursor transport"
    if "S0C" in s0_class:
        return "Wall glycoprotein"
    if "S0D" in s0_class:
        return "Contracted glycan"
    return "Boundary"


def deeploc_group_en(localization: str) -> str:
    if localization == "Extracellular|Golgi apparatus":
        return "Extra+Golgi"
    if localization == "Mitochondrion|Golgi apparatus":
        return "Mito+Golgi"
    if "Extracellular" in localization:
        return "Extracell"
    if "Golgi" in localization:
        return "Golgi"
    if "Mitochondrion" in localization:
        return "Mito"
    if "Plastid" in localization:
        return "Plastid"
    if "Cell membrane" in localization:
        return "Cell memb."
    if "Nucleus" in localization:
        return "Nucleus"
    return "Cyt/other"


def integrated_group_en(tier: str) -> str:
    if tier.startswith("A1"):
        return "Secreted GT"
    if tier.startswith("A2"):
        return "Membrane GT"
    if tier.startswith("B1"):
        return "Secreted wall"
    if tier.startswith("B2"):
        return "Transport membrane"
    if tier.startswith("B3"):
        return "Moderate extracellular"
    if tier.startswith("C+"):
        return "Membrane single"
    if tier.startswith("C_"):
        return "Organelle caution"
    return "Weak intracellular"


def main() -> None:
    required = [
        BASE_TABLE,
        DTU_DIR / "SignalP6_prediction_results.txt",
        DTU_DIR / "TargetP2_output_protein_type.txt",
        DTU_DIR / "DeepLoc2_results_6A49CDFF00304C486001E21F.csv",
        DTU_DIR / "DeepTMHMM_TMRs.gff3",
        DTU_DIR / "DeepTMHMM_predicted_topologies.3line",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("\n".join(missing))

    base = read_tsv(BASE_TABLE)
    signalp = parse_signalp(DTU_DIR / "SignalP6_prediction_results.txt")
    targetp = parse_targetp(DTU_DIR / "TargetP2_output_protein_type.txt")
    deeploc = parse_deeploc(DTU_DIR / "DeepLoc2_results_6A49CDFF00304C486001E21F.csv")
    tmhmm = parse_deeptmhmm_gff(DTU_DIR / "DeepTMHMM_TMRs.gff3")
    tmhmm_class = parse_deeptmhmm_header_class(DTU_DIR / "DeepTMHMM_predicted_topologies.3line")

    tool_ids = set(signalp) | set(targetp) | set(deeploc) | set(tmhmm) | set(tmhmm_class)
    tool_ids_by_gene: dict[str, list[str]] = defaultdict(list)
    for pid in sorted(tool_ids):
        tool_ids_by_gene[norm_gene(pid)].append(pid)

    rows: list[dict[str, object]] = []
    for entry in base:
        gene = entry.get("基因ID", "")
        pid = choose_protein_id(gene, entry.get("蛋白ID", ""), tool_ids_by_gene)
        row: dict[str, object] = {
            "基因ID": gene,
            "蛋白ID": pid,
            "正交家族": entry.get("正交家族", ""),
            "S0类别": entry.get("S0类别", ""),
            "家族变化方向": entry.get("家族变化方向", ""),
            "CAZyme家族": entry.get("CAZyme家族", ""),
            "CAZyme类别": entry.get("CAZyme类别", ""),
            "dbCAN证据等级": entry.get("dbCAN证据等级", ""),
            "HRGP_motif判定": entry.get("HRGP_motif判定", ""),
        }
        row.update(signalp.get(pid, {}))
        row.update(targetp.get(pid, {}))
        row.update(deeploc.get(pid, {}))
        row.update(tmhmm.get(pid, {}))
        row["DeepTMHMM类别"] = tmhmm_class.get(pid, "")
        tier, support, interpretation = evidence_call(row)
        row["定位整合等级"] = tier
        row["定位支持摘要"] = support
        row["Paper1机制写作建议"] = interpretation
        rows.append(row)

    detail_fields = [
        "基因ID",
        "蛋白ID",
        "正交家族",
        "S0类别",
        "家族变化方向",
        "CAZyme家族",
        "CAZyme类别",
        "dbCAN证据等级",
        "HRGP_motif判定",
        "SignalP6预测",
        "SignalP6_SP概率",
        "SignalP6切割位点",
        "DeepTMHMM类别",
        "DeepTMHMM信号肽区间",
        "DeepTMHMM_TMR数量",
        "DeepTMHMM_TMR区间",
        "TargetP2预测",
        "TargetP2_SP概率",
        "TargetP2_mTP概率",
        "TargetP2_cTP概率",
        "TargetP2切割位点",
        "DeepLoc2定位",
        "DeepLoc2信号",
        "DeepLoc2膜类型",
        "DeepLoc2_Extracellular概率",
        "DeepLoc2_Cell membrane概率",
        "DeepLoc2_Golgi apparatus概率",
        "DeepLoc2_Endoplasmic reticulum概率",
        "DeepLoc2_Mitochondrion概率",
        "DeepLoc2_Plastid概率",
        "DeepLoc2_Transmembrane概率",
        "DeepLoc2_Soluble概率",
        "定位整合等级",
        "定位支持摘要",
        "Paper1机制写作建议",
    ]
    write_tsv(OUT_DETAIL, rows, detail_fields)

    counters = {
        "SignalP6预测": Counter(str(r.get("SignalP6预测", "未返回") or "未返回") for r in rows),
        "TargetP2预测": Counter(str(r.get("TargetP2预测", "未返回") or "未返回") for r in rows),
        "DeepTMHMM类别": Counter(str(r.get("DeepTMHMM类别", "未返回") or "未返回") for r in rows),
        "DeepLoc2定位": Counter(str(r.get("DeepLoc2定位", "未返回") or "未返回") for r in rows),
        "定位整合等级": Counter(str(r.get("定位整合等级", "未分级") or "未分级") for r in rows),
    }
    summary_rows: list[dict[str, object]] = []
    for tool, counter in counters.items():
        for label, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
            summary_rows.append({"工具或层级": tool, "类别": label, "数量": count})
    write_tsv(OUT_TOOL_SUMMARY, summary_rows, ["工具或层级", "类别", "数量"])

    figure_rows: list[dict[str, object]] = []
    for row in rows:
        try:
            tmr_n = int(row.get("DeepTMHMM_TMR数量") or 0)
        except ValueError:
            tmr_n = 0
        if tmr_n == 0:
            tmr_band = "0"
        elif tmr_n == 1:
            tmr_band = "1"
        elif tmr_n == 2:
            tmr_band = "2"
        else:
            tmr_band = "8-10"
        cazy_family = str(row.get("CAZyme家族", "") or "")
        if cazy_family in {"", "-"}:
            cazy_family = "No CAZyme"
        figure_rows.append({
            "gene_id": row.get("基因ID", ""),
            "protein_id": row.get("蛋白ID", ""),
            "s0_class": row.get("S0类别", ""),
            "s0_layer": s0_layer_en(str(row.get("S0类别", ""))),
            "cazy_family": cazy_family,
            "signalp6": row.get("SignalP6预测", ""),
            "deeptmhmm": row.get("DeepTMHMM类别", ""),
            "deeptmhmm_signal_range": row.get("DeepTMHMM信号肽区间", ""),
            "tmr_n": tmr_n,
            "tmr_band": tmr_band,
            "targetp2": row.get("TargetP2预测", ""),
            "deeploc_raw": row.get("DeepLoc2定位", ""),
            "deeploc_group": deeploc_group_en(str(row.get("DeepLoc2定位", ""))),
            "extracellular_score": row.get("DeepLoc2_Extracellular概率", ""),
            "golgi_score": row.get("DeepLoc2_Golgi apparatus概率", ""),
            "integrated_group": integrated_group_en(str(row.get("定位整合等级", ""))),
        })
    write_tsv_utf8(OUT_FIGURE_SOURCE, figure_rows, [
        "gene_id",
        "protein_id",
        "s0_class",
        "s0_layer",
        "cazy_family",
        "signalp6",
        "deeptmhmm",
        "deeptmhmm_signal_range",
        "tmr_n",
        "tmr_band",
        "targetp2",
        "deeploc_raw",
        "deeploc_group",
        "extracellular_score",
        "golgi_score",
        "integrated_group",
    ])

    by_tier: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        by_tier[str(row["定位整合等级"])].append(str(row["蛋白ID"]))

    strong_secreted = [
        r for r in rows if str(r["定位整合等级"]).startswith("A1") or str(r["定位整合等级"]).startswith("B1")
    ]
    strong_membrane = [
        r for r in rows if str(r["定位整合等级"]).startswith("A2") or str(r["定位整合等级"]).startswith("B2")
    ]
    cazy_secreted = [
        r for r in rows
        if str(r.get("CAZyme家族", "")).startswith("GT") and str(r["定位整合等级"]).startswith("A1")
    ]
    cazy_membrane = [
        r for r in rows
        if str(r.get("CAZyme家族", "")).startswith("GT") and str(r["定位整合等级"]).startswith("A2")
    ]

    jobs = load_job_meta()
    report = f"""# FACHB-416 Paper1 Round16：S0 候选 DTU 定位预测整合报告

生成时间：{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## 输入与边界

- 本轮在用户明确同意后，将 S0 小批量蛋白序列提交到 DTU 在线服务；交付文件不包含原始蛋白序列。
- 使用工具：DeepTMHMM 1.0、SignalP 6.0、TargetP 2.0（Plant 模式）和 DeepLoc 2.1。
- 本报告只把定位结果作为生信预测证据；不能写成已经实验验证的分泌、跨膜、细胞壁或细胞器定位。
- 功能层证据沿用前序 dbCAN5/Pfam/InterPro/motif 的 S0 证据表，定位层用于提升候选优先级和机制叙事。

## DTU 任务链接

- DeepTMHMM：{jobs["DeepTMHMM"]}
- SignalP6：{jobs["SignalP6"]}
- TargetP2：{jobs["TargetP2"]}
- DeepLoc2：{jobs["DeepLoc2"]}

## 核心结果

- SignalP6 预测 SP：{sum(1 for r in rows if r.get("SignalP6预测") == "SP")}/{len(rows)}。
- TargetP2 预测 SP：{sum(1 for r in rows if r.get("TargetP2预测") == "SP")}/{len(rows)}；mTP：{sum(1 for r in rows if r.get("TargetP2预测") == "mTP")}；cTP：{sum(1 for r in rows if r.get("TargetP2预测") == "cTP")}。
- DeepTMHMM 检出 signal：{sum(1 for r in rows if r.get("DeepTMHMM信号肽区间"))}/{len(rows)}；检出至少 1 个 TMR：{sum(1 for r in rows if int(r.get("DeepTMHMM_TMR数量") or 0) >= 1)}/{len(rows)}。
- DeepLoc2 含 Extracellular 定位：{sum(1 for r in rows if contains(str(r.get("DeepLoc2定位", "")), "Extracellular"))}/{len(rows)}；含 Golgi apparatus 定位：{sum(1 for r in rows if contains(str(r.get("DeepLoc2定位", "")), "Golgi"))}/{len(rows)}。

## 支持 Paper1 主线的候选层级

### A1/B1：分泌或胞外强定位线索

{", ".join(f"{r['蛋白ID']}({r.get('CAZyme家族') or r.get('S0类别')})" for r in strong_secreted) or "无"}

解释：这些候选同时获得 SignalP6、DeepTMHMM、TargetP2 或 DeepLoc2 中多条分泌/胞外证据。特别是 GT49 候选 g651.t1、g1615.t1、g5532.t1 以及 GT32 弱 dbCAN 候选 g12848.t1，可作为“扩张/保留的表面糖链或胞外基质相关模块”的优先候选；S0C 的 g12483.t1 则更适合写作细胞壁糖蛋白成熟或胞外基质蛋白加工相关候选。

### A2/B2：膜系统或高尔基强定位线索

{", ".join(f"{r['蛋白ID']}({r.get('CAZyme家族') or r.get('S0类别')})" for r in strong_membrane) or "无"}

解释：这些候选含跨膜段或 DeepLoc2 高尔基/膜系统定位，适合支撑“糖链装配、糖前体转运和分泌通路辅助模块被保留或强化”的叙事。其中 GT32/GT49 候选 g11399.t1、g9209.t1、g9850.t1 可作为膜系统糖基转移酶候选；S0B 多跨膜候选 g7484.t1、g11586.t1、g12701.t1、g12816.t1、g13195.t1 更适合作为前体转运或膜系统辅助候选。

### C 层：细胞器定位或表面机制需降级

{", ".join(pid for tier, ids in by_tier.items() if tier.startswith("C") for pid in ids) or "无"}

解释：这部分可能参与碳氮代谢、能量或细胞器状态，而不是直接定位到细胞表面。可保留为间接机制或备用候选，不建议放入 EPS/细胞壁主证据链。

## 对 Paper1 假说的增益

本轮结果把原来“S0 家族有糖基化/细胞壁相关功能”的证据推进到“其中一部分候选还具备合理的分泌、胞外、跨膜或高尔基定位预测”。这对第一篇的主线有两个作用：

1. 支撑长期实验室培养后，少数扩张或保留模块并非随机，而是集中在生长维持、膜系统糖基化、胞外/表面结构更新和形态响应相关过程。
2. 将“表面糖链/EPS/细胞壁糖蛋白”从泛泛功能注释，升级为可以按候选等级排序的机制模型：分泌型 GT/壁蛋白候选负责胞外糖链或基质输出，膜系统 GT 与多跨膜候选负责糖链装配和前体转运，细胞器定位候选作为能量或代谢背景调节层。

## 建议写法

- 可以写：S0 候选中存在一组 dbCAN 支持的 GT32/GT49/GT47/GT10 糖基转移酶候选，其中若干成员同时获得 SignalP6、DeepTMHMM、TargetP2 和 DeepLoc2 的分泌/胞外或膜系统定位预测，提示其可能参与表面糖链、EPS 或细胞壁相关糖基化过程。
- 不建议写：这些基因已经被证明定位于细胞壁或直接合成 EPS。
- 主图建议：把 Round16 表中的“功能证据等级、定位等级、S0 家族变化方向、DeepLoc 定位、TMR 数量”合成一个多轨热图或 circos 内圈，用 A1/A2/B1/B2 候选做醒目标注。

## 输出文件

- 整合明细表：`{OUT_DETAIL.name}`
- 工具分类计数：`{OUT_TOOL_SUMMARY.name}`
- R 作图源数据：`{OUT_FIGURE_SOURCE.name}`
- 本报告：`{OUT_REPORT.name}`
- 运行日志：`{OUT_LOG.name}`
"""
    OUT_REPORT.write_text(report, encoding="utf-8-sig")

    log = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "inputs": [str(p) for p in required],
        "outputs": [str(OUT_DETAIL), str(OUT_TOOL_SUMMARY), str(OUT_FIGURE_SOURCE), str(OUT_REPORT), str(OUT_LOG)],
        "row_count": len(rows),
        "tool_id_count": len(tool_ids),
        "counts": {k: dict(v) for k, v in counters.items()},
        "job_urls": jobs,
        "note": "No raw protein sequences were written to outputs by this parser.",
    }
    OUT_LOG.write_text(
        "# Round16 DTU localization parser run log\n\n"
        + "```json\n"
        + json.dumps(log, ensure_ascii=False, indent=2)
        + "\n```\n",
        encoding="utf-8-sig",
    )


if __name__ == "__main__":
    main()
