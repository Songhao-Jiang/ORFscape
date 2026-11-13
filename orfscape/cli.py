# -*- coding: utf-8 -*-
"""ORFscape 命令行入口。

本轮（R2 数据确定，馏分制语义 v1.2）：
  python -m orfscape scan --raw raw目录 --workdir 项目目录 [--no-sha256] [--skip-dump] [--offline]

分组语义：每个（物种×样品类型×酶）组 = 一个样品组；组内每个 raw = 一个馏分（按
文件名采集日期排序编号 F01..Fnn，技术重跑等同额外馏分）；镜像配对在样品级
（mirror_group_id）。多 CV（FAIMS）raw 标记 needs_cv_split=1；用户预拆分产物
（文件名 *_CV.raw）自动识别，其多 CV 原 raw 标记 superseded 跳过。

产物（项目目录 00_input/）：groups.csv / species.tsv / raw_info.tsv /
methods/ / scan_tables/，及 run_state.json 阶段状态机。
"""
import argparse
import csv
import datetime
import os
import sys

from . import __version__
from .core import inputs, rawinfo, taxonomy
from .core.validate import DEFAULT_ISOPSV as V_ISOPSV_DEFAULT


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _log(msg):
    """输出日志：经 i18n 翻译（默认英文，ORFSCAPE_LANG=zh 或 --lang zh 切中文）。

    同时镜像到结构化事件流（``ORFSCAPE_EVENTS``），供 GUI 读取；未设置时无副作用。
    """
    from .core.i18n import L
    from . import progress
    text = f"[{_now()}] {L(msg)}"
    print(text, flush=True)
    progress.log(text)


# ---------------------------------------------------------------- 项目脚手架（字面量子目录）
def ensure_project(workdir):
    """校验并创建项目目录与 12 个字面量子目录；返回项目根路径字符串。

    路径安全：abspath 规范化并拒绝 ".." 成分。
    """
    root = os.path.abspath(str(workdir))
    if ".." in root.replace("/", os.sep).split(os.sep):
        raise ValueError(f"项目路径含穿越成分: {root}")
    os.makedirs(root, exist_ok=True)
    os.makedirs(os.path.join(root, "00_input"), exist_ok=True)
    os.makedirs(os.path.join(root, "01_genome"), exist_ok=True)
    os.makedirs(os.path.join(root, "02_db"), exist_ok=True)
    os.makedirs(os.path.join(root, "02_tx"), exist_ok=True)
    os.makedirs(os.path.join(root, "03_mgf"), exist_ok=True)
    os.makedirs(os.path.join(root, "04_search_annot"), exist_ok=True)
    os.makedirs(os.path.join(root, "04a_qc"), exist_ok=True)
    os.makedirs(os.path.join(root, "04_search_disc"), exist_ok=True)
    os.makedirs(os.path.join(root, "05_novel"), exist_ok=True)
    os.makedirs(os.path.join(root, "06_arbitration"), exist_ok=True)
    os.makedirs(os.path.join(root, "07_annotation"), exist_ok=True)
    os.makedirs(os.path.join(root, "08_report"), exist_ok=True)
    return root


def load_state(root):
    import json
    p = os.path.join(root, "run_state.json")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    return {"project": os.path.basename(root), "created": _now(), "stages": {}}


def save_state(root, state):
    import json
    state["updated"] = _now()
    with open(os.path.join(root, "run_state.json"), "w", encoding="utf-8") as fh:
        json.dump(state, fh, ensure_ascii=False, indent=2)


# ---------------------------------------------------------------- 主流程：scan
def cmd_scan(args):
    _log(f"ORFscape v{__version__} · 数据确定开始")
    # 1) 发现与分组推断（馏分制）
    raws = inputs.discover_raws(args.raw)
    if not raws:
        _log(f"未发现 .raw 文件: {args.raw}")
        return 2
    _log(f"发现 {len(raws)} 个 raw 文件")
    rows, unmatched, mirror_mode = inputs.build_groups(raws)
    if unmatched:
        _log(f"⚠ {len(unmatched)} 个文件未能推断分组（详见 raw_index.unmatched）")
        for p in unmatched[:10]:
            _log(f"   未匹配: {p.name}")
    _log(f"分组推断完成（馏分制）：{'镜像酶切模式' if mirror_mode else '单酶模式'}，"
         f"{len({r['mirror_group_id'] for r in rows if r['mirror_group_id']})} 个样品级镜像组")
    n_pre = sum(1 for r in rows if r.get("pre_split_cv"))
    n_sup = sum(1 for r in rows if r.get("superseded") == "yes")
    if n_pre:
        _log(f"检测到 {n_pre} 个预拆分电压 raw（*_CV.raw）；{n_sup} 个多 CV 原 raw 标记跳过")

    # 2) 物种解析（--species-code / --species-name 可强制覆盖文件名推断）
    if getattr(args, "species_code", None):
        for r in rows:
            r["species_code"] = args.species_code
        _log(f"物种代码覆盖为 {args.species_code}")
    sp_codes = sorted({r["species_code"] for r in rows})
    sp_res = {c: taxonomy.resolve_species(c, online=not args.offline) for c in sp_codes}
    if getattr(args, "species_name", None):
        _name = str(args.species_name).strip()
        try:
            _res = taxonomy.resolve_species(_name, online=not args.offline)
        except Exception as exc:                       # 离线/无命中都不阻断
            _log(f"物种名解析失败，按给定名记录: {exc}")
            _res = {}
        _res = dict(_res or {})
        _res["species"] = _name
        _res["taxid"] = _res.get("taxid", "")
        _res["matched"] = "manual"
        for c in sp_codes:
            sp_res[c] = dict(_res)
        _log(f"物种名覆盖为 {_name}")
    for c, s in sp_res.items():
        _log(f"物种 {c} → {s['species']} (taxid={s['taxid'] or '待定'}, {s['matched']})"
             + (f" [{s['note']}]" if s["note"] else ""))
    for r in rows:
        r["species"] = sp_res[r["species_code"]]["species"]
        r["taxid"] = sp_res[r["species_code"]]["taxid"]

    # 3) 项目目录与产物
    root = ensure_project(args.workdir)
    inp = os.path.join(root, "00_input")
    inputs.save_groups_csv(rows, os.path.join(inp, "groups.csv"))
    with open(os.path.join(inp, "species.tsv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species_code", "species", "taxid", "matched", "hit_term", "note"])
        for c, s in sp_res.items():
            w.writerow([c, s["species"], s["taxid"], s["matched"], s.get("hit_term", ""), s["note"]])

    # 4) raw 信息：rawmethod（方法参数）+ rawdump（扫描表）——按源目录分组 staging
    by_dir = {}
    for r in rows:
        by_dir.setdefault(r["group_dir"], []).append(r)
    method_all, scan_all = {}, {}
    methods_dir = os.path.join(inp, "methods")
    scans_dir = os.path.join(inp, "scan_tables")
    os.makedirs(methods_dir, exist_ok=True)
    os.makedirs(scans_dir, exist_ok=True)
    for i, (gdir, members) in enumerate(sorted(by_dir.items()), 1):
        gid = f"{members[0]['species_code']}_{members[0]['sample_type']}_{members[0]['enzyme']}"
        _log(f"[{i}/{len(by_dir)}] {gid}: {len(members)} raw ← {gdir}")
        staging = os.path.join(root, "_staging", gid)
        rawinfo.stage_group([m["path"] for m in members], staging)
        try:
            if not args.skip_method:
                mt = os.path.join(methods_dir, f"method_{gid}.tsv")
                ok, _ = rawinfo.run_rawmethod(staging, mt)
                if ok:
                    method_all.update(rawinfo.parse_method_tsv(mt))
                    _log(f"   rawmethod ✓ ({len(method_all)} raw 方法已解析)")
                else:
                    _log("   rawmethod ✗（跳过方法解析）")
            if not args.skip_dump:
                ok, tsv = rawinfo.run_rawdump(staging)
                if ok:
                    dst = os.path.join(scans_dir, f"scan_{gid}.tsv")
                    os.replace(tsv, dst)
                    scan_all.update(rawinfo.parse_scan_tsv(dst))
                    _log(f"   rawdump ✓ ({len(scan_all)} raw 扫描表已解析)")
                else:
                    _log("   rawdump ✗")
        finally:
            rawinfo.clean_stage(staging)

    # 5) sha256 与文件信息
    if not args.no_sha256:
        _log("计算 sha256（大文件较慢，可用 --no-sha256 跳过）…")
        for i, r in enumerate(rows, 1):
            r["sha256"] = rawinfo.sha256sum(r["path"])
            if i % 8 == 0:
                _log(f"   sha256 {i}/{len(rows)}")
    for r in rows:
        fi = rawinfo.file_info(r["path"])
        r.update(fi)

    # 6) 汇总 raw_info.tsv
    seg_summary = lambda segs: "; ".join(
        f"CV{s['cv']}:res{s['res1']}/{s['res2']}:IT{s['ms1_it']}/{s['ms2_it']}:"
        f"cyc{s['cycle']}:nce{s['nce']}" for s in (segs or []))
    # 输出目录校验（防穿越）：inp 必须是已验证项目根下的字面量子目录
    inp = os.path.abspath(inp)
    if ".." in inp.replace("/", os.sep).split(os.sep):
        raise ValueError(f"输出目录含穿越成分: {inp}")
    if os.path.commonpath([inp, root]) != root:
        raise ValueError(f"输出目录越界: {inp}")
    fields = ["raw_file", "species_code", "species", "taxid", "sample_type", "condition",
              "enzyme", "enzyme_name", "cleavage", "fraction", "n_fraction_in_group",
              "mirror_group_id", "pattern", "pre_split_cv", "superseded",
              "size_mb", "mtime", "sha256", "n_scan", "n_ms1", "n_ms2",
              "rt_min_first", "rt_min_last", "cvs", "n_cv", "needs_cv_split",
              "resolutions", "cycle_times", "nce", "method_segments", "agc_pcts",
              "med_it_ms2", "it_sat_pct", "med_fill_ms1", "med_fill_ms2"]
    with open(os.path.join(inp, "raw_info.tsv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            mi = method_all.get(r["raw_file"], {})
            sc = scan_all.get(r["raw_file"], {})
            cvs_raw = sc.get("cvs", "")
            n_cv = len([x for x in cvs_raw.split(";") if x]) if cvs_raw else ""
            r2 = dict(r)
            r2.update({
                "cvs": cvs_raw, "n_cv": n_cv,
                "needs_cv_split": ("1" if (n_cv and int(n_cv) > 1 and not r.get("pre_split_cv")) else "0"),
                "n_scan": sc.get("n_scan", ""),
                "n_ms1": sc.get("n_ms1", ""), "n_ms2": sc.get("n_ms2", ""),
                "rt_min_first": sc.get("rt_first", ""), "rt_min_last": sc.get("rt_last", ""),
                "resolutions": mi.get("resolutions", ""), "cycle_times": mi.get("cycle_times", ""),
                "nce": mi.get("nce", ""), "method_segments": seg_summary(mi.get("segments")),
                "agc_pcts": ";".join(mi.get("agc_pcts", [])),
                "med_it_ms2": sc.get("med_it_ms2", ""),
                "it_sat_pct": rawinfo.saturation_pct(sc, mi),
                "med_fill_ms1": sc.get("med_fill_ms1", ""), "med_fill_ms2": sc.get("med_fill_ms2", ""),
                "sha256": r.get("sha256", ""),
            })
            w.writerow(r2)
    if unmatched:
        with open(os.path.join(inp, "raw_index.unmatched.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(str(p) for p in unmatched))

    # 7) 状态机与摘要
    state = load_state(root)
    state.setdefault("stages", {})["data_scan"] = {
        "done": True, "ts": _now(), "n_raw": len(rows), "n_unmatched": len(unmatched),
        "mirror_mode": mirror_mode,
        "n_mirror_groups": len({r["mirror_group_id"] for r in rows if r["mirror_group_id"]}),
        "n_pre_split": n_pre, "n_superseded": n_sup,
        "species": {c: sp_res[c]["species"] for c in sp_codes},
        "orfscape_version": __version__,
    }
    save_state(root, state)
    _log("──────── 摘要 ────────")
    for gid_key in sorted({(r['species_code'], r['sample_type'], r['enzyme']) for r in rows}):
        members = [r for r in rows if (r['species_code'], r['sample_type'], r['enzyme']) == gid_key]
        scans = sum(scan_all.get(m["raw_file"], {}).get("n_scan", 0) for m in members)
        mg = members[0]["mirror_group_id"]
        _log(f"{gid_key[0]}/{gid_key[1]}/{gid_key[2]} [{mg}]: {len(members)} 馏分, "
             f"共 {scans} 扫描")
    _log(f"完成。产物目录: {inp}")
    return 0


def cmd_doctor(args):
    """环境检测器（A 阶段）：工具/引擎自动发现，产出 env_report。"""
    from .core import envdetect
    _log(f"ORFscape v{__version__} · 环境检测开始")
    rows = envdetect.detect_all()
    n_ok, n_all, req_missing, rec_missing = envdetect.summary(rows)
    root = ensure_project(args.workdir)
    out = os.path.join(root, "env_report.tsv")
    envdetect.write_report(rows, out)
    state = load_state(root)
    state.setdefault("stages", {})["env_detect"] = {
        "done": True, "ts": _now(), "n_ok": n_ok, "n_all": n_all,
        "required_missing": req_missing, "recommended_missing": rec_missing,
        "tools": {r["name"]: {"status": r["status"], "version": r["version"]} for r in rows},
    }
    save_state(root, state)
    for r in rows:
        mark = "✓" if r["status"] == "ok" else "✗"
        ver = f" v{r['version']}" if r["version"] else ""
        _log(f"  {mark} [{r['need']}] {r['name']}{ver} — {r['category']}"
             + (f" | {r['note']}" if r["status"] == "missing" and r["note"] else ""))
    _log(f"检测完成：{n_ok}/{n_all} 可用" +
         (f"；⚠ 必需缺失 {req_missing}（阻断后续）" if req_missing else "") +
         (f"；建议缺失 {rec_missing}（按降级链继续）" if rec_missing else ""))
    _log(f"报告: {out}")
    return 0 if not req_missing else 2


def cmd_genome(args):
    """基因组获取（B 阶段）：按项目 species.tsv 逐物种下载三件套并钉扎。"""
    from .core import genome as gmod
    _log(f"ORFscape v{__version__} · 基因组获取开始")
    root = ensure_project(args.workdir)
    inp = os.path.join(root, "00_input")
    sp_tsv = os.path.join(inp, "species.tsv")
    if not os.path.exists(sp_tsv):
        _log(f"缺 {sp_tsv}（先运行 scan）")
        return 2
    import csv as _csv
    metas = {}
    with open(sp_tsv, newline="", encoding="utf-8-sig") as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            code = r.get("species_code", "")
            sp = r.get("species", "")
            if not (code and sp):
                continue
            out_dir = os.path.join(root, "01_genome", code)
            _log(f"[{code}] {sp} → 获取参考装配…")
            meta = gmod.fetch_genome(sp, out_dir, today=_now()[:10])
            metas[code] = meta
            _log(f"  ✓ {meta['accession']}（{meta['refseq_category']}）"
                 f" 蛋白 {meta.get('n_proteins')} 条")
    state = load_state(root)
    state.setdefault("stages", {})["genome_fetch"] = {
        "done": True, "ts": _now(),
        "genomes": {c: {"accession": m["accession"], "category": m["refseq_category"],
                        "n_proteins": m.get("n_proteins")} for c, m in metas.items()}}
    save_state(root, state)
    _log("基因组获取完成")
    return 0


def cmd_split(args):
    """D0 电压拆分全量执行：needs_cv_split=1 的 raw → 分 CV MGF + rawdump 对账。"""
    from .core import cvsplit, envdetect
    _log(f"ORFscape v{__version__} · 电压拆分开始（D0）")
    root = ensure_project(args.workdir)
    inp = os.path.join(root, "00_input")
    raw_info = os.path.join(inp, "raw_info.tsv")
    if not os.path.exists(raw_info):
        _log("缺 raw_info.tsv（先运行 scan）")
        return 2
    msrow = next((r for r in envdetect.detect_all() if r["name"].startswith("msconvert")), None)
    if not msrow or msrow["status"] != "ok":
        _log("✗ 未找到 msconvert（ProteoWizard）")
        return 2
    msc = msrow["path"].split(";")[0].strip()
    _log(f"msconvert: {msc}")
    rows = cvsplit.build_rows(raw_info)
    # path 列在 groups.csv（raw_info 不含）：按 raw_file 关联补上
    groups_csv = os.path.join(inp, "groups.csv")
    import csv as _csv2
    path_of = {}
    with open(groups_csv, newline="", encoding="utf-8-sig") as fh:
        for g in _csv2.DictReader(fh):
            path_of[g.get("raw_file", "")] = g.get("path", "")
    for r in rows:
        r["path"] = path_of.get(r.get("raw_file", ""), "")
    rows = [r for r in rows if r.get("path")]
    _log(f"待拆分 raw：{len(rows)} 个（每 raw × 2 CV）")
    # rawdump 期望计数：{(raw文件名, cv): MS2 数}
    expect = {}
    scandir = os.path.join(inp, "scan_tables")
    for t in os.listdir(scandir) if os.path.isdir(scandir) else []:
        if t.endswith(".tsv"):
            expect.update(cvsplit.scan_ms2_by_cv(os.path.join(scandir, t)))
    _log(f"对账基准：{len(expect)} 个 (raw×CV) 期望计数")
    mgf_dir = os.path.join(root, "03_mgf")
    # ---- 预拆分数据（一个 raw 一个 CV，如 TP）：不走 scanEvent 拆分，直接转 MGF ----
    if not rows:
        pre = []
        with open(raw_info, newline="", encoding="utf-8-sig") as fh:
            for r in _csv2.DictReader(fh, delimiter="\t"):
                if r.get("pre_split_cv") and r.get("superseded", "") != "yes":
                    r["path"] = path_of.get(r.get("raw_file", ""), "")
                    if r.get("path"):
                        pre.append(r)
        if pre:
            _log(f"预拆分 raw：{len(pre)} 个（每 raw 单 CV，直接转 MGF）")
            pre_report = []
            for i, r in enumerate(pre, 1):
                raw_path = os.path.abspath(r["path"])
                if ".." in raw_path.replace("/", os.sep).split(os.sep):
                    _log(f"  ✗ 跳过含穿越成分路径: {r.get('raw_file')}")
                    continue
                if not os.path.isfile(raw_path):
                    _log(f"  ✗ 文件不存在: {r.get('raw_file')}")
                    continue
                fname = os.path.basename(raw_path)
                if not fname.lower().endswith(".raw"):
                    continue
                stem_raw = fname[:-4]
                cv_tag = str(r.get("pre_split_cv", "")).strip()
                if stem_raw.endswith("_" + cv_tag):
                    base = stem_raw[:-(len(cv_tag) + 1)]
                else:
                    base = stem_raw
                out_mgf = os.path.join(mgf_dir, f"{base}_{cv_tag}.mgf")
                if os.path.exists(out_mgf) and cvsplit._count_ions(out_mgf) > 0:
                    st, n = "REUSE", cvsplit._count_ions(out_mgf)
                else:
                    st, n = cvsplit.convert_pre_split(raw_path, out_mgf, cv_tag,
                                                      msc, base)
                pre_report.append((f"{base}_{cv_tag}.mgf", r.get("raw_file", ""),
                                   cv_tag, st, n))
                _log(f"  [{i}/{len(pre)}] {base}_{cv_tag} {st} 谱={n}")
            rep_path = os.path.join(mgf_dir, "split_report.tsv")
            with open(rep_path, "w", newline="", encoding="utf-8-sig") as fh:
                w = _csv2.writer(fh, delimiter="\t")
                w.writerow(["mgf", "raw_file", "cv", "status", "n_spectra"])
                w.writerows(pre_report)
            ok = sum(1 for x in pre_report if x[3] in ("OK", "REUSE"))
            _log(f"预拆分转换完成：{ok}/{len(pre_report)} OK → {rep_path}")
            return 0
    report = []
    for i, r in enumerate(rows, 1):
        # 路径安全：abspath + 拒绝 ".." + 必须真实存在；stem 取 basename 派生
        raw_path = os.path.abspath(r.get("path", ""))
        if ".." in raw_path.replace("/", os.sep).split(os.sep):
            _log(f"  ✗ 跳过含穿越成分路径: {r.get('raw_file')}")
            continue
        if not os.path.isfile(raw_path):
            _log(f"  ✗ 文件不存在: {r.get('raw_file')}")
            continue
        raw_name = os.path.basename(raw_path)
        stem = raw_name[:-4]
        for cv_tag in ("45", "65"):
            ev = 2 if cv_tag == "45" else 4
            out_name = stem + "_" + cv_tag + ".mgf"
            out_mgf = os.path.join(mgf_dir, out_name)
            if os.path.exists(out_mgf):
                n = cvsplit._count_ions(out_mgf)
                status = "OK" if cvsplit.reconcile(n, expect, raw_name, cvsplit.cv_from_tag(cv_tag)) else "MISMATCH"
                report.append((raw_name, cv_tag, status, n, "复用"))
                continue
            status, n = cvsplit.split_one(raw_path, out_mgf, cv_tag, ev, msc, stem)
            if status == "OK" and not cvsplit.reconcile(n, expect, raw_name, cvsplit.cv_from_tag(cv_tag)):
                status = "MISMATCH"
            report.append((raw_name, cv_tag, status, n, ""))
        if i % 4 == 0 or i == len(rows):
            _log(f"  [{i}/{len(rows)}] …{raw_name[-38:]} 完成")
    out_tsv = os.path.join(mgf_dir, "split_report.tsv")
    import csv as _csv
    with open(out_tsv, "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["raw_file", "cv", "status", "n_spectra", "note"])
        w.writerows(report)
    n_ok = sum(1 for x in report if x[2] == "OK")
    n_bad = len(report) - n_ok
    state = load_state(root)
    state.setdefault("stages", {})["cv_split"] = {
        "done": n_bad == 0, "ts": _now(), "n_mgf": len(report),
        "n_ok": n_ok, "n_failed": n_bad}
    save_state(root, state)
    _log(f"拆分完成：{n_ok}/{len(report)} 对账通过"
         + (f"；⚠ {n_bad} 失配（见 split_report.tsv）" if n_bad else ""))
    return 0 if n_bad == 0 else 2


def cmd_dbbuild(args):
    """C 建库：按项目 01_genome 逐物种构建 L0-L5 层与两段漏斗拆分库（路径守卫版）。"""
    from .core import dbbuild
    _log(f"ORFscape v{__version__} . 建库开始（C 阶段）")
    root = ensure_project(args.workdir)
    # 路径安全：基因组根与输出根均在已验证项目根之下（字面量子目录名）
    gdir = os.path.join(root, "01_genome")
    dbroot = os.path.join(root, "02_db")
    for pth in (gdir, dbroot):
        ap = os.path.abspath(pth)
        if ".." in ap.replace("/", os.sep).split(os.sep):
            raise ValueError(f"路径含穿越成分: {ap}")
        if os.path.commonpath([ap, root]) != root:
            raise ValueError(f"路径越界: {ap}")
    if not os.path.isdir(gdir):
        _log("缺 01_genome（先运行 genome）")
        return 2
    all_stats = []
    for code in sorted(os.listdir(gdir)):
        sp_dir = os.path.join(gdir, code)
        out_dir = os.path.join(dbroot, code)
        if not os.path.isdir(sp_dir):
            continue
        _log(f"[{code}] 构建层级库…")
        st = dbbuild.build_species(code, sp_dir, out_dir, log=_log)
        all_stats.append(st)
        _log(f"  L0={st['L0']} L2={st['L2']} L3={st['L3']} | "
             f"db_annot={st['db_annot']} db_aug={st['db_aug']} db_6f={st['db_6f']}")
        _log(f"  声明: {st['notes']}")
    import csv as _csv
    out = os.path.join(dbroot, "db_stats.tsv")
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(all_stats[0].keys()), delimiter="	")
        w.writeheader()
        w.writerows(all_stats)
    state = load_state(root)
    state.setdefault("stages", {})["db_build"] = {
        "done": True, "ts": _now(),
        "stats": {s["code"]: {k: v for k, v in s.items() if k != "code"} for s in all_stats}}
    save_state(root, state)
    _log(f"建库完成。统计表: {out}")
    return 0


def cmd_e1(args):
    """E1 已注释库搜库（短任务名防 MAX_PATH：物种_样品_酶_馏分_CV）。"""
    from .core import envdetect, pfsearch
    from . import progress
    BS = chr(92)  # 反斜杠
    _log(f"ORFscape v{__version__} . E1 已注释库搜库开始")
    root = ensure_project(args.workdir)
    inp = os.path.join(root, "00_input")
    mgf_dir = os.path.join(root, "03_mgf")
    out_root = os.path.join(root, "04_search_annot")
    prow = next((r for r in envdetect.detect_all() if r["name"] == "pFind"), None)
    if not prow or prow["status"] != "ok":
        _log("pFind 未检测到")
        return 2
    pfind_bin = os.path.dirname(prow["path"].split(";")[0].strip())
    _log(f"pFind bin: {pfind_bin}")
    import csv as _csv
    # raw_file -> (species, sample_type, enzyme, fraction) 与 raw_stem -> 短名
    info = {}
    with open(os.path.join(inp, "groups.csv"), newline="", encoding="utf-8-sig") as fh:
        for g in _csv.DictReader(fh):
            info[g.get("raw_file", "")] = (g.get("species_code", ""), g.get("sample_type", ""),
                                           g.get("enzyme", ""), g.get("fraction", ""))
    tasks = []
    for f in sorted(os.listdir(mgf_dir)):
        if not f.endswith(".mgf"):
            continue
        base = f[:-4]
        raw_name, cv = base.rsplit("_", 1)
        raw_name += ".raw"
        sp, st, enz, frac = info.get(raw_name, ("", "", "", ""))
        if not sp:
            # 预拆分数据（raw 名自带 CV，如 ..._T1_45.raw）回退：整段 stem 当 raw 名
            alt = base + ".raw"
            sp, st, enz, frac = info.get(alt, ("", "", "", ""))
            if sp:
                raw_name = alt
        if not sp:
            continue
        short = f"{sp}_{st}_{enz}_{frac}_{cv}"
        tasks.append((f, raw_name, cv, sp, enz, short))
    _log(f"任务数: {len(tasks)}（raw x CV，短任务名）")
    report = []
    for i, (f, raw_name, cv, sp, enz, short) in enumerate(tasks, 1):
        tdir = os.path.join(os.path.join(out_root, sp), short)
        os.makedirs(tdir, exist_ok=True)
        spectra = os.path.join(tdir, "pFind-Filtered.spectra")
        if os.path.exists(spectra):
            n = pfsearch.count_result_rows(spectra)
            report.append((short, sp, enz, cv, "OK", n, "reuse"))
            continue
        fasta = os.path.join(os.path.join(root, "02_db", sp), "db_annot.faa")
        cfg_text = pfsearch.build_cfg_lines(
            fasta.replace(BS, "/"), os.path.join(mgf_dir, f).replace(BS, "/"),
            tdir.replace(BS, "/"), short, enz, thread=args.thread)
        cfg_path = os.path.join(tdir, short + ".cfg")
        with open(cfg_path, "w", encoding="gbk", newline=chr(10)) as fh:
            fh.write(cfg_text)
        rc = pfsearch.run_searcher(cfg_path, pfind_bin)
        n = pfsearch.count_result_rows(spectra)
        status = "OK" if rc == 0 and n > 0 else f"FAIL(rc={rc})"
        report.append((short, sp, enz, cv, status, n, ""))
        progress.tick(i, len(tasks), "task", name=short, status=status, psm=n)
        if i % 4 == 0 or i == len(tasks):
            _log(f"  [{i}/{len(tasks)}] {short} {status} PSM={n}")
    out = os.path.join(out_root, "e1_report.tsv")
    import csv as _csv3
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv3.writer(fh, delimiter=chr(9))
        w.writerow(["task", "species", "enzyme", "cv", "status", "n_psm", "note"])
        w.writerows(report)
    n_ok = sum(1 for r in report if r[4] == "OK")
    state = load_state(root)
    state.setdefault("stages", {})["e1_search"] = {
        "done": n_ok == len(report), "ts": _now(), "n_tasks": len(report), "n_ok": n_ok,
        "total_psm": sum(r[5] for r in report)}
    save_state(root, state)
    _log(f"E1 完成: {n_ok}/{len(report)} 任务，总 PSM {sum(r[5] for r in report)}")
    return 0 if n_ok == len(report) else 2


def cmd_e2(args):
    """E2 发现搜库：db_aug / db_6f（层信息在蛋白表头前缀，事后归层+分开FDR）。"""
    from .core import envdetect, pfsearch
    from . import progress
    BS2 = chr(92)
    _log(f"ORFscape v{__version__} . E2 发现搜库开始（库: {args.dbs}）")
    root = ensure_project(args.workdir)
    inp = os.path.join(root, "00_input")
    mgf_dir = os.path.join(root, "03_mgf")
    out_root = os.path.join(root, "04_search_disc")
    prow = next((r for r in envdetect.detect_all() if r["name"] == "pFind"), None)
    if not prow or prow["status"] != "ok":
        _log("pFind 未检测到")
        return 2
    pfind_bin = os.path.dirname(prow["path"].split(";")[0].strip())
    _log(f"pFind bin: {pfind_bin}")
    import csv as _csv
    info = {}
    with open(os.path.join(inp, "groups.csv"), newline="", encoding="utf-8-sig") as fh:
        for g in _csv.DictReader(fh):
            info[g.get("raw_file", "")] = (g.get("species_code", ""), g.get("sample_type", ""),
                                           g.get("enzyme", ""), g.get("fraction", ""))
    dbs = [d.strip() for d in args.dbs.split(",") if d.strip()]
    tasks = []
    for f in sorted(os.listdir(mgf_dir)):
        if not f.endswith(".mgf"):
            continue
        base = f[:-4]
        raw_name, cv = base.rsplit("_", 1)
        raw_name += ".raw"
        sp, st, enz, frac = info.get(raw_name, ("", "", "", ""))
        if not sp:
            # 预拆分数据（raw 名自带 CV，如 ..._T1_45.raw）回退：整段 stem 当 raw 名
            alt = base + ".raw"
            sp, st, enz, frac = info.get(alt, ("", "", "", ""))
            if sp:
                raw_name = alt
        if not sp:
            continue
        short = f"{sp}_{st}_{enz}_{frac}_{cv}"
        for d in dbs:
            tasks.append((f, sp, enz, cv, short, d))
    _log(f"任务数: {len(tasks)}（raw x CV x 库）")
    report = []
    for i, (f, sp, enz, cv, short, dbk) in enumerate(tasks, 1):
        tname = short + "__" + dbk
        tdir = os.path.join(os.path.join(out_root, sp), tname)
        os.makedirs(tdir, exist_ok=True)
        spectra = os.path.join(tdir, "pFind-Filtered.spectra")
        if os.path.exists(spectra):
            n = pfsearch.count_result_rows(spectra)
            report.append((tname, sp, enz, cv, dbk, "OK", n, "reuse"))
            continue
        fasta = os.path.join(os.path.join(os.path.join(root, "02_db", sp), "db_" + dbk + ".faa"), "")
        fasta = fasta[:-1]
        cfg_text = pfsearch.build_cfg_lines(
            fasta.replace(BS2, "/"), os.path.join(mgf_dir, f).replace(BS2, "/"),
            tdir.replace(BS2, "/"), tname, enz, thread=args.thread)
        cfg_path = os.path.join(tdir, tname + ".cfg")
        with open(cfg_path, "w", encoding="gbk", newline=chr(10)) as fh:
            fh.write(cfg_text)
        rc = pfsearch.run_searcher(cfg_path, pfind_bin)
        n = pfsearch.count_result_rows(spectra)
        status = "OK" if rc == 0 and n > 0 else f"FAIL(rc={rc})"
        report.append((tname, sp, enz, cv, dbk, status, n, ""))
        progress.tick(i, len(tasks), "task", name=tname, status=status, psm=n)
        if i % 8 == 0 or i == len(tasks):
            _log(f"  [{i}/{len(tasks)}] {tname} {status} PSM={n}")
    out = os.path.join(out_root, "e2_report.tsv")
    import csv as _csv3
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = _csv3.writer(fh, delimiter=chr(9))
        w.writerow(["task", "species", "enzyme", "cv", "db", "status", "n_psm", "note"])
        w.writerows(report)
    n_ok = sum(1 for r in report if r[5] == "OK")
    state = load_state(root)
    state.setdefault("stages", {})["e2_search"] = {
        "done": n_ok == len(report), "ts": _now(), "n_tasks": len(report), "n_ok": n_ok,
        "dbs": dbs, "total_psm": sum(r[6] for r in report)}
    save_state(root, state)
    _log(f"E2 完成: {n_ok}/{len(report)} 任务，总 PSM {sum(r[6] for r in report)}")
    return 0 if n_ok == len(report) else 2


def cmd_qc(args):
    """QC（紧跟 E1）：逐任务酶切质控 + 方法/扫描警告 → 04a_qc/。"""
    from .core import qc
    _log(f"ORFscape v{__version__} . QC 开始（E1 后）")
    root = ensure_project(args.workdir)
    annot = os.path.join(root, "04_search_annot")
    qcdir = os.path.join(root, "04a_qc")
    os.makedirs(qcdir, exist_ok=True)
    warns = qc.build_qc(os.path.join(root, "00_input", "raw_info.tsv"), annot,
                        os.path.join(qcdir, "qc_per_task.tsv"),
                        os.path.join(qcdir, "qc_summary.md"))
    for wn in warns[:15]:
        _log(f"  ⚠ {wn}")
    state = load_state(root)
    state.setdefault("stages", {})["qc"] = {"done": True, "ts": _now(), "n_warn": len(warns)}
    save_state(root, state)
    _log(f"QC 完成：{len(warns)} 条警告。产物 04a_qc/")
    return 0


def cmd_validate(args):
    """R8 收口：novel PSM 全量过 isoPSV，Tier 分流门控 → 06_validate/。"""
    from .core import validate as V
    _log(f"ORFscape v{__version__} . R8 novel 收口（isoPSV）开始")
    root = ensure_project(args.workdir)
    subsets = None
    if args.smoke:
        sps = inputs.project_species(root)
        if not sps:
            _log("  --smoke 需要项目内已有物种（先跑 P6）")
            return 2
        subsets = [(sps[0], "TL")]
    source = "research" if getattr(args, "from_research", False) else "novel"
    out_name = "06b_validate_research" if source == "research" else "06_validate"
    if getattr(args, "out", None):
        out_name = str(args.out)
        if out_name not in ("06_validate", "06b_validate_research", "06c_validate_full"):
            _log(f"非法输出目录：{out_name}（只允许 06_validate / 06b_validate_research / 06c_validate_full）")
            return 2
    _log(f"  输出目录：{out_name}")
    reuse, only, ms, mc = None, None, args.min_score, args.min_cos
    if getattr(args, "reuse", None):
        reuse = V.load_reuse_verdicts(root, args.reuse)
        merged = set()
        for sp in inputs.project_species(root):
            for r in V.read_table_in(root, "05b_research", sp,
                                     "novel_after_research.tsv"):
                merged.add(r["peptide"])
        only = merged - set(reuse)
        _log(f"  增量模式：合并集 {len(merged)} 肽 − 已验 {len(reuse)} = 本轮需验 {len(only)}")
        if getattr(args, "cards", None):
            parts = str(args.cards).replace("，", ",").split(",")
            if len(parts) == 2:
                ms, mc = float(parts[0]), float(parts[1])
                _log(f"  沿用卡值 score≥{ms} cos≥{mc}")
    verd, _ps, _tm, cards = V.build_validate(
        root, out_name=out_name, iso_root=args.isopsv, min_score=ms,
        min_cos=mc, subsets=subsets, source=source,
        only_peptides=only, reuse=reuse)
    out_tsv, sumf = V.write_report(root, verd, cards, out_name=out_name,
                                   source=source)
    state = load_state(root)
    state.setdefault("stages", {})["validate"] = {
        "done": True, "ts": _now(), "n_peptides": len(verd),
        "min_score": cards[0], "min_cos": cards[1]}
    save_state(root, state)
    _log(f"R8 完成：{len(verd)} 肽验证，卡值 score≥{cards[0]} cos≥{cards[1]}。"
         f"产物 {out_tsv}")
    return 0


def cmd_novel(args):
    """P6 新肽筛选：E2 层命中 → 四滤 → 坐标定位 → Tier → 05_novel/<sp>/。"""
    from .core import novel as N
    _log(f"ORFscape v{__version__} . P6 新肽筛选")
    root = ensure_project(args.workdir)
    sps = [s.strip() for s in str(args.species).split(",") if s.strip()]
    if not sps:
        sps = list_species_novel(root)
    if not sps:
        _log("05b_research/ 或 05_novel/ 下没有可识别的物种目录，用 --species 指定")
        return 2
    stats = []
    for sp in sps:
        st = N.build_species_novel(
            sp,
            os.path.join(root, "02_db", sp),
            os.path.join(root, "04_search_annot"),
            os.path.join(root, "04_search_disc"),
            os.path.join(root, "05_novel", sp),
            log=_log)
        stats.append(st)
    tot = sum(s.get("n_novel", 0) for s in stats)
    tx = sum(s.get("src_tx", 0) for s in stats)
    l3 = sum(s.get("src_altstart", 0) for s in stats)
    l1 = sum(s.get("src_abinitio", 0) for s in stats)
    l2 = tot - tx - l3 - l1
    t1 = sum(s.get("tier1", 0) for s in stats)
    t2 = sum(s.get("tier2", 0) for s in stats)
    t3 = sum(s.get("tier3", 0) for s in stats)
    state = load_state(root)
    state.setdefault("stages", {})["novel"] = {
        "done": True, "ts": _now(), "n_species": len(stats),
        "n_peptide": tot, "n_tx": tx, "n_l3": l3}
    save_state(root, state)
    _log(f"P6 完成：{len(stats)} 物种，新肽 {tot}"
         f"（L1 ab initio {l1} / L2 六框 {l2} / L3 N 端延伸 {l3} / L6 转录证据 {tx}；"
         f"Tier {t1}/{t2}/{t3}）。产物 05_novel/")
    return 0


def list_species_novel(root):
    """从 05_novel/05b_research 子目录推断物种（不写死物种代码）。"""
    out = []
    for d in ("05_novel", "05b_research"):
        p = os.path.join(root, d)
        if os.path.isdir(p):
            out += [x for x in sorted(os.listdir(p))
                    if os.path.isdir(os.path.join(p, x))]
    seen, uniq = set(), []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq


def cmd_novgene(args):
    """P8 事件分类：可信新肽 → 新编码基因 / 已注释基因校正（含 NG 图表）。"""
    from .core import novgene as G
    if args.from_dir:
        os.environ["ORFSCAPE_VALIDATE_DIR"] = args.from_dir
    _log(f"ORFscape v{__version__} . P8 新基因事件分类")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, args.out or "07_novgene")
    counts, subs, _ = G.build_novgene(root, out_dir, log=_log)
    _log(f"P8 完成：产物 {out_dir}")
    return 0


def cmd_novannot(args):
    """P8 注释：结构/理化/同源（Swiss-Prot + ESM-2 接触图）。"""
    from .core import novannot as A
    if args.from_dir:
        os.environ["ORFSCAPE_VALIDATE_DIR"] = args.from_dir
    _log(f"ORFscape v{__version__} . P8 新 ORF 结构/功能/同源注释")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, args.out or "08_annotation")
    A.build_novannot(root, out_dir, use_esm=not args.no_esm, log=_log)
    _log(f"P8 注释完成：产物 {out_dir}")
    return 0


def cmd_novreport(args):
    """P8 汇总：NG/AN/筛选三表 → Word 报告。"""
    from .core import novreport as R
    _log(f"ORFscape v{__version__} . P8 新基因汇总报告")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, args.out or "07_novgene")
    R.build_novelreport(root, out_dir, log=_log)
    _log(f"P8 报告完成：产物 {out_dir}")
    return 0


def cmd_research(args):
    """P6b 候选 ORF 小库回搜：补新肽证据 → 合并证据重算 Tier → 05b_research/。"""
    from .core import research as R
    _log(f"ORFscape v{__version__} . P6b 候选 ORF 小库回搜开始")
    root = ensure_project(args.workdir)
    species = tuple(s for s in str(args.species).split(",") if s.strip())
    if not species:
        species = tuple(inputs.project_species(root))
        if not species:
            _log("  项目内没有物种（先跑 P6），P6b 无法确定回搜范围")
            return 2
    stats = R.run(root, species=species, thread=args.thread, workers=args.workers,
                  collect_only=args.collect_only, log=_log)
    if not stats:
        _log("P6b 未产出（检查 05_novel/*/novel_peptides.tsv）")
        return 2
    state = load_state(root)
    state.setdefault("stages", {})["research"] = {
        "done": True, "ts": _now(), "n_species": len(stats),
        "n_peptide": sum(s["n_total"] for s in stats),
        "n_new": sum(s["n_new_peptide"] for s in stats)}
    save_state(root, state)
    _log(f"P6b 完成：{len(stats)} 物种，合并后肽 "
         f"{sum(s['n_total'] for s in stats)}，其中新增 "
         f"{sum(s['n_new_peptide'] for s in stats)}。产物 05b_research/")
    return 0


def cmd_l6map(args):
    """L6 坐标回贴：Prodigal 区段内坐标 → 基因组坐标，并按 CDS 重叠过滤。"""
    from .core import dbbuild as B
    _log(f"ORFscape v{__version__} . L6 转录 ORF 坐标回贴")
    root = ensure_project(args.workdir)
    sp = args.species
    out_dir = os.path.join(root, "02_db", sp)
    os.makedirs(out_dir, exist_ok=True)
    cds = args.cds or os.path.join(root, "01_genome", sp, "genomic.gff")
    recs = B.l6_import(args.gff, args.faa, cds, out_dir,
                       min_aa=args.min_aa,
                       max_cds_overlap=args.max_cds_overlap,
                       log=_log)
    _log(f"L6 入层 {len(recs)} 条 → {os.path.join(out_dir, 'L6_tx_orfs.faa')}")
    return 0


def cmd_uniprot(args):
    """按基因组识别物种 → 下载 UniProt 参考蛋白质组（一基因一蛋白）。"""
    from .core import uniprot as U
    _log(f"ORFscape v{__version__} . UniProt 物种参考蛋白质组")
    root = ensure_project(args.workdir) if args.workdir else None
    sp = args.sp
    faa = args.faa or (os.path.join(root, "01_genome", sp, "protein.faa") if root else None)
    fna = args.fna or (os.path.join(root, "01_genome", sp, "genomic.fna") if root else None)
    gff = args.gff or (os.path.join(root, "01_genome", sp, "genomic.gff") if root else None)
    out = args.out or (os.path.join(root, "02_db", sp, "uniprot_ref.faa") if root else None)

    if args.dry_run:
        name, src = ((args.organism, "命令行指定") if args.organism
                     else U.detect_species(faa, fna, gff))
        if not name:
            _log("无法识别物种：请用 --organism 指定")
            return 2
        _log(f"物种：{name}（来源：{src}）")
        if args.proteome:
            _log(f"指定蛋白质组：{args.proteome}")
            return 0
        hit, _used, _all = U.resolve(name, log=_log)
        if not hit:
            _log("UniProt 未找到该物种的蛋白质组")
            return 2
        tax = hit.get("taxonomy", {}) or {}
        _log("→ %s | %s | %s | taxid %s | 蛋白 %s"
             % (hit.get("id"), hit.get("proteomeType"), tax.get("scientificName"),
                tax.get("taxonId"), hit.get("proteinCount")))
        return 0

    stats = U.build(faa=faa, fna=fna, gff=gff, out=out, species=args.organism,
                    proteome=args.proteome, dedupe=not args.no_dedupe, log=_log)
    _log("完成：%s（%d 条蛋白）" % (stats["path"], stats["proteins"]))
    return 0


def cmd_eukbuild(args):
    """真核建库：GFF3/GTF → 剪接蛋白层 + 连接肽层（+ 转录本 ORF 层）。"""
    from .core import eukbuild, eukpresets
    preset = eukpresets.get(args.preset) if args.preset else None
    if args.preset and not preset:
        _log(f"未知预设：{args.preset}（可选 {', '.join(eukpresets.PRESETS)}）")
        return 2
    if preset:
        args.sp = args.sp or preset["code"]
        paths = eukpresets.paths(os.path.abspath(args.workdir), preset["code"])
        args.fna = args.fna or paths["fna"]
        args.gff = args.gff or paths["gff"]
        args.ref_protein = args.ref_protein or paths["ref"]
        _log(f"预设 {preset['code']}：{preset['name']} / {preset['assembly']}"
             f"（来源 {preset['source']}）")
    missing = [n for n, v in (("--fna", args.fna), ("--gff", args.gff)) if not v]
    if missing:
        _log("缺参数：%s（或用 --preset）" % "、".join(missing))
        return 2
    for label, path in (("fna", args.fna), ("gff", args.gff)):
        if not os.path.isfile(path):
            _log(f"{label} 不存在：{path}")
            return 2
    if args.ref_protein and not os.path.isfile(args.ref_protein):
        args.ref_protein = None
    _log(f"ORFscape v{__version__} . 真核建库开始（{args.sp}）")
    root = ensure_project(args.workdir)
    man = eukbuild.build(root, args.sp, args.fna, args.gff,
                         transcripts=args.transcripts, flank_aa=args.flank,
                         min_aa=args.min_aa, ref_protein=args.ref_protein,
                         preset=(preset or {}).get("code"), log=_log)
    c = man["counts"]
    _log(f"真核建库完成：蛋白 {c['proteins']} / 连接肽 {c['junction_peptides']}"
         f"（唯一 {c['junction_unique']}）/ 转录本 ORF {c['tx_orfs']}")
    if "reference_check" in man:
        rc = man["reference_check"]
        _log(f"参考校验：{rc['exact']}/{rc['ours']} 完全一致（{rc['exact_rate']:.1%}），"
             f"含 1 aa 前缀差异 {rc['match_rate']:.1%}")
    return 0


def cmd_braker(args):
    """BRAKER3：生成运行脚本（Linux 执行）+ 把 GTF 解析成蛋白 FASTA 供下游用。"""
    from .core import braker as B
    _log(f"ORFscape v{__version__} . BRAKER3 真核基因模型校正")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, "09_braker")
    os.makedirs(out_dir, exist_ok=True)
    faa = os.path.join(out_dir, "L1_braker_filt.faa")

    if args.gff:
        prot = B.filter_models(B.gtf_to_protein(args.gff, args.genome), [])
        B.write_faa(faa, prot)
        _log(f"GTF 解析：{len(prot)} 条蛋白 → {faa}")
        return 0

    ok, info = B.check_env(args.braker_cmd)
    if not ok:
        _log(f"BRAKER3 不可用：{info}")
    script = os.path.join(out_dir, "run_braker3.sh")
    B.build_run_script(script, args.genome, " ".join(args.bam) if args.bam else "",
                       args.protein, out_dir, threads=args.threads,
                       species=args.species, braker_cmd=args.braker_cmd)
    _log(f"运行脚本已生成：{script}")
    if args.run:
        st, rc = B.run(script, log=_log)
        _log(f"BRAKER3 运行：{st}（rc={rc}）")
        if st != "OK":
            return 3
    _log("提示：在 Linux/WSL 里 `bash run_braker3.sh` 跑完，再用 --gff 解析 braker.gtf")
    return 0


def cmd_qcfigs(args):
    """QC 系列单图：raw 参数质控 5 维度各一张（文件完整性/设计/扫描/方法/IT+AGC）。"""
    from .core import qcfigs
    from .core.i18n import get_lang
    _log(f"ORFscape v{__version__} . QC 系列单图生成")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, "08_report_zh" if get_lang() == "zh" else "08_report")
    qcfigs.build_qcfigs(root, out_dir)
    _log(f"QC1-QC5 完成 → {out_dir}")
    return 0


def cmd_e1figs(args):
    """E1 QC 系列单图：已注释搜库结果蛋白质组质控 7 维度独立出图。"""
    from .core import e1figs
    from .core.i18n import get_lang
    _log(f"ORFscape v{__version__} . E1 QC 系列单图生成")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, "08_report_zh" if get_lang() == "zh" else "08_report")
    e1figs.build_e1figs(root, out_dir)
    _log(f"E1Q1-E1Q7 完成 → {out_dir}")
    return 0


def cmd_novelfigs(args):
    """新肽段筛选漏斗（分样本列）+ 参数比较 → 08_report/。"""
    from .core import novelfigs
    from .core.i18n import get_lang
    if getattr(args, "from_dir", None):
        os.environ["ORFSCAPE_VALIDATE_DIR"] = args.from_dir
    _log(f"ORFscape v{__version__} . 新肽段筛选图与参数比较")
    root = ensure_project(args.workdir)
    out_dir = os.path.join(root, "08_report_zh" if get_lang() == "zh" else "08_report")
    novelfigs.build_novelfigs(root, out_dir)
    from .core import novgene, novannot, novreport
    novgene.build_novgene(root, out_dir, log=_log)
    novannot.build_novannot(root, out_dir, use_esm=not getattr(args, "no_esm", False),
                            log=_log)
    novreport.build_novelreport(root, out_dir, log=_log)
    _log(f"NOV 图/表/结论 完成 → {out_dir}")
    return 0


def cmd_run(args):
    """一键跑：委托给 orchestrator（事件流 + 阶段自检）。"""
    from . import orchestrator
    from .core.i18n import set_lang
    if getattr(args, "lang", None):
        set_lang(args.lang)
    if not args.thread:
        args.thread = max(4, (os.cpu_count() or 8) // 2)
    if not args.workers:
        args.workers = 4
    return orchestrator.run(args)


def cmd_gui(args):
    """图形界面。"""
    try:
        from .gui.app import main as gui_main
    except ImportError as exc:
        _log(f"GUI 不可用（缺 tkinter）: {exc}")
        return 2
    return gui_main(prefill_workdir=getattr(args, "workdir", None))


# ---------------------------------------------------------------- 阶段上下文
_STAGE_OF = {
    "scan": "P0", "genome": "P1", "split": "D0", "dbbuild": "C",
    "e1": "E1", "qc": "QC", "e2": "E2", "novel": "P6", "research": "P6b",
    "validate": "P7", "novgene": "P8", "novannot": "P8", "novreport": "P9",
}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="orfscape", description="ORFscape 原核蛋白质基因组学流程")
    ap.add_argument("--version", action="version", version=__version__)
    ap.add_argument("--lang", default=None, choices=["en", "zh"],
                    help="output language (default: en; or set ORFSCAPE_LANG=zh)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sc = sub.add_parser("scan", help="数据确定：raw 发现/分组/物种/方法参数/扫描统计")
    sc.add_argument("--raw", required=True, help="raw 文件根目录（递归发现）")
    sc.add_argument("--workdir", default=None, help="项目目录（默认 orfscape_projects/proj_时间戳）")
    sc.add_argument("--no-sha256", action="store_true", help="跳过 sha256（提速）")
    sc.add_argument("--skip-dump", action="store_true", help="跳过 rawdump 扫描表")
    sc.add_argument("--skip-method", action="store_true", help="跳过 rawmethod 方法提取")
    sc.add_argument("--offline", action="store_true", help="离线模式（跳过 taxid 在线解析）")
    sc.add_argument("--species-name", default=None,
                    help="物种拉丁名（覆盖文件名推断，如 Treponema pallidum）")
    sc.add_argument("--species-code", default=None, help="物种代码（覆盖文件名推断）")
    sc.set_defaults(func=cmd_scan)
    dc = sub.add_parser("doctor", help="环境检测器：工具/引擎自动发现与版本探测")
    dc.add_argument("--workdir", default="orfscape_projects/envcheck", help="输出目录")
    dc.set_defaults(func=cmd_doctor)
    gc = sub.add_parser("genome", help="基因组获取：按项目物种表下载 RefSeq 三件套并钉扎")
    gc.add_argument("--workdir", required=True, help="项目目录")
    gc.set_defaults(func=cmd_genome)
    sp = sub.add_parser("split", help="D0 电压拆分：多 CV raw → 分 CV MGF（含对账）")
    sp.add_argument("--workdir", required=True, help="项目目录")
    sp.set_defaults(func=cmd_split)
    dbp = sub.add_parser("dbbuild", help="C 建库：L0-L5 层级 + 两段漏斗拆分库")
    dbp.add_argument("--workdir", required=True, help="项目目录")
    dbp.set_defaults(func=cmd_dbbuild)
    e1p = sub.add_parser("e1", help="E1 已注释库搜库（pFind x db_annot x 分CV MGF）")
    e1p.add_argument("--workdir", required=True, help="项目目录")
    e1p.add_argument("--thread", type=int, default=16, help="pFind 线程数")
    e1p.set_defaults(func=cmd_e1)
    e2p = sub.add_parser("e2", help="E2 发现搜库（db_aug/db_6f，层归表头前缀）")
    e2p.add_argument("--workdir", required=True, help="项目目录")
    e2p.add_argument("--dbs", default="aug,6f", help="库后缀逗号表（aug,6f）")
    e2p.add_argument("--thread", type=int, default=12, help="pFind 线程数")
    e2p.set_defaults(func=cmd_e2)
    qcp = sub.add_parser("qc", help="QC（紧跟 E1）：酶切质控+门槛警告")
    qcp.add_argument("--workdir", required=True, help="项目目录")
    qcp.set_defaults(func=cmd_qc)
    qf = sub.add_parser("qcfigs", help="QC 系列单图：raw 参数质控 5 维度独立出图")
    qf.add_argument("--workdir", required=True, help="项目目录")
    qf.set_defaults(func=cmd_qcfigs)
    ef = sub.add_parser("e1figs", help="E1 QC 系列单图：搜库结果蛋白质组质控 7 维度")
    ef.add_argument("--workdir", required=True, help="项目目录")
    ef.set_defaults(func=cmd_e1figs)
    nf = sub.add_parser("novelfigs", help="新肽段筛选漏斗（分样本）+ 参数比较")
    nf.add_argument("--workdir", required=True, help="项目目录")
    nf.add_argument("--no-esm", action="store_true",
                    help="跳过 ESM-2 接触图（无权重或求快）")
    nf.add_argument("--from", dest="from_dir", default=None,
                    help="验收集目录（默认 06_validate；TP 用 06b_validate_research）")
    nf.set_defaults(func=cmd_novelfigs)
    vp = sub.add_parser("validate", help="R8 收口：novel 全量过 isoPSV+Tier 分流门控")
    vp.add_argument("--workdir", required=True, help="项目目录")
    vp.add_argument("--isopsv", default=V_ISOPSV_DEFAULT,
                    help="isoPSV_multiengine 安装根目录")
    vp.add_argument("--min-score", default="auto",
                    help="Tier3 谱图质量卡值（数值或 auto=Tier1 中位标定）")
    vp.add_argument("--min-cos", default="auto",
                    help="Tier3 理论谱图相似性卡值（数值或 auto=Tier1 p25 标定）")
    vp.add_argument("--smoke", action="store_true", help="冒烟：只验第一个物种 × TL 子集")
    vp.add_argument("--from-research", action="store_true",
                    help="用 P6b 合并集（05b_research，tier_new）收口 → 06b_validate_research")
    vp.add_argument("--reuse", default=None,
                    help="增量模式：复用该目录的已验结论（如 06_validate），只验新增肽")
    vp.add_argument("--cards", default=None,
                    help="增量模式卡值 'score,cos'（如 0.28,0.82，沿用旧标定）")
    vp.add_argument("--out", default=None,
                    help="输出目录（06_validate / 06b_validate_research / 06c_validate_full）")
    vp.set_defaults(func=cmd_validate)
    rp = sub.add_parser("research", help="P6b 候选 ORF 小库回搜（补肽证据、重算 Tier）")
    rp.add_argument("--workdir", required=True, help="项目目录")
    rp.add_argument("--species", default="", help="物种逗号表（默认项目内全部物种）")
    rp.add_argument("--thread", type=int, default=5, help="单任务 pFind 线程")
    rp.add_argument("--workers", type=int, default=4, help="并发任务数")
    rp.add_argument("--collect-only", action="store_true",
                    help="只合并已有结果并重算 Tier，不跑搜库")
    rp.set_defaults(func=cmd_research)
    nv = sub.add_parser("novel",
                        help="P6 新肽筛选（E2 层命中 → 四滤 → 坐标 → Tier）")
    nv.add_argument("--workdir", required=True, help="项目目录")
    nv.add_argument("--species", default="",
                    help="物种逗号表（留空则从 05_novel/05b_research 子目录推断）")
    nv.set_defaults(func=cmd_novel)
    for _name, _func, _help in (
            ("novgene", cmd_novgene, "P8 事件分类（新编码基因 / 已注释基因校正）"),
            ("novannot", cmd_novannot, "P8 新 ORF 结构/功能/同源注释"),
            ("novreport", cmd_novreport, "P8 新基因汇总 Word 报告")):
        _p = sub.add_parser(_name, help=_help)
        _p.add_argument("--workdir", required=True, help="项目目录")
        _p.add_argument("--out", default=None, help="输出子目录（默认见各命令）")
        _p.add_argument("--from", dest="from_dir", default=None,
                        help="验收集目录（默认 06_validate；TP 用 06b_validate_research）")
        if _name == "novannot":
            _p.add_argument("--no-esm", action="store_true",
                            help="跳过 ESM-2 接触图（无 GPU/离线时）")
        _p.set_defaults(func=_func)
    bp = sub.add_parser("braker",
                        help="BRAKER3 真核基因模型校正（B/C 线；生成脚本 + 解析 GTF）")
    bp.add_argument("--workdir", required=True, help="项目目录")
    bp.add_argument("--genome", required=True, help="基因组 FASTA")
    bp.add_argument("--bam", nargs="*", default=[], help="RNA-seq BAM（可多个）")
    bp.add_argument("--protein", default=None, help="跨物种蛋白 FASTA（可选，走 EP+ 分支）")
    bp.add_argument("--species", default="arabidopsis", help="BRAKER species 名")
    bp.add_argument("--threads", type=int, default=8)
    bp.add_argument("--braker-cmd", default="braker3", help="BRAKER3 可执行命令")
    bp.add_argument("--run", action="store_true", help="就地执行（需 Linux）")
    bp.add_argument("--gff", default=None, help="已有 GTF：只解析成蛋白 FASTA")
    bp.set_defaults(func=cmd_braker)
    eb = sub.add_parser("eukbuild",
                        help="真核建库：GFF3/GTF → 剪接蛋白层 + 剪接连接肽层（+ 转录本 ORF）")
    eb.add_argument("--workdir", required=True, help="项目目录")
    eb.add_argument("--preset", default=None,
                    help="示例物种预设：CREI（衣藻）/ HSAP（人），自动定位 01_genome/<code>/")
    eb.add_argument("--sp", default=None, help="物种代码（写 02_db/<sp>）")
    eb.add_argument("--fna", default=None, help="基因组 FASTA")
    eb.add_argument("--gff", default=None, help="注释 GFF3 或 GTF（可 .gz）")
    eb.add_argument("--transcripts", default=None,
                    help="转录本 FASTA（可选，给了就建 L2 三框 ORF 层）")
    eb.add_argument("--flank", type=int, default=12,
                    help="连接肽两侧各取多少个 aa（默认 12）")
    eb.add_argument("--min-aa", type=int, default=30, help="最短蛋白/ORF（默认 30）")
    eb.add_argument("--ref-protein", default=None,
                    help="参考蛋白 FASTA（可选，做剪接正确性校验）")
    eb.set_defaults(func=cmd_eukbuild)
    lp = sub.add_parser("l6map",
                        help="L6 转录 ORF：Prodigal 区段内坐标 → 基因组坐标 + CDS 重叠过滤")
    lp.add_argument("--workdir", required=True, help="项目目录")
    lp.add_argument("--species", default="TP", help="物种代码（定位 02_db/<sp>）")
    lp.add_argument("--gff", required=True, help="Prodigal GFF（区段内坐标）")
    lp.add_argument("--faa", required=True, help="Prodigal FAA")
    lp.add_argument("--cds", default=None,
                    help="注释 GFF（默认 01_genome/<sp>/genomic.gff）")
    lp.add_argument("--min-aa", type=int, default=30)
    lp.add_argument("--max-cds-overlap", type=float, default=0.8,
                    help="同链与注释 CDS 重叠比例超过此值判为已注释")
    lp.set_defaults(func=cmd_l6map)
    up = sub.add_parser("uniprot",
                        help="按基因组识别物种并下载 UniProt 参考蛋白质组（一基因一蛋白）")
    up.add_argument("--workdir", default=None, help="项目目录")
    up.add_argument("--sp", default="TP",
                    help="物种代码（定位 01_genome/<sp> 与 02_db/<sp>）")
    up.add_argument("--faa", default=None, help="蛋白 FASTA（默认 01_genome/<sp>/protein.faa）")
    up.add_argument("--fna", default=None, help="基因组 FASTA（默认 01_genome/<sp>/genomic.fna）")
    up.add_argument("--gff", default=None, help="注释 GFF（默认 01_genome/<sp>/genomic.gff）")
    up.add_argument("--organism", default=None, help="手工指定物种名（跳过自动识别）")
    up.add_argument("--proteome", default=None,
                    help="手工指定 UniProt 蛋白质组编号（如 UP000000811，跳过检索）")
    up.add_argument("--out", default=None,
                    help="输出 FASTA（默认 02_db/<sp>/uniprot_ref.faa）")
    up.add_argument("--dry-run", action="store_true",
                    help="只识别物种 + 查蛋白质组，不下载")
    up.add_argument("--no-dedupe", action="store_true", help="不做一基因一蛋白合并")
    up.set_defaults(func=cmd_uniprot)
    rp = sub.add_parser("run", help="一键跑 P0→P9（流式日志 + 事件流 + 阶段自检）")
    rp.add_argument("--workdir", required=True, help="项目目录")
    rp.add_argument("--raw", default=None, help="raw 根目录/单文件（P0 用）")
    rp.add_argument("--species-name", default=None,
                    help="物种拉丁名（覆盖文件名推断，如 Treponema pallidum）")
    rp.add_argument("--species-code", default=None, help="物种代码（默认由文件名推断）")
    rp.add_argument("--desc", default=None, help="自由描述（写进事件流）")
    rp.add_argument("--from", dest="start", default="P0", help="起始步（默认 P0）")
    rp.add_argument("--to", dest="end", default="P8c", help="结束步（默认 P8c）")
    rp.add_argument("--skip", default="", help="跳过步，逗号分隔（如 P6b,QC）")
    rp.add_argument("--thread", type=int, default=0, help="pFind 线程数")
    rp.add_argument("--workers", type=int, default=0, help="P6b 并发任务数")
    rp.add_argument("--dbs", default="aug,6f", help="E2 库后缀")
    rp.add_argument("--min-cos", default="auto", help="Tier3 相似性卡值")
    rp.add_argument("--reuse", default=None, help="P7 增量：复用旧验证目录")
    rp.add_argument("--cards", default=None, help="P7 增量卡值 'score,cos'")
    rp.add_argument("--no-esm", action="store_true", help="novannot 跳过 ESM-2")
    rp.add_argument("--force", action="store_true", help="忽略完成标记")
    rp.add_argument("--dry-run", action="store_true", help="只打印命令")
    rp.add_argument("--check-script", default=None,
                    help="外部阶段自检脚本（如 skill 的 stage_check.py）")
    rp.add_argument("--lang", default=None, choices=["en", "zh"],
                    help="输出语言（默认跟随全局 --lang / ORFSCAPE_LANG）")
    rp.set_defaults(func=cmd_run)
    gp = sub.add_parser("gui", help="图形界面（英文默认，可切中文）")
    gp.add_argument("--workdir", default=None, help="预填项目目录")
    gp.set_defaults(func=cmd_gui)
    args = ap.parse_args(argv)
    # 语言生效（--lang > ORFSCAPE_LANG > 默认 en），并写入项目 run_state 持久化
    from .core.i18n import set_lang, get_lang
    if getattr(args, "lang", None):
        set_lang(args.lang)
    if args.cmd == "scan" and not args.workdir:
        args.workdir = os.path.join("orfscape_projects", f"proj_{datetime.datetime.now():%Y%m%d%H%M%S}")
    # 语言登记到 run_state（workdir 存在时）
    try:
        _root = ensure_project(args.workdir) if args.workdir else None
        if _root:
            _st = load_state(_root)
            _st["language"] = get_lang()
            save_state(_root, _st)
    except Exception:
        pass
    # 阶段上下文：让 _log 的事件流带上阶段号（GUI 用来归行）
    from . import progress
    progress.set_stage(_STAGE_OF.get(args.cmd))
    try:
        return args.func(args)
    finally:
        progress.close()


if __name__ == "__main__":
    sys.exit(main())
