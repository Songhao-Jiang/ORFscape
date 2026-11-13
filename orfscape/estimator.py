# -*- coding: utf-8 -*-
"""资源/体积/耗时估算与数据库面板数据。

三个来源，优先级从高到低：
1. 项目里已经落地的实测值（``00_input/raw_info.tsv``、``02_db`` 文件大小、
   ``01_genome/*/genomic.fna``）——读到就用真值；
2. 本项目上一次跑完的实测阶段耗时（``00_logs/calibration.json``）；
3. 参考项目（TP，36 raw / 6.03 GB / 495,483 MS2）标定的先验系数。

先验只用于"开跑前给个量级"，跑起来以后 ETA 由 :mod:`orfscape.runlog` 用实测
速率接管。所有估算都标注来源，界面区分 ``measured`` 与 ``estimated``。
"""
import json
import os

from .core import netguard

GB = 1024.0 ** 3
MB = 1024.0 ** 2

# ---------------------------------------------------------------- 先验系数
# 参考项目：TP 36 raw / 6.03 GB / 495,483 MS2 / 基因组 1.14 Mb；
# E1 36 任务 26.6 min（44 s/任务，1 MB 库）；E2 144 任务 5.43 h（136 s/任务，6 MB 库）。
PRIOR = {
    "raw_avg_mb": 168.0,
    "ms2_per_gb": 82170.0,        # 实测 495,483 / 6.03 GB
    "mgf_factor": 6.97,           # MGF 总量 / raw 总量
    "db_mb_per_genome_mb": 83.0,  # 全部层+索引 / 基因组 Mb
    "e1_sec_task": 44.0,
    "e1_db_mb": 1.0,
    "e1_out_mb_task": 28.0,
    "e2_sec_task": 136.0,
    "e2_db_mb": 6.0,
    "e2_out_mb_task": 28.5,
    "val_sec_peptide": 2.8,
    "val_out_factor": 0.35,
    "default_genome_mb": 4.0,     # 未知物种的中位细菌基因组
    "db_exp": 0.65,               # 搜库时间对库大小次线性
}

STAGES = ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b", "P7",
          "P8a", "P8b", "P8c"]

# 权重仅用于"没有耗时先验时"的粗粒度进度条
WEIGHTS = {"P0": 1.0, "P1": 3.0, "D0": 10.0, "C": 2.0, "E1": 9.0, "QC": 0.4,
           "E2": 30.0, "P6": 1.5, "P6b": 14.0, "P7": 12.0, "P8a": 0.6,
           "P8b": 4.0, "P8c": 2.0}


# ---------------------------------------------------------------- 格式化
def fmt_bytes(n, decimals=1):
    if n is None:
        return "n/a"
    n = float(n)
    for unit, size in (("TB", 1024.0 ** 4), ("GB", GB), ("MB", MB), ("KB", 1024.0)):
        if abs(n) >= size:
            return "%.*f %s" % (decimals, n / size, unit)
    return "%d B" % n


def fmt_sec(s):
    """人类可读时长：``2 h 14 min`` / ``48 s``。"""
    if s is None:
        return "n/a"
    s = int(max(0, round(s)))
    if s < 10:
        return "%d s" % s
    if s < 90:
        return "%d s" % s
    m, sec = divmod(s, 60)
    if m < 60:
        return "%d min" % m if sec < 30 else "%d min %02d s" % (m, sec)
    h, m = divmod(m, 60)
    return "%d h %02d min" % (h, m)


# ---------------------------------------------------------------- 输入
def scan_raws(paths):
    """把用户给的路径（文件/目录）展开成 raw 列表。只认 .raw，去符号链接。"""
    found = []
    seen = set()
    for p in paths or []:
        p = os.path.abspath(str(p).strip().strip('"'))
        if os.path.isdir(p):
            for root, dirs, files in os.walk(p):
                dirs[:] = [d for d in dirs if not d.startswith(".")]
                for f in files:
                    if f.lower().endswith(".raw"):
                        fp = os.path.realpath(os.path.join(root, f))
                        if fp not in seen:
                            seen.add(fp)
                            found.append(fp)
        elif os.path.isfile(p):
            fp = os.path.realpath(p)
            if fp not in seen:
                seen.add(fp)
                found.append(fp)
    found.sort()
    total = sum(os.path.getsize(f) for f in found if os.path.isfile(f))
    return {"files": found, "n": len(found), "bytes": total}


def read_raw_info(workdir):
    """读已扫描的 raw_info.tsv（有就是真值）。"""
    p = os.path.join(workdir, "00_input", "raw_info.tsv")
    if not os.path.isfile(p):
        return None
    n = 0
    ms2 = 0
    size_mb = 0.0
    try:
        with open(p, encoding="utf-8-sig") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                i_size, i_ms2 = head.index("size_mb"), head.index("n_ms2")
            except ValueError:
                return None
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) <= max(i_size, i_ms2):
                    continue
                n += 1
                try:
                    size_mb += float(c[i_size] or 0)
                    ms2 += int(float(c[i_ms2] or 0))
                except ValueError:
                    pass
    except OSError:
        return None
    return {"n": n, "ms2": ms2, "bytes": size_mb * MB}


def genome_mb(workdir, sp=None):
    """基因组大小（Mb）：优先读 genomic.fna。"""
    base = os.path.join(workdir, "01_genome")
    if not os.path.isdir(base):
        return None
    cands = []
    for name in sorted(os.listdir(base)):
        d = os.path.join(base, name)
        if sp and name != sp:
            continue
        f = os.path.join(d, "genomic.fna")
        if os.path.isfile(f):
            cands.append((name, f))
    if not cands:
        return None
    name, f = cands[0]
    try:
        with open(f, "rb") as fh:
            n = 0
            for line in fh:
                if not line.startswith(b">"):
                    n += len(line.strip())
        return {"code": name, "mb": round(n / 1e6, 2)}
    except OSError:
        return None


def db_files(workdir, sp=None):
    """列出 02_db/<sp> 下的层与库文件（名, 字节）。"""
    base = os.path.join(workdir, "02_db")
    if not os.path.isdir(base):
        return []
    dirs = [d for d in sorted(os.listdir(base))
            if os.path.isdir(os.path.join(base, d)) and (not sp or d == sp)]
    if not dirs:
        return []
    out = []
    d = os.path.join(base, dirs[0])
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if os.path.isfile(p) and f.endswith(".faa"):
            out.append({"name": f, "bytes": os.path.getsize(p), "dir": d})
    return out


# ---------------------------------------------------------------- 标定
def calibration_path(workdir):
    return os.path.join(workdir, "00_logs", "calibration.json")


def load_calibration(workdir):
    p = calibration_path(workdir)
    if not os.path.isfile(p):
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def save_calibration(workdir, state):
    """把一次跑完的实测值写进 calibration.json（供下次预估）。"""
    stages = {}
    for sid, st in (state.stages or {}).items():
        # 只记真正跑完的阶段：跳过的（sec=0）和失败的会污染下次预估
        if getattr(st, "status", "") != "ok":
            continue
        if st.t0 is None or st.t1 is None:
            continue
        sec = st.t1 - st.t0
        if sec < 2:
            continue
        stages[sid] = {"sec": round(sec, 1), "done": st.done, "total": st.total}
    if not stages:
        return None
    p = calibration_path(workdir)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    data = load_calibration(workdir)
    prev = data.get("stages") or {}
    prev.update(stages)
    data["stages"] = prev
    data["updated"] = __import__("time").strftime("%Y-%m-%dT%H:%M:%S")
    try:
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(data, fh, ensure_ascii=False, indent=1)
    except OSError:
        return None
    return p


# ---------------------------------------------------------------- 估算
def _tasks(workdir, dbs):
    """任务数：E1 = MGF 数；E2 = MGF 数 × 库数。"""
    mgf_dir = os.path.join(workdir, "03_mgf")
    n_mgf = 0
    if os.path.isdir(mgf_dir):
        n_mgf = sum(1 for f in os.listdir(mgf_dir) if f.lower().endswith(".mgf"))
    ndb = len([d for d in (dbs or "").split(",") if d.strip()]) or 2
    return n_mgf, n_mgf * ndb


def estimate(workdir, raw_paths=None, threads=16, workers=4, dbs="aug,6f",
             species_code=None):
    """返回估算结果。

    键：inputs / genome / db / tasks / disk{...} / ram / est_sec{stage} /
    weights / warnings / source
    """
    workdir = os.path.abspath(str(workdir)) if workdir else ""
    cal = load_calibration(workdir) if workdir else {}
    cal_stages = (cal.get("stages") or {})

    # ---- 输入规模
    info = read_raw_info(workdir) if workdir else None
    src = "measured" if info and info.get("n") else "estimated"
    if not info or not info.get("n"):
        sc = scan_raws(raw_paths or [])
        n_raw = sc["n"]
        raw_bytes = sc["bytes"]
        info = {"n": n_raw, "bytes": raw_bytes,
                "ms2": int(raw_bytes / GB * PRIOR["ms2_per_gb"]) if raw_bytes else 0}
    raw_bytes = float(info.get("bytes") or 0)
    n_raw = int(info.get("n") or 0)
    ms2 = int(info.get("ms2") or 0)
    if not ms2 and raw_bytes:
        ms2 = int(raw_bytes / GB * PRIOR["ms2_per_gb"])
    ms2_per_file = (ms2 / n_raw) if n_raw else 13800.0

    # ---- 基因组与库
    g = genome_mb(workdir, species_code) if workdir else None
    genome_mb_v = float(g["mb"]) if g else PRIOR["default_genome_mb"]
    files = db_files(workdir, species_code) if workdir else []

    def fsize(name):
        for f in files:
            if f["name"] == name:
                return float(f["bytes"]) / MB
        return None

    db_annot_mb = fsize("db_annot.faa")
    db_e2_mb = None
    if fsize("db_aug.faa") is not None or fsize("db_6f.faa") is not None:
        db_e2_mb = (fsize("db_aug.faa") or 0) + (fsize("db_6f.faa") or 0)
    db_total_mb = sum(f["bytes"] for f in files) / MB if files else None
    db_source = "measured" if db_total_mb else "estimated"
    if db_total_mb is None:
        db_total_mb = PRIOR["db_mb_per_genome_mb"] * genome_mb_v
    if db_annot_mb is None:
        db_annot_mb = max(0.3, 0.9 * genome_mb_v)
    if db_e2_mb is None:
        db_e2_mb = max(1.0, 3.2 * genome_mb_v)

    # ---- 任务数
    n_mgf, n_e2 = _tasks(workdir, dbs) if workdir else (0, 0)
    if not n_mgf:
        n_mgf = n_raw * 1            # 单 CV 情形；多 CV 由 scan 后回填
    if not n_e2:
        n_e2 = n_mgf * max(1, len([d for d in (dbs or "").split(",") if d.strip()]))

    # ---- 时间
    def scaled(base_key, db_mb, db_ref_key, per_task=True):
        base = PRIOR[base_key]
        ref = PRIOR[db_ref_key]
        f_ms2 = max(0.05, ms2_per_file / 13800.0)
        f_db = max(0.05, (db_mb / ref) ** PRIOR["db_exp"])
        return base * f_ms2 * f_db

    e1_task = scaled("e1_sec_task", db_annot_mb, "e1_db_mb")
    e2_task = scaled("e2_sec_task", db_e2_mb, "e2_db_mb")

    est = {}
    def put(sid, sec, why="prior"):
        if sid in cal_stages and cal_stages[sid].get("sec"):
            est[sid] = {"sec": float(cal_stages[sid]["sec"]), "src": "calibrated"}
        else:
            est[sid] = {"sec": float(max(1.0, sec)), "src": why}

    mgf_mb = raw_bytes / MB * PRIOR["mgf_factor"]
    put("P0", 5.0 * n_raw + 5)
    put("P1", 30.0 + 2.0 * genome_mb_v * 0.5)
    put("D0", 3.0 * n_raw + mgf_mb / 60.0)
    put("C", 60.0 + 15.0 * 6)
    put("E1", n_mgf * e1_task)
    put("QC", max(10.0, n_mgf * e1_task * 0.05))
    put("E2", n_e2 * e2_task)
    put("P6", 30.0 + ms2 / 40000.0)
    # P6b：候选 ORF 小库回搜，库小但任务多
    n_p6b = n_mgf * max(1, len([d for d in (dbs or "").split(",") if d.strip()])) // 2
    put("P6b", max(120.0, n_p6b * e2_task * 0.55 / max(1, workers)))
    put("P7", 300.0 + max(0, ms2 // 1000) * PRIOR["val_sec_peptide"] * 0.25)
    put("P8a", 60.0)
    put("P8b", 600.0)
    put("P8c", 120.0)

    # ---- 磁盘
    e1_bytes = n_mgf * PRIOR["e1_out_mb_task"] * MB
    e2_bytes = n_e2 * PRIOR["e2_out_mb_task"] * MB
    disk = {
        "raw": raw_bytes,
        "mgf": raw_bytes * PRIOR["mgf_factor"],
        "db": db_total_mb * MB,
        "e1": e1_bytes,
        "e2": e2_bytes,
        "validate": max(200 * MB, e1_bytes * PRIOR["val_out_factor"]),
        "logs": max(50 * MB, raw_bytes * 0.02),
    }
    disk["total"] = sum(disk.values())
    disk["needed_free"] = disk["total"] * 1.10
    disk["source"] = "measured" if src == "measured" else "estimated"

    # ---- 内存（pFind 峰值随线程线性增长的经验式）
    ram_need = 2.0 + 0.45 * float(threads)
    if workers and workers > 1:
        ram_need += 0.8 * float(workers)     # P6b 多任务叠加
    ram = {"need_gb": round(ram_need, 1), "threads": int(threads),
           "workers": int(workers)}

    # ---- 告警
    warnings = []
    try:
        import shutil
        free = shutil.disk_usage(workdir or os.path.abspath(".")).free / GB
        disk["free_gb"] = round(free, 1)
        if free < disk["needed_free"] / GB:
            warnings.append({"key": "warn.disk", "args": {
                "need": fmt_bytes(disk["needed_free"]), "free": fmt_bytes(free * GB)}})
    except OSError:
        disk["free_gb"] = None
    try:
        from .sysmon import snapshot
        snap = snapshot(workdir or os.getcwd())
        if snap.get("mem_total_gb") and ram_need > snap["mem_total_gb"] * 0.8:
            warnings.append({"key": "warn.ram", "args": {
                "need": "%.1f" % ram_need, "total": "%.1f" % snap["mem_total_gb"]}})
        cpu = snap.get("cpu_n") or 0
        if threads > cpu:
            warnings.append({"key": "warn.threads", "args": {
                "threads": threads, "cpu": cpu}})
        ram["sys_total_gb"] = snap.get("mem_total_gb")
        ram["sys_free_gb"] = round((snap.get("mem_total_gb") or 0)
                                   - (snap.get("mem_used_gb") or 0), 1)
    except Exception:            # 采样失败不影响估算
        pass
    if not n_raw:
        warnings.append({"key": "warn.noraw", "args": {}})

    total_sec = sum(v["sec"] for v in est.values())
    return {
        "inputs": {"n_raw": n_raw, "bytes": raw_bytes, "ms2": ms2,
                   "ms2_per_file": int(ms2_per_file), "source": src},
        "genome": {"mb": genome_mb_v, "code": (g or {}).get("code"),
                   "source": "measured" if g else "assumed"},
        "db": {"annot_mb": round(db_annot_mb, 1), "e2_mb": round(db_e2_mb, 1),
               "total_mb": round(db_total_mb, 1), "source": db_source},
        "tasks": {"mgf": n_mgf, "e1": n_mgf, "e2": n_e2, "p6b": n_p6b},
        "disk": disk, "ram": ram,
        "est_sec": est, "weights": dict(WEIGHTS),
        "total_sec": total_sec, "warnings": warnings,
        "calibrated": bool(cal_stages),
    }


# ---------------------------------------------------------------- 数据库面板
def database_rows(workdir, species_code=None, species_name=None,
                  raw_info=None):
    """返回数据库面板的行（本地已有 + 待下载），带大小与状态。"""
    workdir = os.path.abspath(str(workdir)) if workdir else ""
    rows = []
    base = os.path.join(workdir, "01_genome") if workdir else ""
    codes = [species_code] if species_code else []
    if not codes and os.path.isdir(base):
        codes = [d for d in sorted(os.listdir(base))
                 if os.path.isdir(os.path.join(base, d)) and not d.startswith(".")]
    for code in codes:
        for name, label in (("genomic.fna", "genome.fna"),
                            ("genomic.gff", "genome.gff"),
                            ("protein.faa", "protein.faa")):
            p = os.path.join(base, code, name)
            rows.append({"key": "%s/%s" % (code, name), "group": "refseq",
                         "label": label, "kind": "remote", "path": p,
                         "bytes": os.path.getsize(p) if os.path.isfile(p) else None,
                         "state": "present" if os.path.isfile(p) else "pending"})
        vp = os.path.join(base, code, "versions.json")
        if os.path.isfile(vp):
            try:
                meta = json.load(open(vp, encoding="utf-8"))
                for r in rows:
                    if r["key"].startswith(code + "/") and r["bytes"] is None:
                        nm = r["key"].split("/", 1)[1]
                        sz = ((meta.get("files") or {}).get(nm) or {}).get("size")
                        if sz:
                            r["bytes"] = sz          # gz 大小，已是实测
            except (OSError, ValueError):
                pass
    # UniProt 参考库
    up = os.path.join(workdir, "02_db", species_code or "", "uniprot_ref.faa") \
        if workdir and species_code else ""
    if up and os.path.isfile(up):
        rows.append({"key": "uniprot", "group": "uniprot",
                     "label": "uniprot_ref.faa", "kind": "remote", "path": up,
                     "bytes": os.path.getsize(up), "state": "present"})
    elif workdir:
        rows.append({"key": "uniprot", "group": "uniprot",
                     "label": "UniProt reference proteome", "kind": "remote",
                     "path": None, "bytes": None, "state": "unknown"})
    # 本地层
    for f in db_files(workdir, species_code):
        nm = f["name"]
        if nm.startswith(("L", "db_")):
            rows.append({"key": "db/" + nm, "group": "layers", "label": nm,
                         "kind": "local", "bytes": f["bytes"], "state": "present",
                         "path": os.path.join(f["dir"], nm)})
    return rows


def probe_remote_sizes(rows, species_name=None, code=None, log=print):
    """给 remote 行补大小（HEAD 探测）。

    - RefSeq：先解析装配 accession（:func:`orfscape.core.genome.pick_assembly`），
      再对三件套 .gz 做 HEAD；
    - UniProt：查 proteome，用蛋白条数估（stream 接口通常不给长度）。
    网络失败一律留 None，由界面显示 n/a。
    """
    out = [dict(r) for r in rows]
    by = {(r.get("group"), r.get("label")): r for r in out}
    if species_name:
        try:
            from .core import genome as gmod
            acc, cat, org, ftp = gmod.pick_assembly(species_name)
            base = ftp.rstrip("/")
            stem = base.rsplit("/", 1)[-1]
            for label, suffix in (("genome.fna", "_genomic.fna.gz"),
                                  ("genome.gff", "_genomic.gff.gz"),
                                  ("protein.faa", "_protein.faa.gz")):
                r = by.get(("refseq", label))
                if not r:
                    continue
                r["accession"] = acc
                if r.get("state") == "present" and r.get("bytes"):
                    continue
                url = "%s/%s%s" % (base, stem, suffix)
                try:
                    r["bytes"] = netguard.head_size(url, {"ftp.ncbi.nlm.nih.gov"})
                    r["state"] = "remote"
                    r["url"] = url
                except (netguard.UrlRejected, OSError):
                    r["state"] = "unknown"
                log("probe %s -> %s" % (label, fmt_bytes(r.get("bytes"))))
        except Exception as exc:                  # 解析/网络失败都不致命
            log("probe refseq failed: %s" % exc)
    if species_name:
        try:
            from .core import uniprot as umod
            hit, used, _ = umod.resolve(species_name, log=lambda *_: None)
            n_prot = int(hit.get("proteinCount") or 0)
            r = by.get(("uniprot", "UniProt reference proteome")) or by.get(
                ("uniprot", "uniprot_ref.faa"))
            if r is not None and not r.get("bytes"):
                r["bytes"] = int(n_prot * 0.40 * 1024)   # ~0.40 KB/条
                r["state"] = "estimate"
                r["proteome"] = hit.get("id")
                r["n_proteins"] = n_prot
        except Exception as exc:
            log("probe uniprot failed: %s" % exc)
    return out
