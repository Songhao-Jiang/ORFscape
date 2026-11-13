# -*- coding: utf-8 -*-
"""SRA 转录组数据筛选：物种 → 非基因编辑 RNA-seq run 清单（XML 解析版）。

网络限制：仅 https、域名白名单（eutils.ncbi.nlm.nih.gov）、DNS 公网校验、禁重定向。
"""
import ipaddress
import json
import re
import socket
import urllib.parse
import urllib.request

_ALLOWED = {"eutils.ncbi.nlm.nih.gov"}


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise RuntimeError(f"禁止重定向: {code}")


def _get(url, timeout=30):
    pr = urllib.parse.urlsplit(url)
    if pr.scheme != "https" or pr.hostname not in _ALLOWED:
        raise ValueError(f"非法请求目标: {pr.scheme}://{pr.hostname}")
    for info in socket.getaddrinfo(pr.hostname, 443):
        if not ipaddress.ip_address(info[4][0]).is_global:
            raise ValueError("解析到非公网地址")
    return urllib.request.build_opener(_NoRedirect).open(url, timeout=timeout).read()


EXCLUDE_KWS = ["mutant", "knockout", "knock-out", "overexpress", "crispr",
               "plasmid", "complement", "deletion mutant", "insertion mutant",
               "engineered", "edit", "whib7ko", "ko_tobramycin", "ko_"]


def _parse_exp(xml_text):
    """解析 expxml 字段 → dict(title, sample, bioproject, library, platform)."""
    def grab(tag):
        m = re.search(f"<{tag}>([^<]+)</{tag}>", xml_text)
        return m.group(1) if m else ""
    title = grab("Title")
    sample_m = re.search(r'<Sample[^>]*accession="([^"]+)"', xml_text)
    sample = sample_m.group(1) if sample_m else ""
    bp_m = re.search(r'accession="(PRJ[A-Z]+\d+)"', xml_text)
    bioproject = bp_m.group(1) if bp_m else ""
    library = grab("Library_strategy") or "RNA-Seq"
    layout = "PAIRED" if "paired" in xml_text.lower() else "SINGLE"
    inst_m = re.search(r'<Instrument>([^<]+)</Instrument>', xml_text)
    instrument = inst_m.group(1) if inst_m else ""
    return dict(title=title, sample=sample, bioproject=bioproject,
                library=library, layout=layout, instrument=instrument)


def _parse_runs(runs_xml):
    """解析 runs 字段 → [(run_acc, total_spots, total_bases_mb)]."""
    out = []
    for m in re.finditer(r'<Run\s+acc="(\w+)"\s+total_spots="(\d+)"\s+total_bases="(\d+)"', runs_xml):
        acc, spots, bases = m.group(1), int(m.group(2)), int(m.group(3))
        out.append((acc, spots, round(bases / 1e6, 1)))
    return out


def fetch_runs(taxid, retmax=500):
    """查询 + 解析 + 过滤。返回 (kept_rows, by_bioproject)。"""
    q = f"txid{taxid}[Organism] AND RNA-Seq[Strategy]"
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" +
           urllib.parse.urlencode({"db": "sra", "retmode": "json",
                                   "retmax": retmax, "term": q}))
    ids = json.loads(_get(url))["esearchresult"]["idlist"]
    rows = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        url2 = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" +
                urllib.parse.urlencode({"db": "sra", "retmode": "json",
                                        "id": ",".join(chunk)}))
        d = json.loads(_get(url2))["result"]
        for uid in d["uids"]:
            expxml = d[uid].get("expxml", "")
            runs_xml = d[uid].get("runs", "")
            exp = _parse_exp(expxml)
            for run_acc, spots, mb in _parse_runs(runs_xml):
                full_title = (exp["title"] + " " + run_acc).lower()
                if any(k in full_title for k in EXCLUDE_KWS):
                    continue
                rows.append(dict(run=run_acc, bioproject=exp["bioproject"],
                                 biosample=exp["sample"], layout=exp["layout"],
                                 spots=spots, mb=mb, title=exp["title"][:80],
                                 instrument=exp["instrument"]))
    # 按 BioProject 分组
    by_bp = {}
    for r in rows:
        by_bp.setdefault(r["bioproject"] or f"noBP_{r['run'][:6]}", []).append(r)
    return rows, by_bp
