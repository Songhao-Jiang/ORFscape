# -*- coding: utf-8 -*-
"""基因组自动获取（设计书 B 阶段）：物种 taxid/学名 → RefSeq 最优装配三件套下载与钉扎。

流程（split_test 实战验证过的路线固化）：
1. NCBI FTP 的物种 assembly_summary.txt 解析（表头以 '#assembly_accession' 开头）；
2. 选择 refseq_category=reference genome 优先、其次 representative、再次最新 Complete Genome；
3. 下载 <acc>_<name>_protein.faa.gz / _genomic.gff.gz / _genomic.fna.gz → 解压 → sha256 钉扎
   写 versions.json（accession/类别/下载日/各文件哈希）；
4. 离线回退：目标目录已有 versions.json 且文件齐全 → 直接复用。

网络限制：仅 https 且域名白名单 ftp.ncbi.nlm.nih.gov（复用 taxonomy 的防护风格）。
"""
import gzip
import hashlib
import ipaddress
import json
import os
import shutil
import socket
import urllib.parse
import urllib.request

_ALLOWED_HOSTS = {"ftp.ncbi.nlm.nih.gov"}
_SPECIES_DIR = {  # 学名 → RefSeq 物种目录名（改名物种用现行目录）
    "Mycobacterium abscessus": "Mycobacteroides_abscessus",
    "Mycobacterium chelonae": "Mycobacteroides_chelonae",
}


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise RuntimeError(f"禁止重定向: {code}")


def _guarded_urlopen(url, timeout=60):
    """SSRF 防护：https + 域名白名单 + DNS 公网校验 + 禁重定向。"""
    pr = urllib.parse.urlsplit(url)
    if pr.scheme != "https":
        raise ValueError("仅允许 https")
    if pr.hostname not in _ALLOWED_HOSTS:
        raise ValueError(f"域名不在白名单: {pr.hostname}")
    for info in socket.getaddrinfo(pr.hostname, 443):
        if not ipaddress.ip_address(info[4][0]).is_global:
            raise ValueError("解析到非公网地址")
    return urllib.request.build_opener(_NoRedirect).open(url, timeout=timeout)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pick_assembly(species):
    """解析物种 assembly_summary → (accession, refseq_category, organism, ftp_path)。"""
    sp_dir = _SPECIES_DIR.get(species, species.replace(" ", "_"))
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/{sp_dir}/assembly_summary.txt"
    with _guarded_urlopen(url) as resp:
        txt = resp.read().decode("utf-8")
    import csv
    import io
    lines = [l for l in txt.splitlines() if not l.startswith("##")]
    hdr = lines[0].lstrip("#")
    rows = list(csv.DictReader(io.StringIO("\n".join([hdr] + lines[1:])), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"{sp_dir} 无 RefSeq 装配")
    def rank(r):
        cat = r.get("refseq_category", "na")
        cat_rank = {"reference genome": 0, "representative genome": 1}.get(cat, 2)
        complete = 0 if r.get("assembly_level") == "Complete Genome" else 1
        return (cat_rank, complete)
    best = min(rows, key=rank)
    return best["assembly_accession"], best.get("refseq_category", "na"), \
        best.get("organism_name", ""), best["ftp_path"]


def fetch_genome(species, out_dir, today=""):
    """下载三件套并钉扎。已有齐备 versions.json 时离线复用。返回 versions dict。"""
    os.makedirs(out_dir, exist_ok=True)
    ver_path = os.path.join(out_dir, "versions.json")
    need = ["protein.faa", "genomic.gff", "genomic.fna"]
    if os.path.exists(ver_path):
        meta = json.load(open(ver_path, encoding="utf-8"))
        if meta.get("species") == species and all(
                os.path.exists(os.path.join(out_dir, f)) for f in need):
            return meta  # 离线复用
    acc, cat, org, ftp = pick_assembly(species)
    base = ftp.rstrip("/")
    stem = base.rsplit("/", 1)[-1]
    meta = {"species": species, "accession": acc, "refseq_category": cat,
            "organism": org, "ftp_path": ftp, "downloaded": today, "files": {}}
    for out_name, suffix in (("protein.faa", "_protein.faa.gz"),
                             ("genomic.gff", "_genomic.gff.gz"),
                             ("genomic.fna", "_genomic.fna.gz")):
        url = f"{base}/{stem}{suffix}"
        gz_name = stem + suffix
        with _guarded_urlopen(url) as resp, open(os.path.join(out_dir, gz_name), "wb") as fh:
            shutil.copyfileobj(resp, fh)
        gz_path = os.path.join(out_dir, gz_name)
        meta["files"][out_name] = {"sha256_gz": _sha256(gz_path),
                                   "size": os.path.getsize(gz_path)}
        with gzip.open(gz_path, "rb") as fi, \
                open(os.path.join(out_dir, out_name), "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(gz_path)
    n_prot = sum(1 for l in open(os.path.join(out_dir, "protein.faa"), encoding="utf-8")
                 if l.startswith(">"))
    meta["n_proteins"] = n_prot
    with open(ver_path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    return meta
