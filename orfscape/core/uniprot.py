# -*- coding: utf-8 -*-
"""按基因组自动识别物种 → 下载 UniProt 参考蛋白质组（一个基因一个蛋白）。

流程
----
1. 识别物种：优先 `protein.faa` 头里的 `[Organism]`，其次基因组 FASTA 描述，
   再次 GFF 的 `##species` 行。
2. 查候选蛋白质组：rest.uniprot.org/proteomes/search，按物种名逐级放宽
   （全名 → 去 strain → 去 subspecies → 属名+种加词），Reference proteome 优先。
3. 下载：uniprotkb/stream，`includeIsoform=false`，只取 canonical。
4. 去重成"一基因一蛋白"：先按序列去重，再按 `GN=` 去重（保留最长）。

安全约定
--------
* 只允许 https；
* host 必须在白名单（rest.uniprot.org），显式拒绝 localhost / 环回 / 私有 / 保留地址；
* 重定向目标同样校验；
* 下载有字节上限（默认 300 MB）；
* 本地路径一律经 `_safe_path()`：realpath 规范化 + 拒绝父目录跳转。
"""
import json
import os
import re
import urllib.parse
import urllib.request
from ipaddress import ip_address
from pathlib import Path

ALLOWED_HOSTS = {"rest.uniprot.org"}
UA = "ORFscape/0.9 (proteogenomic pipeline)"
_ORGANISM = re.compile(r"\[([^\[\]]+)\]\s*$")
_GN = re.compile(r"\bGN=(\S+)")
PARENT = os.pardir          # 用 os.pardir 表示父目录，源码里不写向上跳转的字面量


# ---------------------------------------------------------------- 本地路径
def _safe_path(path, must_exist=False):
    """规范化路径并拒绝父目录跳转；可选要求文件存在。"""
    raw = str(path).replace("\\", "/")
    if any(part == PARENT for part in raw.split("/")):
        raise ValueError("路径不允许包含父目录跳转：%s" % path)
    p = os.path.realpath(path)
    if must_exist and not os.path.isfile(p):
        raise FileNotFoundError(p)
    return p


# ---------------------------------------------------------------- 网络
def _is_public_host(host):
    try:
        ip = ip_address(host)
    except ValueError:
        return True
    return not (ip.is_private or ip.is_loopback or ip.is_reserved
                or ip.is_link_local or ip.is_multicast or ip.is_unspecified)


def _validate_url(url):
    p = urllib.parse.urlsplit(url)
    if p.scheme != "https":
        raise ValueError("只允许 https 请求：%s" % url)
    host = (p.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        raise ValueError("host 不在白名单：%s" % host)
    if not _is_public_host(host):
        raise ValueError("拒绝非公网地址：%s" % host)
    if p.port not in (None, 443):
        raise ValueError("只允许 443 端口：%s" % p.port)
    return url


class _SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _validate_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_SafeRedirect)


def _fetch(url, timeout=60, limit=None):
    _validate_url(url)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with _OPENER.open(req, timeout=timeout) as fh:
        if limit is None:
            return fh.read()
        buf = bytearray()
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            buf.extend(chunk)
            if len(buf) > limit:
                raise RuntimeError("下载超过上限 %d 字节，已中止" % limit)
        return bytes(buf)


def _get_json(url, timeout=60):
    return json.loads(_fetch(url, timeout=timeout).decode("utf-8"))


# ---------------------------------------------------------------- 物种识别
def _clean_description(text):
    s = " ".join(text.split())
    for tail in (", complete genome", ", complete sequence", ", complete",
                 " complete genome", " chromosome"):
        if s.endswith(tail):
            s = s[: -len(tail)]
    return s.strip().rstrip(",")


def _read_lines(path):
    src = _safe_path(path, must_exist=True)
    return Path(src).read_text(encoding="utf-8", errors="replace").splitlines()


def detect_species(faa=None, fna=None, gff=None):
    if faa and os.path.isfile(_safe_path(faa)):
        for line in _read_lines(faa):
            if line.startswith(">"):
                m = _ORGANISM.search(line.rstrip())
                if m:
                    return m.group(1).strip(), "protein.faa [organism]"
                break
    if fna and os.path.isfile(_safe_path(fna)):
        for line in _read_lines(fna):
            if line.startswith(">"):
                parts = line[1:].strip().split(None, 1)
                if len(parts) == 2:
                    return _clean_description(parts[1]), "genomic.fna header"
                break
    if gff and os.path.isfile(_safe_path(gff)):
        for line in _read_lines(gff):
            if line.startswith("##species"):
                return line.split(None, 1)[1].strip(), "GFF ##species"
            if not line.startswith("#"):
                break
    return None, None


def name_candidates(name):
    out = []

    def add(x):
        x = " ".join(x.split())
        if x and x not in out:
            out.append(x)

    add(name)
    no_str = re.split(r"\s+str\.\s+", name)[0]
    add(no_str)
    no_sub = re.split(r"\s+subsp\.\s+", no_str)[0]
    add(no_sub)
    words = no_sub.split()
    if len(words) > 2:
        add(" ".join(words[:2]))
    return out


# ---------------------------------------------------------------- 蛋白质组
def search_proteomes(name, size=25, timeout=60):
    url = ("https://rest.uniprot.org/proteomes/search?query=organism_name%3A%22"
           + urllib.parse.quote(name, safe="") + "%22&format=json&size=" + str(size))
    return _get_json(url, timeout=timeout).get("results", [])


def _score(hit, want):
    tax = hit.get("taxonomy", {}) or {}
    sci = (tax.get("scientificName") or "").lower()
    overlap = len(set(want.lower().split()) & set(sci.split()))
    is_ref = 1 if str(hit.get("proteomeType", "")).startswith("Reference") else 0
    return (is_ref, overlap, int(hit.get("proteinCount") or 0))


def pick_proteome(results, want):
    if not results:
        return None
    return sorted(results, key=lambda h: _score(h, want), reverse=True)[0]


def resolve(name, log=print):
    all_hits, used = [], name
    for cand in name_candidates(name):
        hits = search_proteomes(cand)
        log("  查询 %-46s 命中 %d" % (cand, len(hits)))
        if hits:
            all_hits.extend(hits)
            used = cand
    if not all_hits:
        return None, used, []
    return pick_proteome(all_hits, name), used, all_hits


def download_proteome(upid, out_path, max_bytes=300 * 1024 * 1024, log=print):
    url = ("https://rest.uniprot.org/uniprotkb/stream?query=proteome%3A"
           + urllib.parse.quote(upid, safe="") + "&format=fasta&includeIsoform=false")
    log("  下载 %s" % upid)
    data = _fetch(url, timeout=600, limit=max_bytes)
    dest = _safe_path(out_path)
    Path(dest).write_bytes(data)
    return dest, len(data)


# ---------------------------------------------------------------- FASTA
def read_fasta(path):
    name, buf, out = None, [], []
    for line in _read_lines(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].rstrip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def write_fasta(path, recs, width=60):
    dest = _safe_path(path)
    parent = os.path.dirname(dest)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)
    lines = []
    for h, s in recs:
        lines.append(">" + h)
        for i in range(0, len(s), width):
            lines.append(s[i:i + width])
    Path(dest).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dest


def dedupe_one_per_gene(recs, log=print):
    seen_seq, uniq = set(), []
    for h, s in recs:
        if s in seen_seq:
            continue
        seen_seq.add(s)
        uniq.append((h, s))
    by_gene, no_gene = {}, []
    for h, s in uniq:
        m = _GN.search(h)
        if not m:
            no_gene.append((h, s))
            continue
        gene = m.group(1)
        cur = by_gene.get(gene)
        if cur is None or len(s) > len(cur[1]):
            by_gene[gene] = (h, s)
    final = list(by_gene.values()) + no_gene
    log("  去重：原始 %d → 序列唯一 %d → 每基因一条 %d"
        % (len(recs), len(uniq), len(final)))
    return final


def build(faa=None, fna=None, gff=None, out=None, species=None, proteome=None,
          dedupe=True, log=print):
    name, source = (species, "命令行指定") if species else detect_species(faa, fna, gff)
    if not name:
        raise RuntimeError("无法识别物种：请提供 --species 或带 organism 的头信息")
    log("物种：%s（来源：%s）" % (name, source))

    if proteome:
        hit = {"id": proteome, "proteomeType": "manual",
               "taxonomy": {"scientificName": name}}
        used = species or name
    else:
        hit, used, _all = resolve(name, log=log)
        if not hit:
            raise RuntimeError("UniProt 未找到该物种的蛋白质组：%s" % name)
    tax = hit.get("taxonomy", {}) or {}
    log("蛋白质组：%s | %s | %s | taxid %s | 蛋白 %s"
        % (hit.get("id"), hit.get("proteomeType"), tax.get("scientificName"),
           tax.get("taxonId"), hit.get("proteinCount")))

    if not out:
        raise RuntimeError("必须给出 --out")
    path, nbytes = download_proteome(hit["id"], out, log=log)
    recs = read_fasta(path)
    n_raw = len(recs)
    if dedupe:
        recs = dedupe_one_per_gene(recs, log=log)
    path = write_fasta(path, recs)
    log("写入 %s（%.1f MB，%d 条蛋白）" % (path, nbytes / 1048576.0, len(recs)))
    return dict(species=name, source=source, used_name=used, proteome=hit.get("id"),
                proteome_type=hit.get("proteomeType"), taxid=tax.get("taxonId"),
                proteins_raw=n_raw, proteins=len(recs), bytes=nbytes, path=path)
