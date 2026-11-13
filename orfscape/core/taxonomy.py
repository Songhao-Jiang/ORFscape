# -*- coding: utf-8 -*-
"""物种名解析：别名表 → 规范学名 →（在线时）NCBI taxid。

网络访问安全约束（防 SSRF）：协议仅 https、域名白名单、DNS 解析后阻断
私网/环回/链路本地地址、禁止重定向。任一校验失败按离线处理，不阻断主流程。
"""
import ipaddress
import json
import socket
import urllib.parse
import urllib.request

from ..config import SPECIES_ALIAS, TAX_SEARCH_ALIASES

_TAX_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ALLOWED_HOSTS = {"eutils.ncbi.nlm.nih.gov"}
_TAX_CACHE = {}


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise RuntimeError(f"禁止重定向: {code}")


def _guarded_urlopen(url, timeout):
    """SSRF 三重防护后的 urlopen：协议/域名白名单 + 公网 IP 校验 + 禁重定向。"""
    pr = urllib.parse.urlsplit(url)
    if pr.scheme != "https":
        raise ValueError("仅允许 https")
    if pr.hostname not in _ALLOWED_HOSTS:
        raise ValueError(f"域名不在白名单: {pr.hostname}")
    try:
        infos = socket.getaddrinfo(pr.hostname, 443)
    except OSError as e:
        raise ValueError(f"DNS 解析失败: {e}")
    for info in infos:
        ip = ipaddress.ip_address(info[4][0])
        if not ip.is_global:
            raise ValueError(f"解析到非公网地址: {ip}")
    opener = urllib.request.build_opener(_NoRedirect)
    return opener.open(url, timeout=timeout)


def resolve_species(name, online=True, timeout=8):
    """任意写法 → 规范学名 + taxid。

    返回 dict: {input, species, matched(alias|fuzzy|as_is), taxid, note}
    在线失败不影响主流程（taxid 置空，基因组获取阶段重试）。
    """
    key = name.strip()
    if key in _TAX_CACHE:
        return _TAX_CACHE[key]
    low = key.lower()
    if low in SPECIES_ALIAS:
        species, how = SPECIES_ALIAS[low], "alias"
    else:
        # 前缀模糊匹配（如 "myco abscess"）
        hits = {v for k, v in SPECIES_ALIAS.items() if low and (k.startswith(low) or low.startswith(k))}
        if len(hits) == 1:
            species, how = hits.pop(), "fuzzy"
        else:
            species, how = key, "as_is"
    taxid, note, hit_term = "", "", ""
    if online and species:
        # 依次尝试别名链（旧名/新名），首个返回结果者胜出
        for term in TAX_SEARCH_ALIASES.get(species, [species]):
            try:
                # 注：taxonomy 库不支持 [Scientific Name] 字段过滤（实测空返回），用裸名查询
                query = urllib.parse.urlencode({"db": "taxonomy", "retmode": "json",
                                                "term": term})
                url = f"{_TAX_BASE}?{query}"
                with _guarded_urlopen(url, timeout) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                ids = data.get("esearchresult", {}).get("idlist", [])
                if not ids:
                    continue
                taxid = ids[0]
                hit_term = term
                if len(ids) > 1:
                    note = f"多候选取首个（共{len(ids)}）"
                break
            except Exception as e:  # 含安全拦截，此别名按失败继续下一个
                note = f"taxid 在线解析失败: {type(e).__name__}"
        if not taxid and not note:
            note = "taxid 无命中（别名链全部为空）"
    out = {"input": key, "species": species, "matched": how, "taxid": taxid,
           "hit_term": hit_term, "note": note}
    _TAX_CACHE[key] = out
    return out
