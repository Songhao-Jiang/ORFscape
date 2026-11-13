# -*- coding: utf-8 -*-
"""对外请求的安全闸门 + 远程文件大小探测。

只做两件事：
1. :func:`validate_url` —— 仅 http/https、host 白名单、端口 443/80、
   解析出的每个地址都必须是公网地址（拒绝环回/私有/保留/链路本地/组播）。
2. :func:`head_size` —— 用 HEAD（必要时 Range: bytes=0-0）取
   Content-Length，给 GUI 展示"数据库下载大小"。

不跟随重定向（重定向可能指向内网）。任何失败返回 None，由界面显示 n/a。
"""
import http.client
import ipaddress
import socket
import urllib.error
import urllib.parse
import urllib.request

DEFAULT_PORTS = {80, 443}
TIMEOUT = 15


class UrlRejected(ValueError):
    """URL 未通过校验。"""


def _resolve_ok(host, port):
    """解析 host，任一地址不是公网地址即拒绝。"""
    try:
        infos = socket.getaddrinfo(host, port, proto=socket.IPPROTO_TCP)
    except socket.gaierror:
        return False
    if not infos:
        return False
    for info in infos:
        addr = info[4][0]
        try:
            ip = ipaddress.ip_address(addr)
        except ValueError:
            return False
        if (ip.is_private or ip.is_loopback or ip.is_reserved or ip.is_link_local
                or ip.is_multicast or ip.is_unspecified):
            return False
    return True


def validate_url(url, allowed_hosts, schemes=("https",)):
    """校验并返回 (host, port)。不合规抛 UrlRejected。"""
    if not isinstance(url, str) or not url.strip():
        raise UrlRejected("empty url")
    pr = urllib.parse.urlparse(url.strip())
    if pr.scheme not in schemes:
        raise UrlRejected("scheme not allowed: %s" % pr.scheme)
    host = (pr.hostname or "").lower()
    if not host:
        raise UrlRejected("missing host")
    allowed = {h.lower() for h in allowed_hosts}
    if host not in allowed:
        raise UrlRejected("host not allowed: %s" % host)
    port = pr.port or (443 if pr.scheme == "https" else 80)
    if port not in DEFAULT_PORTS:
        raise UrlRejected("port not allowed: %d" % port)
    if not _resolve_ok(host, port):
        raise UrlRejected("host resolves to a non-public address: %s" % host)
    return host, port


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise UrlRejected("redirect not followed: %s" % newurl)


_OPENER = urllib.request.build_opener(_NoRedirect)


def head_size(url, allowed_hosts, timeout=TIMEOUT):
    """返回字节数（int）或 None。

    优先 HEAD；服务端不给 Content-Length 时退化为 ``Range: bytes=0-0``
    读 Content-Range 的总长。
    """
    validate_url(url, allowed_hosts)
    req = urllib.request.Request(url, method="HEAD", headers={
        "User-Agent": "ORFscape/GUI (+size probe)"})
    try:
        with _OPENER.open(req, timeout=timeout) as resp:
            n = resp.headers.get("Content-Length")
            if n and n.isdigit():
                return int(n)
    except UrlRejected:
        raise
    except (urllib.error.URLError, http.client.HTTPException, OSError, ValueError):
        pass
    req = urllib.request.Request(url, headers={
        "User-Agent": "ORFscape/GUI (+size probe)", "Range": "bytes=0-0"})
    try:
        with _OPENER.open(req, timeout=timeout) as resp:
            cr = resp.headers.get("Content-Range") or ""
            if "/" in cr:
                tail = cr.rsplit("/", 1)[1].strip()
                if tail.isdigit():
                    return int(tail)
    except (UrlRejected, urllib.error.URLError, http.client.HTTPException,
            OSError, ValueError):
        return None
    return None
