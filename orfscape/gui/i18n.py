# -*- coding: utf-8 -*-
"""GUI 文案（英文默认，可切中文）。

约定
----
- 键名点分：``区域.元素``（``pipe.eta``、``st.running``、``stage.E2``）。
- 值统一为 ``(en, zh)`` 二元组；英文是交付语言，中文只做界面切换。
- 模板用 ``str.format`` 占位：``"Remaining {t}"``。
- 缺键时返回键本身并加 ``!`` 前缀，便于一眼看出漏翻。
"""
LANG_EN = "en"
LANG_ZH = "zh"
LANGS = [(LANG_EN, "English"), (LANG_ZH, "中文")]

CATALOG = {
    # ---------------------------------------------------------------- 窗口
    "app.title": ("ORFscape — Proteogenomics", "ORFscape — 蛋白质基因组学"),
    "app.subtitle": ("raw → layered database → novel peptide → arbitration → report",
                     "raw → 层级库 → 新肽 → 仲裁 → 报告"),
    "app.lang": ("Language", "语言"),
    "theme.to_dark": ("Dark mode", "深色"),
    "theme.to_light": ("Light mode", "浅色"),

    # ---------------------------------------------------------------- 区块
    "sec.input": ("Input", "输入"),
    "sec.pipeline": ("Pipeline", "流程"),
    "sec.resources": ("Resources", "运行资源"),
    "sec.database": ("Databases", "数据库"),
    "sec.log": ("Log", "日志"),

    # ---------------------------------------------------------------- 输入
    "inp.raw": ("RAW files", "RAW 文件"),
    "inp.add_files": ("Add files…", "添加文件…"),
    "inp.add_folder": ("Add folder…", "添加目录…"),
    "inp.clear": ("Clear", "清空"),
    "inp.remove": ("Remove", "移除"),
    "inp.no_raw": ("No RAW selected", "未选择 RAW"),
    "inp.summary": ("{n} file(s) · {size}", "{n} 个文件 · {size}"),
    "inp.summary_n": ("{n} file(s)", "{n} 个文件"),
    "inp.species": ("Species", "物种"),
    "inp.species_hint": ("Latin name or code; blank = infer from filenames",
                         "拉丁名或代码；留空按文件名推断"),
    "inp.desc": ("Description", "描述"),
    "inp.desc_hint": ("optional note stored in the run log",
                      "可选备注，写入运行日志"),
    "inp.workdir": ("Project folder", "项目目录"),
    "inp.browse": ("Browse…", "浏览…"),
    "inp.threads": ("pFind threads", "pFind 线程"),
    "inp.workers": ("Re-search workers", "回搜并发"),
    "inp.dbs": ("Discovery DBs", "发现库"),
    "inp.font_scale": ("Font size", "字号"),
    "inp.skip_p6b": ("Skip candidate re-search (P6b)", "跳过候选回搜（P6b）"),
    "inp.skip_qc": ("Skip QC", "跳过 QC"),
    "inp.confirm": ("Confirm before start", "开始前确认"),
    "inp.force": ("Re-run completed stages", "重跑已完成阶段"),
    "inp.advanced": ("Advanced", "高级"),
    "inp.attach": ("Load log", "载入日志"),
    "inp.loaded": ("Loaded {n} event(s) from {path}",
                   "已载入 {n} 条事件：{path}"),
    "inp.notfound": ("No progress log in that project", "该项目下没有进度日志"),

    # ---------------------------------------------------------------- 流程
    "pipe.overall": ("Overall progress", "总进度"),
    "pipe.elapsed": ("Elapsed", "已用"),
    "pipe.eta": ("Remaining", "剩余"),
    "pipe.eta_unknown": ("estimating…", "估算中…"),
    "pipe.current": ("Current", "当前"),
    "pipe.idle": ("Idle — choose inputs and press Start", "空闲 — 选好输入后点开始"),
    "pipe.empty": ("Add a RAW folder or load an existing project",
                   "添加 RAW 目录，或载入已有项目"),
    "pipe.stage_n": ("stage {i}/{n}", "第 {i}/{n} 步"),
    "pipe.col.stage": ("Stage", "阶段"),
    "pipe.col.status": ("Status", "状态"),
    "pipe.col.progress": ("Progress", "进度"),
    "pipe.col.elapsed": ("Elapsed", "耗时"),
    "pipe.col.eta": ("ETA", "预计"),
    "pipe.col.check": ("Check", "自检"),
    "pipe.step": ("{done}/{total} {unit}", "{done}/{total} {unit}"),
    "pipe.psm": ("PSM {n}", "谱 {n}"),
    "pipe.stages": ("{done}/{total} stages", "{done}/{total} 阶段"),

    # ---------------------------------------------------------------- 状态
    "st.pending": ("Pending", "待跑"),
    "st.running": ("Running", "运行中"),
    "st.ok": ("Done", "完成"),
    "st.skip": ("Skipped", "跳过"),
    "st.fail": ("Failed", "失败"),
    "st.check-fail": ("Check failed", "自检失败"),
    "st.stopped": ("Stopped", "已停止"),
    "st.aborted": ("Not run", "未运行"),
    "st.warn": ("Warning", "警告"),

    # ---------------------------------------------------------------- 阶段名
    "stage.P0": ("Data confirmation", "数据确定"),
    "stage.P1": ("Genome fetch", "基因组获取"),
    "stage.D0": ("Spectra / CV split", "谱图/电压拆分"),
    "stage.C": ("Layer database", "层级建库"),
    "stage.E1": ("Annotated search", "已注释库搜库"),
    "stage.QC": ("QC", "质控"),
    "stage.E2": ("Discovery search", "发现搜库"),
    "stage.P6": ("Novel peptides", "新肽筛选"),
    "stage.P6b": ("Candidate re-search", "候选回搜"),
    "stage.P7": ("isoPSV arbitration", "isoPSV 仲裁"),
    "stage.P8a": ("Event classification", "事件分类"),
    "stage.P8b": ("Functional annotation", "功能注释"),
    "stage.P8c": ("Report", "报告"),

    # 看板用短名（宽度只够 12 字符）
    "stage_short.P0": ("Confirm", "数据确定"),
    "stage_short.P1": ("Genome", "基因组"),
    "stage_short.D0": ("Spectra", "谱图"),
    "stage_short.C": ("Database", "建库"),
    "stage_short.E1": ("Annotated", "已注释"),
    "stage_short.QC": ("QC", "质控"),
    "stage_short.E2": ("Discovery", "发现"),
    "stage_short.P6": ("Novel", "新肽"),
    "stage_short.P6b": ("Re-search", "回搜"),
    "stage_short.P7": ("isoPSV", "仲裁"),
    "stage_short.P8a": ("Events", "事件"),
    "stage_short.P8b": ("Annotation", "注释"),
    "stage_short.P8c": ("Report", "报告"),

    # ---------------------------------------------------------------- 资源
    "res.cpu": ("CPU", "CPU"),
    "res.mem": ("System memory", "内存"),
    "res.disk": ("Project disk", "项目磁盘"),
    "res.rss": ("Pipeline RSS", "流程占用"),
    "res.cores": ("{n} cores", "{n} 核"),
    "res.mem_val": ("{used} / {total} GB", "{used} / {total} GB"),
    "res.disk_val": ("{free} free of {total}", "剩余 {free} / 共 {total}"),
    "res.na": ("n/a", "无数据"),

    # ---------------------------------------------------------------- 估算
    "est.title": ("Estimate", "预估"),
    "est.raw": ("RAW input", "RAW 规模"),
    "est.ms2": ("MS2 spectra", "MS2 谱数"),
    "est.genome": ("Genome", "基因组"),
    "est.tasks": ("Tasks (E1 / E2)", "任务（E1 / E2）"),
    "est.time": ("Total time", "总耗时"),
    "est.disk": ("Disk required", "所需磁盘"),
    "est.ram": ("Peak RAM", "峰值内存"),
    "est.src.measured": ("measured", "实测"),
    "est.src.estimated": ("estimated", "估算"),
    "est.src.assumed": ("assumed", "假设"),
    "est.src.calibrated": ("calibrated", "已标定"),
    "est.value_src": ("{v} ({src})", "{v}（{src}）"),
    "est.disk_val": ("{need} · {free} free", "{need} · 剩余 {free}"),
    "est.recalc": ("Recalculate", "重新估算"),
    "est.unknown": ("—", "—"),

    # ---------------------------------------------------------------- 数据库
    "db.col.name": ("Database", "数据库"),
    "db.col.size": ("Size", "大小"),
    "db.col.id": ("ID", "编号"),
    "db.col.state": ("Status", "状态"),
    "db.total": ("Total", "合计"),
    "db.open": ("Databases…", "数据库…"),
    "db.close": ("Close", "关闭"),
    "db.refresh": ("Refresh", "刷新"),
    "db.copy_path": ("Copy path", "复制路径"),
    "db.probe": ("Probe remote sizes", "探测远程大小"),
    "db.probing": ("Probing…", "探测中…"),
    "db.hint": ("Local sizes are exact; RefSeq uses HEAD requests; UniProt is "
                "estimated from protein count.", "本地为实测；RefSeq 用 HEAD 探测；"
                "UniProt 按蛋白条数估算。"),
    "db.state.present": ("present", "已就绪"),
    "db.state.pending": ("to download", "待下载"),
    "db.state.remote": ("remote", "远程"),
    "db.state.estimate": ("estimated", "估算"),
    "db.state.unknown": ("unknown", "未知"),

    # ---------------------------------------------------------------- 日志
    "log.level": ("Level", "级别"),
    "log.filter": ("Filter", "过滤"),
    "log.autoscroll": ("Auto-scroll", "自动滚动"),
    "log.clear": ("Clear", "清空"),
    "log.empty": ("No output yet.", "暂无输出。"),
    "log.save": ("Save log…", "导出日志…"),
    "log.all": ("All", "全部"),
    "log.warn_plus": ("Warnings+", "警告以上"),
    "log.errors": ("Errors only", "仅错误"),
    "log.find": ("Find", "查找"),
    "log.find_prev": ("Prev", "上一个"),
    "log.find_next": ("Next", "下一个"),
    "menu.copy": ("Copy", "复制"),
    "menu.copy_all": ("Copy all", "复制全部"),
    "menu.filter": ("Filter by this text", "按此文本过滤"),
    "menu.open_log": ("Open log file", "打开日志文件"),

    # ---------------------------------------------------------------- 按钮
    "btn.start": ("Start", "开始"),
    "btn.resume": ("Resume", "续跑"),
    "btn.stop": ("Stop", "停止"),
    "btn.open": ("Open project folder", "打开项目目录"),
    "btn.results": ("Open results", "打开结果"),
    "btn.export": ("Export summary", "导出摘要"),
    "btn.cancel": ("Cancel", "取消"),

    # ---------------------------------------------------------------- 提示
    "msg.running": ("Running · {stage}", "运行中 · {stage}"),
    "msg.finished": ("Finished in {t}", "完成，用时 {t}"),
    "msg.failed": ("Failed at {stage}", "在 {stage} 失败"),
    "msg.checkfail": ("Stage check failed at {stage}", "{stage} 自检未通过"),
    "msg.stopped": ("Stopped after {t}", "已停止，用时 {t}"),
    "msg.stop_req": ("Stopping…", "正在停止…"),
    "msg.need_raw": ("Select at least one RAW file or folder", "请至少选择一个 RAW 文件或目录"),
    "msg.no_raw_found": ("No .raw file found under the selection", "所选路径下没有 .raw 文件"),
    "msg.resync": ("Progress log was rewritten — state reloaded",
                   "进度日志被重写，状态已重新加载"),
    "msg.need_dir": ("Choose a project folder", "请选择项目目录"),
    "msg.bad_dir": ("Project folder does not exist", "项目目录不存在"),
    "msg.started": ("Started · {t}", "已启动 · {t}"),
    "msg.attached": ("Attached to running log", "已接管运行中的日志"),
    "msg.done_prev": ("Loaded previous run ({n} events)", "已载入上一轮运行（{n} 条事件）"),
    "msg.no_log": ("—", "—"),
    "msg.exported": ("Saved {path}", "已保存 {path}"),
    "msg.copied": ("Copied {what}", "已复制 {what}"),
    "tip.status": ("Status: {v}", "状态：{v}"),
    "tip.progress": ("Progress: {v}", "进度：{v}"),
    "tip.elapsed": ("Elapsed {v} · est {e}", "已用 {v} · 预计 {e}"),
    "tip.check": ("Check {f} fail / {w} warn", "自检 {f} 失败 / {w} 警告"),
    "plan.title": ("Run plan", "运行计划"),
    "plan.hint": ("Stages already completed are skipped automatically.",
                  "已完成的阶段会自动跳过。"),
    "plan.est": ("Estimated {t} · disk {d} · peak RAM {r}",
                 "预计 {t} · 磁盘 {d} · 峰值内存 {r}"),
    "plan.skip": ("skip", "跳过"),
    "plan.run": ("run", "执行"),

    # ---------------------------------------------------------------- 告警
    "warn.disk": ("Not enough free disk: need {need}, {free} free",
                  "磁盘不足：需 {need}，剩余 {free}"),
    "warn.ram": ("Estimated peak RAM {need} GB exceeds available memory ({total} GB)",
                 "预计峰值内存 {need} GB 超过可用内存（{total} GB）"),
    "warn.threads": ("{threads} threads requested, machine has {cpu} cores",
                     "请求 {threads} 线程，机器只有 {cpu} 核"),
    "warn.noraw": ("No RAW file selected", "未选择 RAW 文件"),
}


def tr(key, lang=LANG_EN, **kw):
    """取文案；缺键返回 ``!key``。"""
    pair = CATALOG.get(key)
    if not pair:
        return "!" + key
    text = pair[1] if lang == LANG_ZH else pair[0]
    if kw:
        try:
            return text.format(**kw)
        except (KeyError, IndexError, ValueError):
            return text
    return text


def missing_keys(used):
    """测试用：返回 used 里没有定义的键。"""
    return sorted(k for k in used if k not in CATALOG)


def all_keys():
    return sorted(CATALOG)
