# -*- coding: utf-8 -*-
"""Append a LysargiNase definition to a MaxQuant enzymes.xml file.

    python tools/add_enzyme.py "C:/path/to/MaxQuant/bin/conf/enzymes.xml"

MaxQuant does not ship LysargiNase (cleaves N-terminal to K/R); ORFscape's
mirror-digest support expects it when MaxQuant is used as a cross-check engine.
"""
import sys
from pathlib import Path


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip())
        return 2
    path = Path(argv[1]).expanduser()
    if not path.is_file():
        print("not found: %s" % path)
        return 2
    data = path.read_text(encoding="utf-8")
    if "LysargiNase" in data:
        print("already present")
        return 0
    aas = "ARNDCQEGHILKMFPSTWYV"
    spec = sorted({a + k for a in aas for k in "KR"})
    spec_lines = "\n".join("         <string>%s</string>" % s for s in spec)
    block = (
        '   <enzyme title="LysargiNase" description="Cleaves before lysine and '
        'arginine (N-terminal to KR)" create_date="2026-11-01T16:00:00" '
        'last_modified_date="2026-11-01T16:00:00" user="ORFscape">\n'
        "      <specificity>\n" + spec_lines + "\n      </specificity>\n   </enzyme>\n"
    )
    path.write_text(data.replace("</enzymes>", block + "</enzymes>"), encoding="utf-8")
    print("LysargiNase appended (%d specificity entries)" % len(spec))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
