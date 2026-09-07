#!/usr/bin/env python3
"""fix_math.py FILE.docx

Repair pandoc's LaTeX -> OMML accent conversion.

  \\bar{x} becomes <m:acc> carrying U+203E OVERLINE as an accent character, so
  the bar renders offset to one side instead of centred over the symbol. The
  correct construct is <m:bar m:pos="top">.

  \\overline{...} becomes <m:bar> but keeps a stray m:chr (a circumflex, from
  whatever was nested inside), which Word has no meaning for.
"""
import re, sys, zipfile, os

OVERLINE = "‾"
path = sys.argv[1]

def fix(xml):
    n_acc = n_bar = 0

    def acc_to_bar(m):
        nonlocal n_acc
        blk = m.group(0)
        chr_ = re.search(r'<m:chr m:val="([^"]*)"\s*/>', blk)
        if not chr_ or chr_.group(1) != OVERLINE:
            return blk
        n_acc += 1
        blk = blk.replace('<m:acc>', '<m:bar>').replace('</m:acc>', '</m:bar>')
        blk = blk.replace('<m:accPr>', '<m:barPr>').replace('</m:accPr>', '</m:barPr>')
        blk = re.sub(r'<m:chr m:val="[^"]*"\s*/>', '<m:pos m:val="top"/>', blk, count=1)
        return blk

    xml = re.sub(r'<m:acc>.*?</m:acc>', acc_to_bar, xml, flags=re.S)

    def clean_bar(m):
        nonlocal n_bar
        blk = m.group(0)
        if '<m:chr' not in blk:
            return blk
        n_bar += 1
        return re.sub(r'<m:chr m:val="[^"]*"\s*/>', '', blk, count=1)

    xml = re.sub(r'<m:bar>.*?</m:bar>', clean_bar, xml, flags=re.S)
    return xml, n_acc, n_bar

tmp = path + ".tmp"
zin  = zipfile.ZipFile(path)
zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
a = b = 0
for it in zin.infolist():
    data = zin.read(it.filename)
    if it.filename == "word/document.xml":
        t, a, b = fix(data.decode("utf-8"))
        data = t.encode("utf-8")
    zout.writestr(it, data)
zout.close(); zin.close(); os.replace(tmp, path)
print("  %d offset bars corrected, %d stray accent chars removed" % (a, b))
