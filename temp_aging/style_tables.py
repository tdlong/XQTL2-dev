#!/usr/bin/env python3
"""style_tables.py FILE.docx [target_twips]

Pandoc hardcodes docx tables to 7920 twips (5.5 in) whatever the reference
document's page size, so every column is squeezed and cells wrap. This rescales
each table to the real text width (default 9360 twips = 6.5 in at 1 in margins),
preserving pandoc's relative widths, and adds what pandoc omits: cell padding so
columns do not touch, and rules above the table, under the header and below the
last row.
"""
import re, sys, zipfile, shutil, os

path = sys.argv[1]
target = int(sys.argv[2]) if len(sys.argv) > 2 else 9360

def widen(xml):
    out, pos, n = [], 0, 0
    for m in re.finditer(r"<w:tbl>.*?</w:tbl>", xml, re.S):
        tbl = m.group(0)
        cols = [int(x) for x in re.findall(r'<w:gridCol w:w="(\d+)"', tbl)]
        if cols and sum(cols):
            f = target / sum(cols)
            new = tbl
            new = re.sub(r'(<w:gridCol w:w=")(\d+)(")',
                         lambda g: g.group(1) + str(round(int(g.group(2)) * f)) + g.group(3), new)
            new = re.sub(r'(<w:tcW w:w=")(\d+)(")',
                         lambda g: g.group(1) + str(round(int(g.group(2)) * f)) + g.group(3), new)
            new = re.sub(r'<w:tblW w:w="\d+" w:type="\w+"\s*/>',
                         '<w:tblW w:w="%d" w:type="dxa"/>' % target, new)
            # rules top / under header / bottom, and padding so columns breathe
            borders = ('<w:tblBorders>'
                       '<w:top w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
                       '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
                       '<w:insideH w:val="nil"/><w:insideV w:val="nil"/>'
                       '<w:left w:val="nil"/><w:right w:val="nil"/>'
                       '</w:tblBorders>'
                       '<w:tblCellMar>'
                       '<w:left w:w="108" w:type="dxa"/><w:right w:w="108" w:type="dxa"/>'
                       '<w:top w:w="40" w:type="dxa"/><w:bottom w:w="40" w:type="dxa"/>'
                       '</w:tblCellMar>')
            new = re.sub(r'(<w:tblPr>)', r'\1' + borders, new, count=1)
            # underline the header row
            first = re.search(r'<w:tr\b.*?</w:tr>', new, re.S)
            if first:
                hdr = first.group(0)
                hdr2 = hdr.replace('<w:tcPr>',
                    '<w:tcPr><w:tcBorders>'
                    '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
                    '</w:tcBorders>')
                new = new.replace(hdr, hdr2, 1)
            out.append(xml[pos:m.start()]); out.append(new); pos = m.end(); n += 1
    out.append(xml[pos:])
    return "".join(out), n

tmp = path + ".tmp"
zin = zipfile.ZipFile(path)
zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
total = 0
for it in zin.infolist():
    data = zin.read(it.filename)
    if it.filename == "word/document.xml":
        t, total = widen(data.decode("utf-8"))
        data = t.encode("utf-8")
    zout.writestr(it, data)
zout.close(); zin.close()
os.replace(tmp, path)
print("  widened %d tables to %.2f in" % (total, target / 1440))
