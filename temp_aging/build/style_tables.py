#!/usr/bin/env python3
"""style_tables.py FILE.docx [target_twips]

Turn pandoc's docx tables into scientific tables, following the booktabs
conventions (Fear, _Publication quality tables in LaTeX_):

  - never a vertical rule
  - never a double rule
  - rules of graded weight: heavy above and below the table, lighter under the
    column headers, lighter still over a column group
  - a group rule spans only the columns it refers to
  - generous space above and below each rule

and two things pandoc gets wrong for docx output:

  - it hardcodes tables to 7920 twips (5.5 in) whatever the reference document's
    page size, so cells wrap; rescaled here to the real text width
  - it emits no cell padding, so columns touch

GROUP HEADERS. pandoc allows one header row, so put the group labels there and
the column labels in the first body row:

    |      |     | h^2  |      |        |
    |:-----|----:|-----:|-----:|-------:|
    | Design | n | 2%   | 1%   | 0.5%   |

A non-empty cell in the group row spans forward over the empty cells after it
and takes a rule over exactly that span; leading empty cells stay bare.
"""
import re, sys, zipfile, os

path   = sys.argv[1]
target = int(sys.argv[2]) if len(sys.argv) > 2 else 9360

# w:sz is eighths of a point. booktabs: heavy .08em, light .05em, cmidrule .03em
HEAVY, LIGHT, GROUP = 12, 8, 6

def rule(sz):
    return '<w:bottom w:val="single" w:sz="%d" w:space="0" w:color="000000"/>' % sz

BORDERS = ('<w:tblBorders>'
           '<w:top w:val="single" w:sz="%d" w:space="0" w:color="000000"/>'
           '<w:bottom w:val="single" w:sz="%d" w:space="0" w:color="000000"/>'
           '<w:insideH w:val="nil"/><w:insideV w:val="nil"/>'
           '<w:left w:val="nil"/><w:right w:val="nil"/>'
           '</w:tblBorders>') % (HEAVY, HEAVY)
CELLMAR = ('<w:tblCellMar>'
           '<w:left w:w="108" w:type="dxa"/><w:right w:w="108" w:type="dxa"/>'
           '<w:top w:w="60" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/>'
           '</w:tblCellMar>')

cells   = lambda tr: re.findall(r'<w:tc>.*?</w:tc>', tr, re.S)
# <m:t> as well as <w:t>: a header that is an equation, such as $h^2$, puts its
# text in the math run and would otherwise read as an empty cell
text_of = lambda tc: "".join(re.findall(r'<[wm]:t[^>]*>([^<]*)</[wm]:t>', tc)).strip()

def with_props(tc, extra):
    """pandoc emits no <w:tcPr>, so insert one rather than assume it exists."""
    if '<w:tcPr>' in tc:
        return tc.replace('<w:tcPr>', '<w:tcPr>' + extra, 1)
    return tc.replace('<w:tc>', '<w:tc><w:tcPr>' + extra + '</w:tcPr>', 1)

def centre(tc):
    """Centre every paragraph in a cell, inserting <w:pPr> where absent."""
    def fix(m):
        p = m.group(0)
        if '<w:pPr>' in p:
            return p.replace('<w:pPr>', '<w:pPr><w:jc w:val="center"/>', 1)
        return p.replace('<w:p>', '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>', 1)
    return re.sub(r'<w:p\b.*?</w:p>', fix, tc, flags=re.S)

def rebuild(tr, new_cells):
    return re.sub(r'(<w:tr\b[^>]*>).*(</w:tr>)',
                  lambda m: m.group(1) + "".join(new_cells) + m.group(2), tr, flags=re.S)

def style(tbl):
    rows = re.findall(r'<w:tr\b.*?</w:tr>', tbl, re.S)
    if not rows:
        return tbl
    hdr, cs = rows[0], cells(rows[0])
    txt = [text_of(c) for c in cs]

    if any(t == "" for t in txt) and any(t for t in txt):     # group-header row
        out, i = [], 0
        while i < len(cs):
            if txt[i] == "":
                out.append(cs[i]); i += 1; continue
            j = i + 1
            while j < len(cs) and txt[j] == "":
                j += 1
            span  = j - i
            extra = ('<w:gridSpan w:val="%d"/>' % span if span > 1 else '') \
                    + '<w:tcBorders>' + rule(GROUP) + '</w:tcBorders>'
            cell = with_props(cs[i], extra)
            cell = centre(cell)          # a group label sits over its columns
            out.append(cell)
            i = j
        tbl = tbl.replace(hdr, rebuild(hdr, out), 1)
        label_row = rows[1] if len(rows) > 1 else None
    else:
        label_row = rows[0]

    if label_row:                    # rule under the labels, labels centred
        tbl = tbl.replace(label_row, rebuild(label_row,
              [centre(with_props(c, '<w:tcBorders>' + rule(LIGHT) + '</w:tcBorders>'))
               for c in cells(label_row)]), 1)

    cols = [int(x) for x in re.findall(r'<w:gridCol w:w="(\d+)"', tbl)]
    if cols and sum(cols):
        f = target / sum(cols)
        tbl = re.sub(r'(<w:gridCol w:w=")(\d+)(")',
                     lambda g: g.group(1) + str(round(int(g.group(2)) * f)) + g.group(3), tbl)
        tbl = re.sub(r'(<w:tcW w:w=")(\d+)(")',
                     lambda g: g.group(1) + str(round(int(g.group(2)) * f)) + g.group(3), tbl)
    tbl = re.sub(r'<w:tblW w:w="\d+" w:type="\w+"\s*/>',
                 '<w:tblW w:w="%d" w:type="dxa"/>' % target, tbl)
    return re.sub(r'(<w:tblPr>)', r'\1' + BORDERS + CELLMAR, tbl, count=1)

tmp, n = path + ".tmp", 0
zin  = zipfile.ZipFile(path)
zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
for it in zin.infolist():
    data = zin.read(it.filename)
    if it.filename == "word/document.xml":
        xml, parts, pos = data.decode("utf-8"), [], 0
        for m in re.finditer(r"<w:tbl>.*?</w:tbl>", xml, re.S):
            parts += [xml[pos:m.start()], style(m.group(0))]; pos = m.end(); n += 1
        parts.append(xml[pos:])
        data = "".join(parts).encode("utf-8")
    zout.writestr(it, data)
zout.close(); zin.close(); os.replace(tmp, path)
print("  styled %d tables at %.2f in" % (n, target / 1440))
