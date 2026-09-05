#!/bin/bash
# make_power_table_html.sh -- the power table as HTML, for pasting into Google Docs.
#
#   bash temp_aging/make_power_table_html.sh
#
# Open temp_aging/power_table.html in a browser, select all, copy, paste into the
# document. Docs keeps table structure, colspan, rules and alignment from an HTML
# paste and turns it into a native Docs table, so nothing HTML reaches the
# journal -- this is transport only.
set -euo pipefail
cd "$(dirname "$0")/.."
cat > temp_aging/power_table.html <<'HTML'
<!doctype html><meta charset="utf-8">
<body style="font-family:Georgia,'Times New Roman',serif;font-size:11pt;margin:40px">

<table style="border-collapse:collapse;width:100%;border-top:2px solid #000;border-bottom:2px solid #000">
  <tr>
    <td style="padding:3px 10px"></td>
    <td style="padding:3px 10px"></td>
    <td colspan="3" style="padding:3px 10px;text-align:center;border-bottom:1px solid #000">h&sup2;</td>
  </tr>
  <tr>
    <td style="padding:3px 10px;border-bottom:1.5px solid #000"></td>
    <td style="padding:3px 10px;text-align:center;border-bottom:1.5px solid #000">n</td>
    <td style="padding:3px 10px;text-align:center;border-bottom:1.5px solid #000">2%</td>
    <td style="padding:3px 10px;text-align:center;border-bottom:1.5px solid #000">1%</td>
    <td style="padding:3px 10px;text-align:center;border-bottom:1.5px solid #000">0.5%</td>
  </tr>
  <tr><td style="padding:3px 10px">Outbred panel</td>
      <td style="padding:3px 10px;text-align:right">1,000</td>
      <td style="padding:3px 10px;text-align:right">0.35</td>
      <td style="padding:3px 10px;text-align:right">0.04</td>
      <td style="padding:3px 10px;text-align:right">0.004</td></tr>
  <tr><td style="padding:3px 10px">Inbred lines</td>
      <td style="padding:3px 10px;text-align:right">200</td>
      <td style="padding:3px 10px;text-align:right">0.02</td>
      <td style="padding:3px 10px;text-align:right">0.00</td>
      <td style="padding:3px 10px;text-align:right">0.000</td></tr>
  <tr><td style="padding:3px 10px">Eight-way RILs</td>
      <td style="padding:3px 10px;text-align:right">750</td>
      <td style="padding:3px 10px;text-align:right">0.70</td>
      <td style="padding:3px 10px;text-align:right">0.16</td>
      <td style="padding:3px 10px;text-align:right">0.022</td></tr>
  <tr><td style="padding:3px 10px">X-QTL</td>
      <td style="padding:3px 10px;text-align:right">114,000</td>
      <td style="padding:3px 10px;text-align:right">1.00</td>
      <td style="padding:3px 10px;text-align:right">1.00</td>
      <td style="padding:3px 10px;text-align:right">0.894</td></tr>
</table>

<p style="margin-top:8px">Note. &mdash; Power to detect a locus of heritability h&sup2; at the size each
design is run at. Thresholds are &minus;log<sub>10</sub>P of 6 for the outbred and
inbred-line panels, 4 for the RIL scan and 7.5 for X-QTL.</p>

</body>
HTML
echo "Wrote temp_aging/power_table.html"
