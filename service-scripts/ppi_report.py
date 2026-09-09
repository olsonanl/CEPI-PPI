#!/usr/bin/env python3
"""Build the single-page HTML summary for a CEPI-PPI job.

Reads the three files predict_ppi.py writes into the output directory plus the
two input FASTAs named in config.json, and emits a self-contained HTML report.
predict_ppi.py needs no change: everything here is derived from its output.

The sequence names come from the FASTAs rather than from the `ids` column of
the output, because that column is built as `<name1>_<name2>` and BV-BRC
feature ids contain underscores themselves (`511145.12.peg.648_511145.12.peg.246`
has no recoverable split point). Sequences are matched by their text instead.

Usage:
    ppi_report config.json                  # the job case
    ppi_report --output-dir D --query Q.fasta --target T.fasta --threshold 0.5
"""

import argparse
import ast
import csv
import datetime
import glob
import html
import json
import os
import sys
from collections import defaultdict

BLUE = "#196E9C"
LINK = "#0d78ef"

#
# Rendering caps. A job pairing every query against every target grows as the
# product of the two input counts, so each section that is O(pairs) or
# O(sequences) needs a ceiling or a 200x200 job writes a report no browser will
# open. The ranked table and the source TSV remain complete.
#
MAX_MATRIX_CELLS = 2500
MAX_TOP_PAIRS = 25
MAX_SEQUENCE_MAPS = 50

# csv fields hold a whole sequence and two per-residue vectors; the default
# 128 KiB limit is reached by a pair of ~1200-residue proteins.
csv.field_size_limit(min(2**31 - 1, sys.maxsize))


def read_fasta(path):
    """[(name, sequence)] in file order. Name is the first whitespace token."""
    name, buf, out = None, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                parts = line[1:].split()
                name, buf = (parts[0] if parts else ""), []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def name_index(seqs):
    """sequence text -> name, first occurrence wins.

    Returns (index, n_duplicates). Two records carrying the same sequence are
    indistinguishable downstream -- the model scored the sequence, not the
    name -- so the report says so rather than silently picking one.
    """
    idx, dups = {}, 0
    for name, seq in seqs:
        if seq in idx:
            dups += 1
        else:
            idx[seq] = name
    return idx, dups


def ramp(frac, full=0.35):
    """White to BV-BRC blue. `full` is the fraction that saturates the ramp."""
    if frac <= 0:
        return "#f7f7f7"
    a = min(1.0, frac / full) if full > 0 else 1.0
    r = 255 + (25 - 255) * a
    g = 255 + (110 - 255) * a
    b = 255 + (156 - 255) * a
    return "#%02x%02x%02x" % (int(r), int(g), int(b))


def runs(bits):
    """Lengths of the contiguous 1-runs in a 0/1 vector."""
    out, n = [], 0
    for v in list(bits) + [0]:
        if v:
            n += 1
        elif n:
            out.append(n)
            n = 0
    return out


def newest(outdir, suffix):
    hits = sorted(glob.glob(os.path.join(outdir, "*" + suffix)))
    return hits[-1] if hits else None


def load_pairs(class_labels, qidx, tidx):
    """Per-pair records from the class-labels file."""
    pairs, unresolved = [], 0
    with open(class_labels, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s1, s2 = row["Prot1"], row["Prot2"]
            try:
                c1 = ast.literal_eval(row["Classes1"])
                c2 = ast.literal_eval(row["Classes2"])
                l1 = ast.literal_eval(row["Labels1"])
                l2 = ast.literal_eval(row["Labels2"])
            except (ValueError, SyntaxError) as e:
                print("ppi_report: skipping unparseable row %s: %s" %
                      (row.get("ids", "?"), e), file=sys.stderr)
                continue
            if not c1 and not c2:
                continue
            q = qidx.get(s1)
            t = tidx.get(s2)
            if q is None or t is None:
                unresolved += 1
                ids = (row.get("ids") or "").strip()
                q = q or (ids or "query")
                t = t or (ids or "target")
            n1, n2 = sum(c1), sum(c2)
            pairs.append(dict(
                q=q, t=t, s1=s1, s2=s2, c1=c1, c2=c2, n1=n1, n2=n2,
                frac=(n1 + n2) / (len(c1) + len(c2)),
                mx=max(max(l1, default=0.0), max(l2, default=0.0)),
                patch=max(runs(c1) + runs(c2) + [0]),
            ))
    return pairs, unresolved


def consensus(pairs, seqs, role):
    """Per input sequence, how it behaves across all of its partners."""
    key = "c1" if role == "query" else "c2"
    skey = "s1" if role == "query" else "s2"
    out = {}
    for name, seq in seqs:
        mine = [p for p in pairs if p[skey] == seq]
        if not mine:
            continue
        hits = defaultdict(int)
        for p in mine:
            for i, v in enumerate(p[key]):
                if v:
                    hits[i] += 1
        out[name] = dict(
            seq=seq, role=role, partners=len(mine), hits=hits,
            constitutive=sum(1 for i in hits if hits[i] == len(mine)),
            union=len(hits),
            mean=sum(sum(p[key]) / len(p[key]) for p in mine) / len(mine),
        )
    return out


STYLE = """
body{font-family:'Work Sans',Helvetica,Arial,sans-serif;color:#3d4448;
 margin:2rem auto;max-width:1100px;padding:0 1.25rem;line-height:1.55}
h1{font-size:2em;font-weight:600;margin-bottom:.2rem}
h2{font-size:1.4em;font-weight:500;border-bottom:2px solid %(blue)s;
 padding-bottom:.2rem;margin-top:2.4rem}
h3{font-weight:500;margin-bottom:.3rem}
a{color:%(link)s;text-decoration:none} a:hover{text-decoration:underline}
table{border-collapse:collapse;margin:.7rem 0;font-size:.9em;
 font-variant-numeric:tabular-nums}
th{color:#fff;background:%(blue)s;padding:6px 10px;text-align:left;font-weight:500}
td{padding:5px 10px;border-bottom:1px solid #eceff1}
td.n{text-align:right}
.kv td:first-child{color:#6d7a82;padding-right:1.75rem;white-space:nowrap}
.note{background:#fbf7e8;border-left:4px solid #d8b13a;padding:.75rem 1rem;
 margin:1rem 0;font-size:.92em}
.wrap{overflow-x:auto;max-width:100%%}
.mat td{text-align:center;border:1px solid #fff;min-width:64px;line-height:1.25}
.mat th{font-size:.82em;white-space:nowrap}
.mat th.rh{background:#eef3f7;color:#333;text-align:right}
.seqmap{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
 font-size:12px;line-height:1.6;white-space:pre;overflow-x:auto;background:#fafbfc;
 padding:.7rem;border:1px solid #e6eaed;border-radius:2px}
.seqmap b{font-weight:600;border-radius:2px;padding:0 1px}
.bar{display:inline-block;height:9px;background:%(blue)s;vertical-align:middle;
 min-width:1px;border-radius:1px}
details{margin:.4rem 0} summary{cursor:pointer;color:%(link)s}
.small{font-size:.85em;color:#6d7a82}
""" % {"blue": BLUE, "link": LINK}


def render(pairs, queries, targets, qdups, tdups, unresolved,
           threshold, model, query_file, target_file, files):
    E = html.escape
    out = []
    A = out.append

    npairs = len(pairs)
    withhit = sum(1 for p in pairs if p["n1"] or p["n2"])
    res_tot = sum(len(p["c1"]) + len(p["c2"]) for p in pairs)
    res_hit = sum(p["n1"] + p["n2"] for p in pairs)
    allpatch = [x for p in pairs for x in runs(p["c1"]) + runs(p["c2"])]
    singletons = sum(1 for x in allpatch if x == 1)

    qcons = consensus(pairs, queries, "query")
    tcons = consensus(pairs, targets, "target")

    A("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
      "<title>Protein-Protein Interface Report</title>\n"
      "<link href=\"https://fonts.googleapis.com/css?family=Work+Sans:300,400,500,600,700\""
      " rel=\"stylesheet\">\n<style>%s</style>\n</head>\n<body>\n" % STYLE)

    A("<h1>Protein&#8211;Protein Interface Prediction</h1>\n"
      '<p class="small">Per-residue interface calls for every query &times; target pair, '
      "from a fine-tuned ESM2 token-classification model. Generated %s.</p>\n"
      % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

    # ---------------- Run ----------------
    A("<h2>Run</h2>\n<table class=\"kv\">\n")
    A("<tr><td>Query sequences</td><td><b>%d</b> &nbsp; <span class=\"small\">%s</span></td></tr>\n"
      % (len(queries), E(os.path.basename(query_file))))
    A("<tr><td>Target sequences</td><td><b>%d</b> &nbsp; <span class=\"small\">%s</span></td></tr>\n"
      % (len(targets), E(os.path.basename(target_file))))
    A("<tr><td>Pairs scored</td><td><b>%d</b></td></tr>\n" % npairs)
    A("<tr><td>Interface threshold</td><td>%s</td></tr>\n" % E(str(threshold)))
    A("<tr><td>Model</td><td class=\"small\">%s</td></tr>\n" % E(model))
    A("</table>\n")

    if not npairs:
        A('<div class="note"><b>No scored pairs were found in the job output.</b> '
          "The prediction step produced no rows; nothing can be summarized.</div>\n"
          "</body>\n</html>\n")
        return "".join(out)

    # ---------------- Summary ----------------
    A("<h2>Summary</h2>\n<table>\n"
      "<tr><th>Pairs with &ge;1 predicted interface residue</th><th>Residues scored</th>"
      "<th>Residues predicted interface</th><th>Overall rate</th></tr>\n")
    A('<tr><td class="n">%d / %d (%.0f%%)</td><td class="n">%s</td>'
      '<td class="n">%s</td><td class="n">%.1f%%</td></tr>\n'
      % (withhit, npairs, 100.0 * withhit / npairs, format(res_tot, ","),
         format(res_hit, ","), 100.0 * res_hit / res_tot if res_tot else 0.0))
    A("</table>\n")

    if withhit == 0:
        A('<div class="note"><b>No interface residues were predicted for any pair</b> '
          "at threshold %s. This is a normal outcome, not an error: it means no "
          "residue scored above the cutoff. The per-residue probabilities are in "
          "the results TSV if you want to look below it.</div>\n" % E(str(threshold)))

    if unresolved:
        A('<div class="note">%d scored pair(s) could not be matched back to a '
          "record in the input FASTAs and are labelled with their raw pair id.</div>\n"
          % unresolved)
    if qdups or tdups:
        A('<div class="note">%d input record(s) repeat a sequence that appears '
          "earlier in the same file. Identical sequences are scored identically, "
          "so only the first name of each is shown.</div>\n" % (qdups + tdups))

    # ---------------- Interaction matrix ----------------
    A("<h2>Interaction matrix</h2>\n")
    if len(queries) * len(targets) <= MAX_MATRIX_CELLS:
        A('<p class="small">Percent of residues called across both partners. '
          "Cell text is the raw count (query&nbsp;+&nbsp;target); hover for the "
          "highest per-residue probability in the pair.</p>\n")
        idx = {(p["q"], p["t"]): p for p in pairs}
        A('<div class="wrap"><table class="mat">\n<tr><th class="rh">query \\ target</th>')
        for name, _ in targets:
            A("<th>%s</th>" % E(name))
        A("</tr>\n")
        for qname, _ in queries:
            A('<tr><th class="rh">%s</th>' % E(qname))
            for tname, _ in targets:
                p = idx.get((qname, tname))
                if p is None:
                    A('<td style="background:#f7f7f7"></td>')
                    continue
                A('<td style="background:%s" title="max p = %.2f">%.0f%%<br>'
                  '<span class="small">%d+%d</span></td>'
                  % (ramp(p["frac"]), p["mx"], 100 * p["frac"], p["n1"], p["n2"]))
            A("</tr>\n")
        A("</table></div>\n")
    else:
        A('<p class="small">Matrix omitted: %d&times;%d exceeds the %s-cell '
          "rendering limit. The ranked table below and the results TSV cover "
          "the same data.</p>\n"
          % (len(queries), len(targets), format(MAX_MATRIX_CELLS, ",")))

    # ---------------- Top pairs ----------------
    A("<h2>Top pairs</h2>\n")
    A('<p class="small">Ranked by the fraction of residues called across both '
      "partners.</p>\n")
    A("<div class=\"wrap\"><table>\n<tr><th>Query</th><th>Target</th><th>Length</th>"
      "<th>Called</th><th>% called</th><th>Max p</th><th>Longest patch</th></tr>\n")
    for p in sorted(pairs, key=lambda p: -p["frac"])[:MAX_TOP_PAIRS]:
        A('<tr><td>%s</td><td>%s</td><td class="n">%d&nbsp;/&nbsp;%d</td>'
          '<td class="n">%d&nbsp;/&nbsp;%d</td><td class="n">%.1f%%</td>'
          '<td class="n">%.2f</td><td class="n">%d</td></tr>\n'
          % (E(p["q"]), E(p["t"]), len(p["c1"]), len(p["c2"]), p["n1"], p["n2"],
             100 * p["frac"], p["mx"], p["patch"]))
    A("</table></div>\n")
    if npairs > MAX_TOP_PAIRS:
        A('<p class="small">Showing %d of %d pairs. The results TSV has them all.</p>\n'
          % (MAX_TOP_PAIRS, npairs))

    # ---------------- Per-sequence ----------------
    A("<h2>Per-sequence summary</h2>\n")
    A('<p class="small">How each input sequence behaves across <i>all</i> of its '
      "partners. &ldquo;Constitutive&rdquo; residues are called against every "
      "partner; &ldquo;union&rdquo; residues are called against at least one. A "
      "sequence with many union but no constitutive residues is calling a "
      "different surface for each partner.</p>\n")
    for role, cons in (("Query", qcons), ("Target", tcons)):
        if not cons:
            continue
        A("<h3>%s sequences</h3>\n<div class=\"wrap\"><table>\n<tr><th>%s</th>"
          "<th>Length</th><th>Partners</th><th>Mean %% called</th>"
          "<th>Constitutive</th><th>Union</th><th></th></tr>\n" % (role, role))
        for name, c in cons.items():
            A('<tr><td>%s</td><td class="n">%d</td><td class="n">%d</td>'
              '<td class="n">%.1f%%</td><td class="n">%d</td><td class="n">%d</td>'
              '<td><span class="bar" style="width:%dpx"></span></td></tr>\n'
              % (E(name), len(c["seq"]), c["partners"], 100 * c["mean"],
                 c["constitutive"], c["union"],
                 int(min(1.0, c["mean"] / 0.35) * 160)))
        A("</table></div>\n")

    # ---------------- Consensus maps ----------------
    A("<h2>Consensus interface residues</h2>\n")
    mapped = [(n, c) for n, c in list(qcons.items()) + list(tcons.items()) if c["union"]]
    mapped.sort(key=lambda kv: -kv[1]["union"])
    if not mapped:
        A('<p class="small">No sequence had a residue called against any partner, '
          "so there is nothing to map.</p>\n")
    else:
        A('<p class="small">Each sequence with its called positions highlighted; '
          "darker means called against more partners. 60 residues per line, "
          "numbered from 1.</p>\n")
        shown = mapped[:MAX_SEQUENCE_MAPS]
        for name, c in shown:
            A("<details><summary>%s <b>%s</b> &mdash; %d residue(s) called across "
              "%d partner(s), %d constitutive</summary>\n<div class=\"seqmap\">"
              % (c["role"].capitalize(), E(name), c["union"], c["partners"],
                 c["constitutive"]))
            seq = c["seq"]
            for off in range(0, len(seq), 60):
                A('<span class="small">%5d </span>' % (off + 1))
                for i, ch in enumerate(seq[off:off + 60]):
                    n = c["hits"].get(off + i, 0)
                    if n:
                        A('<b style="background:%s">%s</b>'
                          % (ramp(0.35 * n / c["partners"]), E(ch)))
                    else:
                        A(E(ch))
                A('<span class="small"> %d</span>\n' % min(off + 60, len(seq)))
            A("</div></details>\n")
        if len(mapped) > MAX_SEQUENCE_MAPS:
            A('<p class="small">Showing the %d sequences with the most called '
              "residues, of %d with any.</p>\n" % (MAX_SEQUENCE_MAPS, len(mapped)))

    # ---------------- Caveats ----------------
    A("<h2>Interpreting these results</h2>\n"
      '<div class="note"><b>All-vs-all pairing is not evidence of interaction.</b> '
      "Every query is scored against every target whether or not the two proteins "
      "ever meet. A high score means &ldquo;these residues look like interface "
      "residues&rdquo;, not &ldquo;these two proteins bind&rdquo;.</div>\n<ul>\n")
    A("<li><b>Patch structure.</b> %s contiguous run(s) of called residues; "
      "%s (%.0f%%) are single isolated residues. Real interfaces are contiguous, "
      "so isolated calls are the least trustworthy part of the output.</li>\n"
      % (format(len(allpatch), ","), format(singletons, ","),
         100.0 * singletons / len(allpatch) if allpatch else 0.0))
    A("<li><b>The threshold is a dial, not a fact.</b> Everything above is at "
      "p &ge; %s. The per-residue probabilities are in the results TSV if you "
      "want to re-cut it.</li>\n" % E(str(threshold)))
    A("<li><b>The TP/TN/FP/FN columns in the results TSV are not meaningful "
      "here.</b> They compare against a ground-truth label vector this "
      "application fills with zeros, because prediction inputs have no known "
      "answer.</li>\n</ul>\n")

    # ---------------- Files ----------------
    A("<h2>Files</h2>\n<table>\n<tr><th>File</th><th>Contents</th></tr>\n")
    for fname, desc in files:
        A("<tr><td>%s</td><td>%s</td></tr>\n" % (E(fname), E(desc)))
    A("</table>\n</body>\n</html>\n")
    return "".join(out)


def parse_args(argv):
    p = argparse.ArgumentParser(
        prog="ppi_report",
        description="Build the HTML summary report for a CEPI-PPI job.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("config", nargs="?",
                   help="predict_ppi config.json; supplies the output directory, "
                        "the two input FASTAs, the threshold and the model")
    p.add_argument("--output-dir", help="directory holding the predict_ppi output "
                                        "(overrides the config)")
    p.add_argument("--query", help="query FASTA (overrides the config)")
    p.add_argument("--target", help="target FASTA (overrides the config)")
    p.add_argument("--threshold", type=float,
                   help="interface probability cutoff used for the run "
                        "(overrides the config)")
    p.add_argument("--model", help="model description for the Run table")
    p.add_argument("-o", "--outfile",
                   help="report path (default: <output-dir>/<prefix>_interface_report.html)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="do not print the report path on success")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])

    outdir = args.output_dir
    query_file, target_file = args.query, args.target
    threshold, model = args.threshold, args.model

    if args.config:
        with open(args.config) as fh:
            cfg = json.load(fh)
        params = cfg.get("params") or {}
        outdir = outdir or cfg.get("output_data_dir")
        query_file = query_file or cfg.get("query")
        target_file = target_file or cfg.get("target")
        if threshold is None:
            threshold = cfg.get("threshold", params.get("threshold"))
        if model is None:
            model = " + ".join(x for x in (params.get("pt_model"),
                                           os.path.basename(cfg.get("model_path") or ""))
                               if x)

    missing = [n for n, v in (("output directory", outdir), ("query FASTA", query_file),
                              ("target FASTA", target_file)) if not v]
    if missing:
        print("ppi_report: no %s given (pass config.json or the matching option)"
              % ", ".join(missing), file=sys.stderr)
        return 2

    class_labels = newest(outdir, "_class_labels.txt")
    inference = newest(outdir, "_inference_results.tsv")
    if class_labels is None:
        print("ppi_report: no *_class_labels.txt in %s; nothing to report on" % outdir,
              file=sys.stderr)
        return 1

    threshold = threshold if threshold is not None else "unspecified"
    model = model or "unspecified"

    queries = read_fasta(query_file)
    targets = read_fasta(target_file)
    qidx, qdups = name_index(queries)
    tidx, tdups = name_index(targets)

    pairs, unresolved = load_pairs(class_labels, qidx, tidx)

    files = [(os.path.basename(class_labels),
              "Per pair: per-residue interface call (0/1) at the threshold")]
    if inference:
        files.append((os.path.basename(inference),
                      "Per pair: both sequences, per-residue probability, counts"))
    ids_file = newest(outdir, "_sequence_ids.txt")
    if ids_file:
        files.append((os.path.basename(ids_file),
                      "Pair id and the two input sequences"))

    doc = render(pairs, queries, targets, qdups, tdups, unresolved,
                 threshold, model, query_file, target_file, files)

    dest = args.outfile
    if not dest:
        base = os.path.basename(class_labels)
        prefix = base[:-len("_class_labels.txt")]
        dest = os.path.join(outdir, (prefix + "_" if prefix else "") +
                            "interface_report.html")
    with open(dest, "w") as fh:
        fh.write(doc)
    if not args.quiet:
        print("ppi_report: wrote %s" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
