#!/usr/bin/env python3
"""Generates docs/stories.html and docs/stories/*.csv: six end-to-end stories of increasing complexity, each a small bakery question answered by REAL runs of build/stats. Nothing on the page is typed by hand except the
narration: every number, every excerpt and every figure is drawn from the JSON the tool printed on the data of docs/stories/. The data are invented (a seeded generator below); the runs are real.

    python3 scripts/gen_stories.py [--bin build/stats] [--check]

`--check` regenerates in memory and fails if docs/stories.html or the CSV files differ: the tool's answers are the same bytes everywhere, so a page that is stale is a page that is wrong."""
import argparse
import hashlib
import html
import json
import pathlib
import random
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


# ------------------------------------------------------------------------------------------------ the invented data
def make_bakery():
    """80 batches, 40 per oven. Oven B bakes faster and rises higher, but three of its batches ran long; the crust and the weight do not differ; the salt does not differ; B has fewer defects."""
    r = random.Random(5)
    cols = []
    for oven, mu, sd in (("A", 43.0, 3.5), ("B", 39.0, 3.0)):
        for i in range(40):
            t = r.gauss(mu, sd)
            if oven == "B" and i in (7, 23, 31):
                t += 15
            cols.append({"oven": oven, "bake_min": round(t, 1)})
    r2 = random.Random(23)
    for row in cols:
        b = row["oven"] == "B"
        row["height_cm"] = round(r2.gauss(11.9 if b else 11.0, 0.8), 1)
        row["crust"] = round(min(10, max(1, r2.gauss(6.4 if b else 6.2, 1.3))), 1)
        row["weight_g"] = round(r2.gauss(480, 12), 0)
        row["moisture_pct"] = round(r2.gauss(36.1 if b else 37.2, 1.4), 1)
        row["salt_pct"] = round(r2.gauss(1.8, 0.1), 2)
        row["defect"] = "y" if r2.random() < (0.025 if b else 0.3) else "n"
    r.shuffle(cols)
    return cols


def make_bio():
    """log2 expression of one gene in the liver of six wild-type and six knockout mice; one knockout animal is high"""
    r = random.Random(4)
    rows = [("wildtype", round(r.gauss(8.0, 0.45), 2)) for _ in range(6)] + [("knockout", round(r.gauss(9.0, 0.9), 2)) for _ in range(6)]
    rows[-1] = ("knockout", round(rows[-1][1] + 1.8, 2))
    return rows


def make_chem():
    """caffeine (mg per serving) in twelve samples by HPLC and by UV absorbance; the UV method reads 2% high"""
    r = random.Random(3)
    rows = []
    for i in range(12):
        t = round(r.uniform(20, 80), 1)
        rows.append((i + 1, round(t + r.gauss(0, 0.5), 2), round(t * 1.02 + r.gauss(0, 0.5), 2)))
    return rows


def make_physics():
    """photopeak energies (keV) of the same source on two detectors: one clean, one whose gain drifts so that the peak splits in two (same mean)"""
    r = random.Random(5)
    rows = [("det1", round(r.gauss(662, 16), 2)) for _ in range(30)]
    rows += [("det2", round(r.gauss(662 + (22 if r.random() < 0.5 else -22), 6), 2)) for _ in range(30)]
    return rows


def make_clock():
    """twenty timestamps (seconds since the epoch) of one repeated event, microseconds of jitter on a large offset"""
    r = random.Random(8)
    return [round(1700000000.5 + r.gauss(0, 3e-6), 7) for _ in range(20)]


def make_town():
    """150 residents of three districts: support for a policy, satisfaction on a 1-5 scale (higher further south), education"""
    r = random.Random(1)
    rows = []
    for d, pm, lik in (("North", 0.45, [0.15, 0.25, 0.3, 0.2, 0.1]), ("Centre", 0.6, [0.08, 0.17, 0.3, 0.28, 0.17]), ("South", 0.72, [0.04, 0.1, 0.26, 0.32, 0.28])):
        for _ in range(50):
            x = r.random()
            c = 0
            level = 5
            for k, pk in enumerate(lik):
                c += pk
                if x <= c:
                    level = k + 1
                    break
            rows.append((d, "yes" if r.random() < pm else "no", level, r.choice(["school", "college", "degree"])))
    return rows


COLS = ["oven", "bake_min", "height_cm", "crust", "weight_g", "moisture_pct", "salt_pct", "defect"]


def csv_text(rows, cols):
    def fmt(v):
        return ("%d" % v) if isinstance(v, float) and v == int(v) and abs(v) >= 100 else str(v)
    return ",".join(cols) + "\n" + "\n".join(",".join(fmt(row[c]) for c in cols) for row in rows) + "\n"


# ------------------------------------------------------------------------------------------------ real runs
class Runner:
    def __init__(self, binary, folder):
        self.binary = binary
        self.folder = folder

    def run(self, args, file="bakery.csv"):
        p = subprocess.run([self.binary] + args + ["--root", str(self.folder), file], capture_output=True, text=True, timeout=600)
        if p.returncode != 0:
            raise SystemExit("gen_stories: `%s` exited %d: %s" % (" ".join(args), p.returncode, p.stdout[:300]))
        return json.loads(p.stdout)["data"], p.stdout


def refusal(binary, folder, args, file):
    """the real refusal: the document the tool prints with a nonzero exit"""
    p = subprocess.run([binary] + args + ["--root", str(folder), file], capture_output=True, text=True, timeout=600)
    if p.returncode == 0:
        raise SystemExit("gen_stories: `%s` was expected to be refused" % " ".join(args))
    return p.returncode, json.loads(p.stdout)["error"]


def claim(cond, what):
    """the narration says something about the numbers: if the numbers stop saying it, the page is not generated"""
    if not cond:
        raise SystemExit("gen_stories: the narration no longer matches the runs: " + what)


def shown(args, file="bakery.csv"):
    return "stats " + " ".join(a if re.fullmatch(r"[\w.,:=/-]+", a) else "'%s'" % a for a in args) + " " + file


def pretty(obj, width=78, indent=0):
    """indented JSON in which every object or list that fits on one line is kept on one line"""
    flat = json.dumps(obj, ensure_ascii=False)
    pad = "  " * indent
    if len(flat) + len(pad) <= width or not isinstance(obj, (dict, list)) or not obj:
        return flat
    if isinstance(obj, dict):
        items = ["%s  %s: %s" % (pad, json.dumps(k, ensure_ascii=False), pretty(v, width, indent + 1)) for k, v in obj.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    items = ["%s  %s" % (pad, pretty(v, width, indent + 1)) for v in obj]
    return "[\n" + ",\n".join(items) + "\n" + pad + "]"


def num(x, k=4):
    s = ("%." + str(k) + "g") % x
    return s


# ------------------------------------------------------------------------------------------------ figures (inline SVG, colours from the page's variables)
W = 640


def esc(s):
    return html.escape(str(s), quote=True)


def axis(x0, x1, y, lo, hi, ticks, label, fmt=lambda v: "%g" % v):
    sx = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
    out = ['<line class="ax" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (x0, y, x1, y)]
    for t in ticks:
        out.append('<line class="ax" x1="%g" y1="%g" x2="%g" y2="%g"/><text class="tk" x="%g" y="%g" text-anchor="middle">%s</text>' % (sx(t), y, sx(t), y + 5, sx(t), y + 19, esc(fmt(t))))
    out.append('<text class="lb" x="%g" y="%g" text-anchor="middle">%s</text>' % ((x0 + x1) / 2, y + 38, esc(label)))
    return "\n".join(out), sx


def svg(h, body, label):
    return '<figure class="fig"><svg viewBox="0 0 %d %d" role="img" aria-label="%s">%s</svg></figure>' % (W, h, esc(label), body)


def fig_boxes(groups):
    """horizontal box plots of bake_min per oven from the quantiles describe printed: min .. max whiskers, q1-q3 box, median bar, mean diamond"""
    lo, hi = 30, 62
    ax, sx = axis(60, 610, 150, lo, hi, [30, 35, 40, 45, 50, 55, 60], "minutes in the oven")
    body = [ax]
    for k, (name, c, cls) in enumerate(groups):
        y = 40 + k * 55
        q = {x["p"]: x["value"] for x in c["quantiles"]}
        body.append('<text class="nm" x="52" y="%g" text-anchor="end">oven %s</text>' % (y + 5, esc(name)))
        body.append('<line class="wh %s" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (cls, sx(q["0"]), y, sx(q["1"]), y))
        body.append('<line class="wh %s" x1="%g" y1="%g" x2="%g" y2="%g"/><line class="wh %s" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (cls, sx(q["0"]), y - 7, sx(q["0"]), y + 7, cls, sx(q["1"]), y - 7, sx(q["1"]), y + 7))
        body.append('<rect class="bx %s" x="%g" y="%g" width="%g" height="22" rx="3"/>' % (cls, sx(q["0.25"]), y - 11, sx(q["0.75"]) - sx(q["0.25"])))
        body.append('<line class="md" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (sx(q["0.5"]), y - 11, sx(q["0.5"]), y + 11))
        mx = sx(c["mean"])
        body.append('<path class="mn" d="M%g %g l5 5 l-5 5 l-5 -5 z"/>' % (mx, y - 21 + 8))
        body.append('<text class="tk" x="%g" y="%g" text-anchor="middle">max %s</text>' % (sx(q["1"]), y - 14, esc("%g" % q["1"])))
    body.append('<text class="tk" x="610" y="14" text-anchor="end">box: quartiles &#183; bar: median &#183; diamond: mean &#183; whiskers: min and max</text>')
    return svg(200, "\n".join(body), "Box plots of the baking time of each oven")


def fig_interval(est, lo_, hi_, p, xlo, xhi, xlabel):
    ticks = [t for t in range(int(xlo), int(xhi) + 1)]
    ax, sx = axis(60, 610, 110, xlo, xhi, ticks, xlabel)
    body = [ax]
    body.append('<line class="zero" x1="%g" y1="30" x2="%g" y2="110"/><text class="tk" x="%g" y="24" text-anchor="middle">no difference</text>' % (sx(0), sx(0), sx(0)))
    body.append('<line class="ciline" x1="%g" y1="70" x2="%g" y2="70"/><line class="ciline" x1="%g" y1="62" x2="%g" y2="78"/><line class="ciline" x1="%g" y1="62" x2="%g" y2="78"/>' % (sx(lo_), sx(hi_), sx(lo_), sx(lo_), sx(hi_), sx(hi_)))
    body.append('<circle class="pt" cx="%g" cy="70" r="6"/>' % sx(est))
    body.append('<text class="tk" x="%g" y="56" text-anchor="middle">%s</text>' % (sx(est), esc("%.2f" % est)))
    body.append('<text class="tk" x="%g" y="94" text-anchor="end">%s</text><text class="tk" x="%g" y="94" text-anchor="middle">%s</text>' % (sx(lo_) - 4, esc("%.2f" % lo_), sx(hi_), esc("%.2f" % hi_)))
    return svg(160, "\n".join(body), "The difference of mean baking times with its 95% interval and the line of no difference")


def fig_strip(values):
    """one dot per batch, by oven, with the long-running batches ringed"""
    lo, hi = 30, 62
    ax, sx = axis(60, 610, 170, lo, hi, [30, 35, 40, 45, 50, 55, 60], "minutes in the oven")
    body = [ax]
    for k, (name, vs, cls) in enumerate(values):
        y0 = 55 + k * 60
        body.append('<text class="nm" x="52" y="%g" text-anchor="end">oven %s</text>' % (y0 + 5, esc(name)))
        body.append('<line class="band" x1="60" y1="%g" x2="610" y2="%g"/>' % (y0, y0))
        for i, v in enumerate(sorted(vs)):
            jitter = ((i * 7919) % 11 - 5) * 2.0
            body.append('<circle class="dot %s" cx="%g" cy="%g" r="3.6"/>' % (cls, sx(v), y0 + jitter))
            if v > 53:
                body.append('<circle class="ring" cx="%g" cy="%g" r="8"/>' % (sx(v), y0 + jitter))
    body.append('<text class="tk" x="610" y="14" text-anchor="end">ringed: the three batches that ran past 53 minutes</text>')
    return svg(215, "\n".join(body), "Every batch as a dot, by oven, with the batches that ran long ringed")


def fig_dots(values, lo, hi, ticks, label):
    """one dot per animal or sample, by group"""
    ax, sx = axis(120, 610, 40 + 60 * len(values), lo, hi, ticks, label)
    body = [ax]
    for k, (name, vs, cls) in enumerate(values):
        y0 = 40 + k * 60
        body.append('<text class="nm" x="112" y="%g" text-anchor="end">%s</text><line class="band" x1="120" y1="%g" x2="610" y2="%g"/>' % (y0 + 5, esc(name), y0, y0))
        for i, v in enumerate(sorted(vs)):
            body.append('<circle class="dot %s" cx="%g" cy="%g" r="5"/>' % (cls, sx(v), y0 + ((i * 7919) % 5 - 2) * 3.0))
    return svg(40 + 60 * len(values) + 40, "\n".join(body), "Every animal as a dot, by genotype")


def fig_pbars(items):
    """p-values on a log axis, 0.05 marked"""
    import math
    lo, hi = math.log10(0.0005), math.log10(1.0)
    y = 40 + 34 * len(items)
    ax, sx = axis(210, 520, y, lo, hi, [math.log10(v) for v in (0.001, 0.01, 0.1, 1.0)], "p-value (log scale)", fmt=lambda v: "%g" % (10 ** v))
    body = [ax]
    body.append('<line class="zero" x1="%g" y1="22" x2="%g" y2="%g"/><text class="tk" x="%g" y="16" text-anchor="middle">0.05</text>' % (sx(math.log10(0.05)), sx(math.log10(0.05)), y, sx(math.log10(0.05))))
    for k, (name, p, cls) in enumerate(items):
        yy = 40 + 34 * k
        body.append('<text class="nm" x="200" y="%g" text-anchor="end">%s</text>' % (yy + 5, esc(name)))
        w = max(2, sx(max(p, 0.0005)) - 210)
        body.append('<rect class="pbar %s" x="210" y="%g" width="%g" height="16" rx="3"/>' % (cls, yy - 8, w))
        body.append('<text class="tk" x="%g" y="%g">%s</text>' % (210 + w + 6, yy + 5, esc("%.4f" % p if p >= 0.001 else "%.1e" % p)))
    return svg(y + 52, "\n".join(body), "The p-values of the tests on a logarithmic axis against the line at 0.05")


def fig_adjust(rows):
    """per measure: raw p (hollow) and adjusted p (filled) on a log axis, 0.05 marked"""
    import math
    lo, hi = math.log10(0.0005), math.log10(1.0)
    y = 40 + 34 * len(rows)
    ax, sx = axis(210, 600, y, lo, hi, [math.log10(v) for v in (0.001, 0.01, 0.1, 1.0)], "p-value (log scale)", fmt=lambda v: "%g" % (10 ** v))
    body = [ax]
    body.append('<line class="zero" x1="%g" y1="22" x2="%g" y2="%g"/><text class="tk" x="%g" y="16" text-anchor="middle">0.05</text>' % (sx(math.log10(0.05)), sx(math.log10(0.05)), y, sx(math.log10(0.05))))
    for k, (name, raw, adj) in enumerate(rows):
        yy = 40 + 34 * k
        body.append('<text class="nm" x="200" y="%g" text-anchor="end">%s</text>' % (yy + 5, esc(name)))
        xr, xa = sx(math.log10(max(raw, 0.0005))), sx(math.log10(max(adj, 0.0005)))
        body.append('<line class="link" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (xr, yy, xa, yy))
        body.append('<circle class="raw" cx="%g" cy="%g" r="6"/><circle class="adj" cx="%g" cy="%g" r="6"/>' % (xr, yy, xa, yy))
    body.append('<text class="tk" x="600" y="%g" text-anchor="end">&#9675; p as computed &#183; &#9679; p after Holm&#8217;s correction</text>' % (y + 52))
    return svg(y + 62, "\n".join(body), "For each measure the p-value before and after the correction for six tests")


def fig_arrangements(total, extreme, cols=14):
    cells = []
    rows = (total + cols - 1) // cols
    for i in range(total):
        x, y = 60 + (i % cols) * 24, 14 + (i // cols) * 24
        cells.append('<rect class="cell %s" x="%g" y="%g" width="19" height="19" rx="3"/>' % ("hit" if i < extreme else "", x, y))
    cells.append('<text class="nm" x="60" y="%g">%d of the %d ways to split the eight batches (red) are as extreme as the real split</text>' % (14 + rows * 24 + 14, extreme, total))
    return svg(14 + rows * 24 + 28, "\n".join(cells), "%d of the %d ways to split the pilot batches between the ovens are as extreme as the real one" % (extreme, total))


def fig_two_by_two(a_yes, a_n, b_yes, b_n):
    body = []
    for k, (name, yes, n, cls) in enumerate((("oven A", a_yes, a_n, "pa"), ("oven B", b_yes, b_n, "pb"))):
        y = 30 + k * 50
        body.append('<text class="nm" x="90" y="%g" text-anchor="end">%s</text>' % (y + 6, name))
        body.append('<rect class="seg" x="100" y="%g" width="%g" height="22" rx="3"/>' % (y - 8, 360 * (n - yes) / max(a_n, b_n)))
        body.append('<rect class="seg bad" x="%g" y="%g" width="%g" height="22" rx="3"/>' % (100 + 360 * (n - yes) / max(a_n, b_n), y - 8, 360 * yes / max(a_n, b_n)))
        body.append('<text class="tk" x="%g" y="%g">%d of %d defective</text>' % (100 + 360 * n / max(a_n, b_n) + 10, y + 7, yes, n))
    return svg(130, "\n".join(body), "Defective loaves out of the batches of each oven")


def fig_scatter(pairs):
    """one dot per sample, method 1 against method 2, with the line of agreement"""
    lo, hi = 15, 85
    ax, sx = axis(60, 340, 290, lo, hi, [20, 40, 60, 80], "HPLC, mg per serving")
    sy = lambda v: 270 - (v - lo) / (hi - lo) * 250
    body = [ax, '<line class="ax" x1="60" y1="20" x2="60" y2="290"/>']
    for t in (20, 40, 60, 80):
        body.append('<text class="tk" x="52" y="%g" text-anchor="end">%d</text>' % (sy(t) + 4, t))
    body.append('<text class="lb" x="14" y="150" transform="rotate(-90 14 150)" text-anchor="middle">UV, mg per serving</text>')
    body.append('<line class="zero" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (sx(lo), sy(lo), sx(hi), sy(hi)))
    for h, u in pairs:
        body.append('<circle class="dot pa" cx="%g" cy="%g" r="4"/>' % (sx(h), sy(u)))
    body.append('<text class="tk" x="340" y="14" text-anchor="end">dashed: agreement</text>')
    return body, sx


def fig_chem(pairs, est, lo_, hi_):
    left, _ = fig_scatter(pairs)
    # the differences, against the mean of the two methods, with the interval of the mean difference as a band
    mlo, mhi = 15, 85
    dlo, dhi = -2.5, 1.0
    sxm = lambda v: 400 + (v - mlo) / (mhi - mlo) * 220
    syd = lambda v: 270 - (v - dlo) / (dhi - dlo) * 250
    right = ['<line class="ax" x1="400" y1="290" x2="620" y2="290"/><line class="ax" x1="400" y1="20" x2="400" y2="290"/>']
    for t in (20, 40, 60, 80):
        right.append('<line class="ax" x1="%g" y1="290" x2="%g" y2="295"/><text class="tk" x="%g" y="309" text-anchor="middle">%d</text>' % (sxm(t), sxm(t), sxm(t), t))
    for t in (-2, -1, 0, 1):
        right.append('<text class="tk" x="392" y="%g" text-anchor="end">%d</text>' % (syd(t) + 4, t))
    right.append('<text class="lb" x="510" y="328" text-anchor="middle">mean of the two methods, mg</text>')
    right.append('<text class="lb" x="352" y="150" transform="rotate(-90 352 150)" text-anchor="middle">HPLC minus UV, mg</text>')
    right.append('<rect class="bandci" x="400" y="%g" width="220" height="%g"/>' % (syd(hi_), syd(lo_) - syd(hi_)))
    right.append('<line class="zero" x1="400" y1="%g" x2="620" y2="%g"/>' % (syd(0), syd(0)))
    right.append('<line class="ciline" x1="400" y1="%g" x2="620" y2="%g"/>' % (syd(est), syd(est)))
    for h, u in pairs:
        right.append('<circle class="dot pb" cx="%g" cy="%g" r="4"/>' % (sxm((h + u) / 2), syd(h - u)))
    right.append('<text class="tk" x="620" y="345" text-anchor="end">band: 95%% interval of the mean difference, %.2f to %.2f</text>' % (lo_, hi_))
    return svg(358, "\n".join(left + right), "Left: the two methods against each other. Right: their differences with the interval of the mean difference")


def fig_ecdf(a, b, d):
    """the two empirical distribution functions and the largest gap between them"""
    lo, hi = 610, 720
    ax, sx = axis(60, 610, 230, lo, hi, [620, 640, 660, 680, 700, 720], "photopeak energy, keV")
    sy = lambda v: 210 - v * 190
    body = [ax, '<line class="ax" x1="60" y1="20" x2="60" y2="230"/>']
    for t in (0, 0.5, 1):
        body.append('<text class="tk" x="52" y="%g" text-anchor="end">%g</text>' % (sy(t) + 4, t))
    best = (0, None)
    pts = sorted(a + b)
    def F(v, xs): return sum(1 for x in xs if x <= v) / len(xs)
    for name, xs, cls in (("det1", a, "ecdfa"), ("det2", b, "ecdfb")):
        xs = sorted(xs)
        path = "M%g %g" % (sx(lo), sy(0))
        for i, v in enumerate(xs):
            path += " L%g %g L%g %g" % (sx(v), sy(i / len(xs)), sx(v), sy((i + 1) / len(xs)))
        path += " L%g %g" % (sx(hi), sy(1))
        body.append('<path class="%s" d="%s"/>' % (cls, path))
    for v in pts:
        g = abs(F(v, a) - F(v, b))
        if g > best[0]:
            best = (g, v)
    body.append('<line class="gap" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (sx(best[1]), sy(F(best[1], a)), sx(best[1]), sy(F(best[1], b))))
    body.append('<text class="tk" x="%g" y="%g">largest gap D = %.3f</text>' % (sx(best[1]) + 8, (sy(F(best[1], a)) + sy(F(best[1], b))) / 2, best[0]))
    body.append('<text class="nm" x="76" y="30" style="fill:var(--accent)">detector 1</text><text class="nm" x="76" y="48" style="fill:var(--good)">detector 2</text>')
    return svg(275, "\n".join(body), "The empirical distribution functions of the two detectors and their largest gap"), best[0]


def fig_prop(items):
    """intervals of a proportion, with the line at one half"""
    ax, sx = axis(60, 610, 100, 0.3, 0.8, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8], "share of residents who support the policy")
    body = [ax, '<line class="zero" x1="%g" y1="26" x2="%g" y2="100"/><text class="tk" x="%g" y="20" text-anchor="middle">one half</text>' % (sx(0.5), sx(0.5), sx(0.5))]
    for k, (name, est, lo_, hi_) in enumerate(items):
        y = 46 + k * 30
        body.append('<text class="nm" x="170" y="%g" text-anchor="end">%s</text>' % (y + 4, esc(name)))
        body.append('<line class="ciline" x1="%g" y1="%g" x2="%g" y2="%g"/><circle class="pt" cx="%g" cy="%g" r="5"/>' % (sx(lo_), y, sx(hi_), y, sx(est), y))
        body.append('<text class="tk" x="%g" y="%g">%.3f &#8211; %.3f</text>' % (sx(hi_) + 8, y + 4, lo_, hi_))
    return svg(150, "\n".join(body), "Two 95% intervals of the share of supporters")


def fig_likert(groups):
    """stacked bars of the five answers, by district"""
    shades = ["lk1", "lk2", "lk3", "lk4", "lk5"]
    body = []
    for k, (name, counts) in enumerate(groups):
        y = 36 + k * 44
        tot = sum(counts)
        x = 110
        body.append('<text class="nm" x="102" y="%g" text-anchor="end">%s</text>' % (y + 5, esc(name)))
        for lv, c in enumerate(counts):
            w = 480 * c / tot
            body.append('<rect class="%s" x="%g" y="%g" width="%g" height="24"/>' % (shades[lv], x, y - 8, w))
            if c >= 4:
                body.append('<text class="lk" x="%g" y="%g" text-anchor="middle">%d</text>' % (x + w / 2, y + 8, c))
            x += w
    body.append('<text class="tk" x="110" y="16">1 very dissatisfied</text><text class="tk" x="590" y="16" text-anchor="end">5 very satisfied</text>')
    return svg(36 + 44 * len(groups) + 8, "\n".join(body), "How the residents of each district answered the satisfaction question, 1 to 5")


def fig_hashes(rows):
    body = []
    for k, (label, h) in enumerate(rows):
        y = 28 + k * 34
        body.append('<text class="nm" x="205" y="%g" text-anchor="end">%s</text><rect class="hash" x="215" y="%g" width="395" height="22" rx="4"/><text class="hx" x="225" y="%g">%s</text>' % (y + 5, esc(label), y - 10, y + 5, esc(h)))
    return svg(28 + 34 * len(rows) + 10, "\n".join(body), "The digests of the answers of the same question asked with the rows in different orders")


# ------------------------------------------------------------------------------------------------ the page
STORY_CSS = """
.hero.small { padding: 3rem 0 1rem; }
.story { padding: 3rem 0; }
.story .wrap { max-width: 50rem; }
.story h2 { text-align: left; margin-bottom: .3rem; }
.badge { display: inline-flex; gap: .4rem; flex-wrap: wrap; margin: 0 0 1rem; }
.badge span { font: .78rem var(--mono); color: var(--muted); border: 1px solid var(--line); border-radius: 99px; padding: .1rem .6rem; background: var(--surface); }
.badge .lvl { color: var(--accent); border-color: var(--accent); }
.said { border-left: 3px solid var(--accent); padding: .2rem 0 .2rem 1rem; margin: 1rem 0; color: var(--muted); font-style: italic; }
.story .term { max-width: none; margin: 1rem 0; box-shadow: none; }
.story .term pre { min-height: 0; white-space: pre; word-break: normal; font-size: .8rem; }
.story pre .cm { color: var(--term-acc); }
.read { background: var(--surface); border: 1px solid var(--line); border-radius: 12px; padding: 1rem 1.2rem; margin: 1rem 0; }
.read b { color: var(--accent); }
figure.fig { margin: 1rem 0; }
figure.fig svg { width: 100%; height: auto; display: block; background: var(--surface); border: 1px solid var(--line); border-radius: 12px; }
.fig .ax { stroke: var(--muted); stroke-width: 1; }
.fig .tk { fill: var(--muted); font: 12px var(--font); }
.fig .lb { fill: var(--muted); font: 13px var(--font); }
.fig .nm { fill: var(--fg); font: 600 13px var(--font); }
.fig .wh { stroke-width: 2; } .fig .wh.pa, .fig .bx.pa { stroke: var(--accent); } .fig .wh.pb, .fig .bx.pb { stroke: var(--good); }
.fig .bx { fill: var(--soft); stroke-width: 2; }
.fig .md { stroke: var(--fg); stroke-width: 3; }
.fig .mn { fill: var(--warn); }
.fig .zero { stroke: var(--bad); stroke-width: 1.5; stroke-dasharray: 5 4; }
.fig .ciline { stroke: var(--accent); stroke-width: 3; }
.fig .pt { fill: var(--accent); }
.fig .band { stroke: var(--line); stroke-width: 1; }
.fig .dot.pa { fill: var(--accent); opacity: .8; } .fig .dot.pb { fill: var(--good); opacity: .8; }
.fig .ring { fill: none; stroke: var(--bad); stroke-width: 2; }
.fig .pbar.pa { fill: var(--accent); } .fig .pbar.pb { fill: var(--good); } .fig .pbar.warn { fill: var(--warn); }
.fig .link { stroke: var(--muted); stroke-width: 1.5; }
.fig .raw { fill: var(--surface); stroke: var(--accent); stroke-width: 2.5; } .fig .adj { fill: var(--accent); }
.fig .cell { fill: var(--soft); stroke: var(--line); } .fig .cell.hit { fill: var(--bad); stroke: var(--bad); }
.fig .seg { fill: var(--good); } .fig .seg.bad { fill: var(--bad); }
.fig .hash { fill: var(--soft); stroke: var(--good); stroke-width: 1.5; }
.fig .hx { fill: var(--fg); font: 12px var(--mono); }
.fig .bandci { fill: var(--accent); opacity: .18; }
.fig .ecdfa { fill: none; stroke: var(--accent); stroke-width: 2.5; } .fig .ecdfb { fill: none; stroke: var(--good); stroke-width: 2.5; }
.fig .gap { stroke: var(--bad); stroke-width: 3; }
.fig .lk1 { fill: #b42d27; } .fig .lk2 { fill: #d97a52; } .fig .lk3 { fill: #c9c9c9; } .fig .lk4 { fill: #5a9bd4; } .fig .lk5 { fill: #1f5fa6; }
.fig .lk { fill: #fff; font: 600 12px var(--font); }
.part { padding: 2.5rem 0 .5rem; text-align: center; border-top: 1px solid var(--line); }
.part h2 { text-align: center; }
.toc { display: flex; flex-wrap: wrap; gap: .6rem; justify-content: center; margin-top: 1.2rem; }
.toc a { text-decoration: none; border: 1px solid var(--line); border-radius: 99px; padding: .25rem .9rem; background: var(--surface); font-size: .9rem; }
.data { color: var(--muted); font-size: .9rem; text-align: center; margin-top: 1rem; }
"""


def wrap_cmd(command, width=84):
    out, line = [], ""
    for word in command.split(" "):
        if line and len(line) + 1 + len(word) > width:
            out.append(line + " \\")
            line = "    " + word
        else:
            line = (line + " " + word) if line else word
    out.append(line)
    return "\n".join(out)


def term(command, output):
    command = wrap_cmd(command)
    return '<div class="term"><div class="bar"><i></i><i></i><i></i> a real run</div><pre><span class="cm">$ %s</span>\n%s</pre></div>' % (esc(command), esc(output))


def story(n, title, level, verbs, said, parts):
    body = "".join(parts)
    badges = '<span class="lvl">%s</span>' % esc(level) + "".join("<span>%s</span>" % esc(v) for v in verbs)
    return '<section class="%s story" id="s%d"><div class="wrap"><h2>%d &middot; %s</h2><div class="badge">%s</div><p class="said">%s</p>%s</div></section>' % ("alt" if n % 2 == 0 else "", n, n, esc(title), badges, said, body)


def build(binary):
    rows = make_bakery()
    files = {"bakery.csv": csv_text(rows, COLS)}
    a_first = [r for r in rows if r["oven"] == "A"][:4]
    b_first = [r for r in rows if r["oven"] == "B"][:4]
    files["pilot.csv"] = csv_text(a_first + b_first, ["oven", "bake_min"])
    with tempfile.TemporaryDirectory() as t:
        folder = pathlib.Path(t)
        for name, text in files.items():
            (folder / name).write_text(text)
        R = Runner(binary, folder)
        sections = []
        # 1 -------------------------------------------------------------------------------------------------
        args = ["describe", "--by", "oven", "--cols", "bake_min", "--bins", "6"]
        d, raw = R.run(args)
        g = {grp["value"]: grp["columns"][0] for grp in d["groups"]}
        excerpt = pretty({"by": d["by"], "groups": [{"value": k, "bake_min": {x: v[x] for x in ("count", "mean", "sd", "skewness", "quantiles")}} for k, v in g.items()]})
        sk = g["B"]["skewness"]
        qm = lambda c: [x["value"] for x in c["quantiles"] if x["p"] == "0.5"][0]
        claim(qm(g["B"]) < qm(g["A"]) and g["B"]["skewness"] > g["A"]["skewness"] and g["B"]["max"] > g["A"]["max"], "story 1: B has the lower median, the larger skewness and the slower slowest batch")
        sections.append(story(1, "Before anything else, look", "start here", ["describe"],
            "&#8220;I bought a second oven. Is it any good?&#8221; &#8212; Marta, who bakes for a street of cafes",
            ["<p>Marta kept a note of every batch for two months: which oven, how many minutes, how high it rose, how the crust scored, whether a loaf came out defective. Eighty batches, forty per oven. Before she asks anything, she looks at the one number she cares about most, the baking time.</p>",
             term(shown(args), excerpt), fig_boxes([("A", g["A"], "pa"), ("B", g["B"], "pb")]),
             '<div class="read"><b>Reading it.</b> The mean of oven A is %s minutes (sd %s); oven B&#8217;s is %s (sd %s). B&#8217;s median is lower, and its box sits left of A&#8217;s. But B has a long right whisker (its slowest batch took %g minutes) and a skewness of %s against %s for A. The sums behind these means are exact and rounded once, so the digits are the same on any machine. Something about oven B&#8217;s tail is worth a second look, later.</div>' % (num(g["A"]["mean"], 5), num(g["A"]["sd"], 4), num(g["B"]["mean"], 5), num(g["B"]["sd"], 4), g["B"]["max"], num(sk, 3), num(g["A"]["skewness"], 3))]))
        # 2 -------------------------------------------------------------------------------------------------
        args = ["test", "t2", "--cols", "bake_min", "--group", "oven"]
        d, raw = R.run(args)
        t = d["tests"][0]
        iv = t["interval"]
        args2 = ["ci", "diff", "--cols", "bake_min", "--group", "oven"]
        d2, _ = R.run(args2)
        iv2 = d2["intervals"][0]["interval"]
        claim(iv["low"] < 0 < iv["high"] and t["p_value"] > 0.05 and iv["estimate"] > 0, "story 2: A slower, the interval reaches across zero, p above 0.05")
        sections.append(story(2, "Is the new oven faster?", "one question", ["test t2", "ci diff"],
            "&#8220;Eyeballing it, B is quicker. Show me it is not my imagination.&#8221;",
            ["<p>The textbook answer is Welch&#8217;s t test: compare the two means, allow the spreads to differ. The tool prints the statistic, the degrees of freedom, the p-value, an effect size and the interval of the difference. The same interval is available alone with <code>ci diff</code>, byte for byte.</p>",
             term(shown(args), pretty({"levels": d["levels"], "tests": [{k: t[k] for k in ("statistic", "df", "p_value", "effect_size", "interval", "n")}]})),
             fig_interval(iv["estimate"], iv["low"], iv["high"], t["p_value"], -2, 7, "mean of A minus mean of B, in minutes"),
             '<div class="read"><b>Reading it.</b> A is %s minutes slower on average, with a 95%% interval from %s to %s and p = %s. The interval reaches across zero, so at the 5%% level this does not show a difference. (<code>ci diff</code> prints the same interval: %s &#8211; %s.) Marta is not convinced, and she is not wrong to be careful &#8212; but the pictures in story 1 suggested something. What is it that a mean cannot see?</div>' % (num(iv["estimate"], 3), num(iv["low"], 3), num(iv["high"], 3), num(t["p_value"], 3), num(iv2["low"], 3), num(iv2["high"], 3))]))
        # 3 -------------------------------------------------------------------------------------------------
        args = ["test", "mannwhitney", "--cols", "bake_min", "--group", "oven"]
        dm, _ = R.run(args)
        tm = dm["tests"][0]
        argsp = ["resample", "perm", "--cols", "bake_min", "--group", "oven", "--seed", "7", "--permutations", "9999"]
        dp, _ = R.run(argsp)
        tp = dp["tests"][0]
        pooled, _ = R.run(["test", "t2", "--cols", "bake_min", "--group", "oven", "--variance", "pooled"])
        pb = [r["bake_min"] for r in rows if r["oven"] == "B"]
        pa = [r["bake_min"] for r in rows if r["oven"] == "A"]
        claim(tm["p_value"] < 0.05 and t["p_value"] > 0.05 and pooled["tests"][0]["p_value"] > 0.05 and tp["p"] > 0.05 and len([v for v in pb if v > 53]) == 3 and not [v for v in pa if v > 53], "story 3: the rank test shows a difference, the mean-based tests do not, three batches of B ran past 53 minutes and none of A")
        sections.append(story(3, "Three loaves ran long", "when the data misbehave", ["test mannwhitney", "resample perm"],
            "&#8220;Three batches in oven B took over fifty minutes. The timer stuck, I think. I did not throw them away.&#8221;",
            ["<p>Every batch as a dot shows what the means hid: oven B&#8217;s typical batch is faster than A&#8217;s, and three batches &#8212; ringed &#8212; are far to the right. A mean is pulled by such points; a test that compares <em>ranks</em> is not. The Mann&#8211;Whitney test asks whether a batch from one oven tends to beat a batch from the other. The permutation test below asks the question of the means directly, by shuffling which batch went to which oven.</p>",
             fig_strip([("A", pa, "pa"), ("B", pb, "pb")]),
             term(shown(args), pretty({"tests": [{k: tm[k] for k in ("statistic", "p_value", "p_method", "effect_size", "n")}]})),
             term(shown(argsp), pretty({"p_method": dp["p_method"], "seed": dp["seed"], "tests": [{k: tp[k] for k in ("statistic", "p", "extreme", "permutations", "mc_error")}]})),
             fig_pbars([("Welch t test", t["p_value"], "pa"), ("pooled t test", pooled["tests"][0]["p_value"], "pa"), ("permutation, means", tp["p"], "warn"), ("Mann\u2013Whitney", tm["p_value"], "pb")]),
             '<div class="read"><b>Reading it.</b> The tests that use the <em>mean</em> (both t tests, and the permutation test of the mean difference) all say &#8220;not shown&#8221;, because three long batches inflate B&#8217;s mean and its spread. The rank test, which looks at the bulk of the batches, says p = %s. Neither is &#8220;the right one&#8221;: they answer different questions, and the answer says which one it answered. The permutation test\'s count is exact (%d of %d shuffles were as extreme; the Monte Carlo error is printed, %s), and every one of those numbers can be reproduced by an independent program from the seed. What Marta now knows: <em>the typical batch of B is faster; the average is dragged by three exceptions she can investigate.</em></div>' % (num(tm["p_value"], 3), tp["extreme"], tp["permutations"], num(tp["mc_error"], 2))]))
        # 4 -------------------------------------------------------------------------------------------------
        measures = ["bake_min", "height_cm", "crust", "weight_g", "moisture_pct", "salt_pct"]
        args = ["test", "t2", "--cols", ",".join(measures), "--group", "oven", "--correct", "holm"]
        d, raw = R.run(args)
        raw_p = {x["column"]: x["p_value"] for x in d["tests"]}
        adj_p = {x["column"]: pa_ for x, pa_ in zip(d["tests"], d["family"]["p_adjusted"])}
        claim({m for m in measures if adj_p[m] < 0.05} == {"height_cm", "moisture_pct"} and all(adj_p[m] >= raw_p[m] for m in measures), "story 4: height and moisture are the only measures that survive Holm")
        sections.append(story(4, "Six questions at once", "many tests, one decision", ["test t2", "--correct holm"],
            "&#8220;While I am at it: height, crust, weight, moisture, salt. Did the new oven change any of them?&#8221;",
            ["<p>Ask six questions at the 5% level and the chance that at least one comes out &#8220;significant&#8221; by luck alone is about one in four. The honest way to run a family of tests is to say so up front. <code>--correct holm</code> declares the six tests a family and prints each p-value both as computed and after Holm&#8217;s step-down correction; the tool refuses a family larger than <code>--max-tests</code>, so a fishing trip has to be declared.</p>",
             term(shown(args), pretty({"adjustment": d["adjustment"], "family": d["family"], "tests": [{"column": x["column"], "p_value": x["p_value"]} for x in d["tests"]]})),
             fig_adjust([(m, raw_p[m], adj_p[m]) for m in measures]),
             '<div class="read"><b>Reading it.</b> Height and moisture still stand out after the correction; the crust, the weight and the salt do not, and neither does the baking-time mean (the story 3 problem again). The hollow and filled circles are the same test before and after the family is accounted for. Marta can now say which of her improvements are real and which could be six coin flips.</div>' % ()]))
        # 5 -------------------------------------------------------------------------------------------------
        argsq = ["resample", "perm", "--cols", "bake_min", "--group", "oven", "--permutations", "all"]
        dq, _ = R.run(argsq, "pilot.csv")
        tq = dq["tests"][0]
        rc_big, err_big = refusal(binary, folder, argsq, "bakery.csv")
        argsf = ["test", "fisher", "--cols", "defect", "--vs", "oven"]
        df, _ = R.run(argsf)
        tf = df["tests"][0]
        a_n = sum(1 for r in rows if r["oven"] == "A")
        b_n = sum(1 for r in rows if r["oven"] == "B")
        a_yes = sum(1 for r in rows if r["oven"] == "A" and r["defect"] == "y")
        b_yes = sum(1 for r in rows if r["oven"] == "B" and r["defect"] == "y")
        pilot_a = [r["bake_min"] for r in a_first]
        pilot_b = [r["bake_min"] for r in b_first]
        claim(tf["p_value"] < 0.05 and tq["p"] > 0.5 and tq["extreme"] >= 1, "story 5: Fisher's p is below 0.05 and the pilot's exact p is large")
        sections.append(story(5, "Is it luck? Count it", "exact answers to small questions", ["resample perm --permutations all", "test fisher"],
            "&#8220;In the first week I only had four batches from each oven. Could four and four have told me anything?&#8221;",
            ["<p>With eight batches, &#8220;by luck&#8221; can be counted. Hand any four of the eight batches to oven A and the other four to oven B: there are %d ways. For each, compute the difference of the means; the p-value is the fraction of ways at least as far from zero as the real one. <code>--permutations all</code> does exactly that, and the comparison is made on exact integers, with no tolerance. The first four batches of each oven were %s minutes in A and %s in B: not obviously different, and the count will say how different &#8220;not obviously&#8221; is.</p>" % (tq["permutations"], ", ".join("%g" % v for v in pilot_a), ", ".join("%g" % v for v in pilot_b)),
             term(shown(argsq, "pilot.csv"), pretty({"p_method": dq["p_method"], "tests": [{k: tq[k] for k in ("statistic", "p", "extreme", "permutations", "n")}]})),
             fig_arrangements(tq["permutations"], tq["extreme"]),
             "<p>Why not do the same on all eighty batches? Because there are about 10<sup>23</sup> ways to split them, and the tool says so instead of trying:</p>",
             term(shown(argsq), "exit %d\n%s" % (rc_big, pretty({k: err_big[k] for k in ("code", "rule", "message", "hint", "detail")}))),
             '<div class="read"><b>Reading it.</b> %d of %d arrangements are as extreme as the real one (the red squares), so p = %s &#8212; the exact permutation p-value, not an approximation. Four against four cannot show a difference this small, and the count says so in the plainest way there is. There is no seed in this answer because nothing was drawn.</div>' % (tq["extreme"], tq["permutations"], num(tq["p"], 4)),
             "<p>The same idea decides the yes/no question. Over all eighty batches, %d of %d loaves from oven A had a defect and %d of %d from B. Fisher&#8217;s exact test enumerates every table with the same margins:</p>" % (a_yes, a_n, b_yes, b_n),
             term(shown(argsf), pretty({"tests": [{k: tf[k] for k in ("p_value", "p_method", "effect_size", "assumptions")}]})),
             fig_two_by_two(a_yes, a_n, b_yes, b_n),
             '<div class="read"><b>Reading it.</b> p = %s, from enumerating the tables, not from a formula that is only good for large counts: oven B has far fewer defective loaves, and this is the question on which the evidence is clearest. The odds ratio is printed too.</div>' % num(tf["p_value"], 3)]))
        # 6 -------------------------------------------------------------------------------------------------
        shuffled = list(rows)
        random.Random(99).shuffle(shuffled)
        reversed_rows = list(reversed(rows))
        alt = {"as recorded": rows, "rows shuffled": shuffled, "rows reversed": reversed_rows}
        q_args = ["test", "t2", "--cols", "bake_min,height_cm", "--group", "oven", "--correct", "holm"]
        digests = []
        for label, rr in alt.items():
            (folder / "variant.csv").write_text(csv_text(rr, COLS))
            _, rawtext = R.run(q_args, "variant.csv")
            digests.append((label + "  t2", hashlib.sha256(rawtext.encode()).hexdigest()[:48]))
        # exact permutation of the pilot, rows of the pilot in another order
        pil_rows = a_first + b_first
        pil_alt = {"as recorded": pil_rows, "rows shuffled": random.Random(3).sample(pil_rows, len(pil_rows)), "rows reversed": list(reversed(pil_rows))}
        pil_digests = []
        for label, rr in pil_alt.items():
            (folder / "variant.csv").write_text(csv_text(rr, ["oven", "bake_min"]))
            _, rawtext = R.run(argsq, "variant.csv")
            pil_digests.append((label + "  perm", hashlib.sha256(rawtext.encode()).hexdigest()[:48]))
        same_t = len({h for _, h in digests}) == 1
        same_p = len({h for _, h in pil_digests}) == 1
        sections.append(story(6, "Would anyone get the same numbers?", "the point of the tool", ["test t2", "resample perm", "sha256"],
            "&#8220;My accountant will run this on her laptop, and the cafe&#8217;s analyst on a server. Will they get my numbers?&#8221;",
            ["<p>The sums are exact, so the order in which rows are added cannot change a bit; the exact permutation test depends only on which values are in each group. Here is the proof you can run yourself: the same questions on the same data with the rows in three orders, and the SHA-256 of everything the tool printed.</p>",
             fig_hashes(digests + pil_digests),
             '<div class="read"><b>Reading it.</b> %s Across platforms it is the same: the project&#8217;s CI runs these answers on Linux x86-64, Linux arm64 and macOS arm64 and compares their hashes to the file <code>docs/xplat_fp.txt</code>. (One thing that <em>does</em> depend on row order, by design: a bootstrap resamples row positions, so <code>ci boot</code> is a function of the file&#8217;s bytes and the seed. The answer says so.)</div>' % ("All three digests of each question are identical." if same_t and same_p else "THE DIGESTS DIFFER: this is a bug, not a story.")]))
        if not (same_t and same_p):
            raise SystemExit("gen_stories: the digests differ across row orders")
        sections.append('<section class="part" id="lab"><div class="wrap"><h2>Part II &middot; In the lab and in the field</h2><p class="lead">The same tool on four small datasets from the sciences: a gene in mice, a chemistry method comparison, two particle detectors, a town survey. The data are invented; the runs are real; each story is built around a different habit of the field.</p></div></section>')
        from fractions import Fraction
        # 7 biology -----------------------------------------------------------------------------------------
        bio = make_bio()
        files["bio.csv"] = "genotype,expression\n" + "\n".join("%s,%s" % x for x in bio) + "\n"
        (folder / "bio.csv").write_text(files["bio.csv"])
        a1 = ["test", "t2", "--cols", "expression", "--group", "genotype"]
        a2 = ["test", "mannwhitney", "--cols", "expression", "--group", "genotype"]
        a3 = ["resample", "perm", "--cols", "expression", "--group", "genotype", "--permutations", "all"]
        d1, _ = R.run(a1, "bio.csv")
        d2, _ = R.run(a2, "bio.csv")
        d3, _ = R.run(a3, "bio.csv")
        t1, t2_, t3 = d1["tests"][0], d2["tests"][0], d3["tests"][0]
        claim(t1["p_value"] > 0.05 and t2_["p_method"] == "exact" and t2_["p_value"] < 0.05 and t3["p"] < 0.05 and t3["permutations"] == 924, "story 7: Welch does not show, the exact rank test and the exact permutation test do, 924 arrangements")
        wt = [v for g, v in bio if g == "wildtype"]
        ko = [v for g, v in bio if g == "knockout"]
        sections.append(story(7, "Six mice against six", "biology: small samples, exact answers", ["test t2", "test mannwhitney", "resample perm --permutations all"],
            "&#8220;Six wild-type and six knockout mice, one liver sample each. The knockout looks higher. One animal is very high. Is the gene doing something?&#8221; &#8212; a postdoc, late on a Friday",
            ["<p>Twelve animals is a small sample, and an assumption that the data are bell-shaped cannot be checked on twelve points. Three questions, from the most assumption-laden to the least: Welch&#8217;s t test (needs the means to be roughly normal), the Mann&#8211;Whitney test (ranks only; with no ties and so few animals the tool computes the <em>exact</em> p-value, <code>p_method: exact</code>, not the normal approximation), and a permutation test of the means over all %d ways to hand six of the twelve animals to the knockout group.</p>" % t3["permutations"],
             fig_dots([("wild-type", wt, "pa"), ("knockout", ko, "pb")], 6.5, 12.5, [7, 8, 9, 10, 11, 12], "log2 expression of the gene"),
             term(shown(a2, "bio.csv"), pretty({"tests": [{k: t2_[k] for k in ("statistic", "p_value", "p_method", "effect_size", "n")}]})),
             term(shown(a3, "bio.csv"), pretty({"p_method": d3["p_method"], "tests": [{k: t3[k] for k in ("statistic", "p", "extreme", "permutations", "n")}]})),
             fig_pbars([("Welch t test", t1["p_value"], "pa"), ("Mann\u2013Whitney, exact", t2_["p_value"], "pb"), ("permutation, exact", t3["p"], "pb")]),
             '<div class="read"><b>Reading it.</b> Welch&#8217;s p is %s: the one very high knockout animal inflates the knockout&#8217;s variance, and the test is cautious. The exact rank test gives p = %s and the exact permutation test p = %s (%d of %d arrangements as extreme). None of these is a Monte Carlo estimate: there is no seed because nothing was drawn. In a small experiment the honest sentence is the one with the method in it: <em>by an exact permutation test, p = %s</em>.</div>' % (num(t1["p_value"], 3), num(t2_["p_value"], 3), num(t3["p"], 3), t3["extreme"], t3["permutations"], num(t3["p"], 3))]))
        # 8 chemistry ---------------------------------------------------------------------------------------
        chem = make_chem()
        files["chem.csv"] = "sample,hplc,uv\n" + "\n".join("%d,%s,%s" % x for x in chem) + "\n"
        (folder / "chem.csv").write_text(files["chem.csv"])
        c1 = ["test", "paired", "--pairs", "hplc:uv"]
        c2 = ["ci", "corr", "--cols", "hplc,uv"]
        dc1, _ = R.run(c1, "chem.csv")
        dc2, _ = R.run(c2, "chem.csv")
        pt = dc1["tests"][0]
        pc = dc2["intervals"][0]
        claim(pc["estimate"] > 0.999 and pt["p_value"] < 0.05 and pt["interval"]["high"] < 0 and pc["interval"]["low"] > 0.99, "story 8: r above 0.999 and a paired interval that excludes zero")
        sections.append(story(8, "Two methods, one sample set", "chemistry: correlation is not agreement", ["test paired", "ci corr"],
            "&#8220;The UV method is ten times cheaper and it correlates almost perfectly with the HPLC. Can we switch?&#8221; &#8212; the lab manager",
            ["<p>Twelve caffeine samples were measured by both methods. A correlation coefficient asks whether the two move together; it says nothing about whether they give the <em>same number</em>. The tool answers both questions on the same file: <code>ci corr</code> with a Fisher-z interval for the correlation, and a paired test of the differences, whose mean is the bias between the methods.</p>",
             term(shown(c2, "chem.csv"), pretty({"intervals": [{k: pc[k] for k in ("x", "y", "n", "estimate", "interval")}]})),
             term(shown(c1, "chem.csv"), pretty({"tests": [{k: pt[k] for k in ("pair", "statistic", "df", "p_value", "interval", "n")}]})),
             fig_chem([(h, u) for _, h, u in chem], pt["interval"]["estimate"], pt["interval"]["low"], pt["interval"]["high"]),
             '<div class="read"><b>Reading it.</b> The correlation is %s, with a 95%% interval %s to %s: the methods move together almost perfectly. The paired differences have a mean of %s mg with a 95%% interval from %s to %s, which excludes zero (p = %s): the UV method reads systematically higher, here by %s per cent of the mean amount. A calibration against the HPLC would fix a bias like this; a correlation of 0.999 could never have shown it. The figure on the right is the same finding as a picture.</div>' % (num(pc["estimate"], 5), num(pc["interval"]["low"], 5), num(pc["interval"]["high"], 5), num(pt["interval"]["estimate"], 3), num(pt["interval"]["low"], 3), num(pt["interval"]["high"], 3), num(pt["p_value"], 3), num(-100 * pt["interval"]["estimate"] / (sum(h for _, h, _ in chem) / len(chem)), 2))]))
        # 9 physics -----------------------------------------------------------------------------------------
        phys = make_physics()
        files["detectors.csv"] = "detector,energy_keV\n" + "\n".join("%s,%s" % x for x in phys) + "\n"
        (folder / "detectors.csv").write_text(files["detectors.csv"])
        clock = make_clock()
        files["clock.csv"] = "t_s\n" + "\n".join(repr(v) for v in clock) + "\n"
        (folder / "clock.csv").write_text(files["clock.csv"])
        k1 = ["test", "t2", "--cols", "energy_keV", "--group", "detector"]
        k2 = ["test", "ks", "--cols", "energy_keV", "--group", "detector", "--exact", "yes"]
        dk1, _ = R.run(k1, "detectors.csv")
        dk2, _ = R.run(k2, "detectors.csv")
        tk1, tk2 = dk1["tests"][0], dk2["tests"][0]
        e1 = [v for g, v in phys if g == "det1"]
        e2 = [v for g, v in phys if g == "det2"]
        ecdf, dpy = fig_ecdf(e1, e2, tk2["statistic"])
        claim(tk1["p_value"] > 0.5 and tk2["p_value"] < 0.05 and tk2["p_method"].startswith("exact") and abs(dpy - tk2["statistic"]) < 1e-12, "story 9: equal means, a shape the KS test sees; the gap drawn is the statistic printed")
        k3 = ["describe", "--cols", "t_s"]
        dk3, _ = R.run(k3, "clock.csv")
        cv = dk3["columns"][0]
        fr = [Fraction(v) for v in clock]
        mean_fr = sum(fr, Fraction(0)) / len(fr)
        var_exact = float(sum(((x - mean_fr) ** 2 for x in fr), Fraction(0)) / (len(fr) - 1))
        naive = (sum(v * v for v in clock) / len(clock) - (sum(clock) / len(clock)) ** 2) * len(clock) / (len(clock) - 1)
        claim(cv["variance"] == var_exact and abs(naive - var_exact) > 1000 * var_exact, "story 9: the tool's variance is the exact one and the textbook one-pass formula in doubles is far from it")
        sections.append(story(9, "Same source, two detectors", "physics: shape, and digits that cancel", ["test t2", "test ks --exact yes", "describe"],
            "&#8220;Both detectors see the 662 keV line of the same caesium source. Their average energies agree. Why does the second one look broader?&#8221; &#8212; a detector physicist",
            ["<p>Thirty counts on each detector. The second detector&#8217;s gain drifted during the run, so its line is split into two close peaks. A comparison of means cannot see that; the two-sample Kolmogorov&#8211;Smirnov test compares the whole distributions, and with 60 values and no ties the tool computes its p-value exactly (<code>--exact yes</code> refuses if it cannot).</p>",
             term(shown(k1, "detectors.csv"), pretty({"tests": [{k: tk1[k] for k in ("statistic", "df", "p_value", "interval", "n")}]})),
             term(shown(k2, "detectors.csv"), pretty({"tests": [{k: tk2[k] for k in ("statistic", "p_value", "p_method", "n")}]})),
             ecdf,
             '<div class="read"><b>Reading it.</b> The t test sees nothing (p = %s): the two means differ by %s keV out of 662. The largest vertical gap between the two distribution functions is D = %s, and the exact p-value is %s. The drawn gap is the printed statistic.</div>' % (num(tk1["p_value"], 3), num(abs(tk1["interval"]["estimate"]), 2), num(tk2["statistic"], 4), num(tk2["p_value"], 3)),
             "<p>A second habit of the field: timestamps. Twenty arrivals of the same repeated event, recorded as seconds since 1970 &#8212; a number near 1.7&times;10<sup>9</sup> carrying microseconds of jitter. What is the spread? The textbook one-pass formula, <code>(&Sigma;x&sup2;/n &minus; mean&sup2;) n/(n&minus;1)</code>, is what many spreadsheets and quick scripts compute:</p>",
             term(shown(k3, "clock.csv"), pretty({"columns": [{k: cv[k] for k in ("count", "mean", "variance", "sd")}]})),
             '<div class="read"><b>Reading it.</b> The tool&#8217;s variance is %s s&sup2;, which is the exact variance of the twenty recorded numbers rounded once (checked against exact fractions when this page is generated). The one-pass formula in doubles, run on the same twenty numbers, answers %s: the two terms it subtracts are about 2.9&times;10<sup>18</sup> and agree to the last digit that a double holds, so what is left is rounding noise, here about 10<sup>%d</sup> times too large. The sums in <code>stats</code> are exact, so the offset costs nothing.</div>' % (num(cv["variance"], 6), num(naive, 6), int(round(__import__("math").log10(abs(naive / var_exact)))))]))
        # 10 social science ---------------------------------------------------------------------------------
        town = make_town()
        files["town.csv"] = "district,support,satisfaction,education\n" + "\n".join("%s,%s,%d,%s" % x for x in town) + "\n"
        (folder / "town.csv").write_text(files["town.csv"])
        s1 = ["ci", "prop", "--cols", "support", "--success", "yes"]
        s1b = s1 + ["--method", "clopper-pearson"]
        s2 = ["test", "kruskal", "--cols", "satisfaction", "--group", "district"]
        s3 = ["test", "chisq", "--cols", "support", "--vs", "education"]
        ds1, _ = R.run(s1, "town.csv")
        ds1b, _ = R.run(s1b, "town.csv")
        ds2, _ = R.run(s2, "town.csv")
        ds3, _ = R.run(s3, "town.csv")
        i1, i1b = ds1["intervals"][0], ds1b["intervals"][0]
        ts2, ts3 = ds2["tests"][0], ds3["tests"][0]
        mr = dict(zip(ds2["levels"], ts2["mean_rank"]))
        claim(mr["North"] < mr["Centre"] < mr["South"] and ts3["assumptions"]["expected_below_5"] == 0, "story 10: mean ranks rise from North to South and no expected count is below 5")
        claim(i1["interval"]["low"] < 0.5 and i1["interval"]["high"] > 0.5 and (i1b["interval"]["high"] - i1b["interval"]["low"]) > (i1["interval"]["high"] - i1["interval"]["low"]) and ts2["p_value"] < 0.01 and ts3["p_value"] > 0.2, "story 10: the interval of the share contains one half, Clopper-Pearson is wider than Wilson, Kruskal-Wallis shows a district effect and chi-square shows no education effect")
        likert = [(d, [sum(1 for r in town if r[0] == d and r[2] == lv) for lv in range(1, 6)]) for d in ("North", "Centre", "South")]
        sections.append(story(10, "What does the town think?", "social science: ordinal answers and honest intervals", ["ci prop", "test kruskal", "test chisq"],
            "&#8220;A hundred and fifty residents answered. Does the town support the new bus line? Do the districts feel the same?&#8221; &#8212; a council analyst",
            ["<p>Survey data are counts and ordinal scales, and the usual mistakes are treating a Likert score as a measurement on a ruler and quoting a share without its uncertainty. First the share: the Wilson interval is the default, and Clopper&#8211;Pearson (the exact binomial interval) is one flag away and always the wider.</p>",
             term(shown(s1, "town.csv"), pretty({"intervals": [{"column": i1["column"], "estimate": i1["estimate"], "interval": i1["interval"], "n": i1["n"]}]})),
             term(shown(s1b, "town.csv"), pretty({"intervals": [{"column": i1b["column"], "estimate": i1b["estimate"], "interval": i1b["interval"]}]})),
             fig_prop([("Wilson", i1["estimate"], i1["interval"]["low"], i1["interval"]["high"]), ("Clopper\u2013Pearson", i1b["estimate"], i1b["interval"]["low"], i1b["interval"]["high"])]),
             '<div class="read"><b>Reading it.</b> %s of the sample support the line. Both intervals contain one half, so &#8220;a majority supports it&#8221; is not something this sample can claim; the exact interval is wider, as it must be.</div>' % num(i1["estimate"], 3),
             "<p>Second, do the districts feel the same? Satisfaction is on a 1&#8211;5 scale: an order, not a ruler. The Kruskal&#8211;Wallis test uses only the order, and prints the number of ties and the mean rank of each district. Third, is support related to education? Counts in a table call for the chi-square test of independence; its answer carries the smallest expected count and how many are below 5, so that the reader can see the test applies.</p>",
             fig_likert(likert),
             term(shown(s2, "town.csv"), pretty({"levels": ds2["levels"], "tests": [{k: ts2[k] for k in ("statistic", "df", "p_value", "effect_size", "n", "mean_rank", "ties")}]})),
             term(shown(s3, "town.csv"), pretty({"tests": [{k: ts3[k] for k in ("statistic", "df", "p_value", "effect_size", "n", "assumptions")}]})),
             '<div class="read"><b>Reading it.</b> The districts differ clearly (Kruskal&#8211;Wallis p = %s): satisfaction rises from the North to the South, as the bars show. Support for the line shows no sign of depending on education (chi-square p = %s). Two different questions, two answers, and each with the test that suits the kind of data.</div>' % (num(ts2["p_value"], 3), num(ts3["p_value"], 3))]))
    return sections, files


def page(sections):
    index = (DOCS / "index.html").read_text()
    style = re.search(r"<style>(.*?)</style>", index, re.S).group(1)
    nav = ' '.join('<a class="opt" href="#s%d">%d</a>' % (i, i) for i in range(1, 11))
    toc = "".join('<a href="#s%d">%d &middot; %s</a>' % (i, i, esc(tt)) for i, tt in enumerate(("Look", "One question", "Outliers", "Six tests", "Exact counts", "Same bytes", "Mice", "Chemistry", "Detectors", "Survey"), 1))
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>cancho-stats: ten stories, from a bakery to the lab</title>
<meta name="description" content="Ten end-to-end examples of cancho-stats, in increasing complexity: a bakery, then biology, chemistry, physics and a social survey. Every number and every figure is drawn from a real run.">
<meta name="color-scheme" content="light dark">
<style>%s%s</style>
</head>
<body>
<header class="top"><div class="wrap">
  <a class="brand" href=".">cancho-stats</a>
  <nav>%s<a href=".">Overview</a><a href="https://github.com/alpibrusl/cancho-stats">GitHub</a></nav>
</div></header>
<div class="wrap hero small">
  <h1>Ten stories, from a bakery to the lab.</h1>
  <p class="sub">Marta has bought a second oven and wants to know whether it was a good idea; six stories follow her, each asking a little more of the tool than the last. Then four short ones from the sciences: biology, chemistry, physics and a social survey. Every number and every figure below is drawn from a real run of <code>stats</code> on the data you can download. The data are invented, with effects planted on purpose (B bakes faster and rises higher, three of its batches run long, the crust, weight and salt do not differ, A has more defects); the runs are real, and the narration is checked against the numbers whenever the page is generated.</p>
  <div class="toc">%s</div>
  <p class="data">Data (invented, seeded): <a href="stories/bakery.csv">bakery.csv</a> (80 batches), <a href="stories/pilot.csv">pilot.csv</a> (the first four batches of each oven), and for part II <a href="stories/bio.csv">bio.csv</a>, <a href="stories/chem.csv">chem.csv</a>, <a href="stories/detectors.csv">detectors.csv</a>, <a href="stories/clock.csv">clock.csv</a>, <a href="stories/town.csv">town.csv</a>. Rebuild this page: <code>python3 scripts/gen_stories.py</code>.</p>
</div>
%s
<footer><div class="wrap">cancho-stats &middot; EUPL-1.2 &middot; the data are invented; the runs, the digits and the figures are real &middot; <a href=".">overview</a></div></footer>
</body>
</html>
""" % (style, STORY_CSS, nav, toc, "\n".join(sections))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "build" / "stats"))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    sections, files = build(a.bin)
    out = {DOCS / "stories.html": page(sections)}
    for name, text in files.items():
        out[DOCS / "stories" / name] = text
    if a.check:
        stale = [str(p.relative_to(ROOT)) for p, text in out.items() if not p.exists() or p.read_text() != text]
        if stale:
            print("gen_stories: stale: " + ", ".join(stale))
            return 1
        print("gen_stories: up to date (%d files)" % len(out))
        return 0
    for p, text in out.items():
        p.parent.mkdir(exist_ok=True)
        p.write_text(text)
    print("gen_stories: wrote %d files" % len(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
