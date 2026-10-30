#!/usr/bin/env python3
"""The parts of the stories pages that do not depend on the tool: the terminal block, the figure frame and axes, the narration check, the page frame. cancho-stats, cancho-table and cancho-ml each carry the same copy of this
file (scripts/storylib.py) so that each repository builds its own page alone; change the three together."""
import html
import json
import re

W = 640


def claim(cond, what):
    """the narration says something about the numbers: if the numbers stop saying it, the page is not generated"""
    if not cond:
        raise SystemExit("gen_stories: the narration no longer matches the runs: " + what)


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
