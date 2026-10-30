"""Part III of docs/stories.html: marketing, and the one question under all of it, incrementality (what would have happened without the campaign?). Three stories on invented, seeded data: a randomized holdout, the dashboard that
would have said otherwise, and a geo experiment. As in gen_stories.py, every number is drawn from a real run, and `claim` stops the build when the narration stops matching them."""
import random

from storylib import claim, pretty, esc, axis, svg, term, story


def num(x, k):
    """k decimals, never an exponent"""
    return "%.*f" % (k, x)


BIG = {}  # the data files above the 2,000-line limit of a committed file: regenerated on demand (`gen_stories.py --data-out DIR`), never committed


# ------------------------------------------------------------------------------------------------ the invented data
def make_holdout():
    """8000 users, alternately shown the ads or held out (a coin the marketer controls); 6.2% and 5.0% convert; a converter spends a lognormal amount"""
    r = random.Random(11)
    rows = []
    for u in range(1, 8001):
        arm = "ads" if u % 2 == 0 else "holdout"
        c = 1 if r.random() < (0.062 if arm == "ads" else 0.05) else 0
        rows.append((u, arm, c, ("%.2f" % r.lognormvariate(3.7, 0.6)) if c else "0"))
    return rows


def make_observed():
    """6000 users the campaign was NOT randomized on: the busier a user, the likelier the ad reached them (ad targeting follows activity) and the likelier they buy anyway; the ad itself adds 1.2 points"""
    r = random.Random(5)
    rows = []
    for u in range(1, 6001):
        s = 1 + int(r.expovariate(1 / 3.2))
        e = r.random() < min(0.85, 0.08 + 0.08 * s)
        c = r.random() < min(0.9, 0.004 + 0.016 * s + (0.012 if e else 0))
        rows.append((u, "exposed" if e else "none", s, 1 if c else 0))
    return rows


def make_geo():
    """twelve sales regions, six drawn at random to get a TV campaign; weekly sales in thousands of euros the month before and the month after"""
    r = random.Random(15)
    arms = ["tv"] * 6 + ["none"] * 6
    r.shuffle(arms)
    rows = []
    for i, a in enumerate(arms, 1):
        b = round(r.gauss(120, 15), 1)
        ch = round(r.gauss(4.5 if a == "none" else 8.5, 2.8), 1)
        rows.append(("R%02d" % i, a, b, round(b + ch, 1), ch))
    return rows


def csv_rows(header, rows):
    return header + "\n" + "\n".join(",".join(str(x) for x in r) for r in rows) + "\n"


# ------------------------------------------------------------------------------------------------ figures
def fig_rows(items, xlo, xhi, ticks, xlabel, zero_label, fmt="%.1f"):
    """one interval per row; the first column of an item is its label, then estimate, low, high"""
    h = 80 + 34 * len(items)
    ax, sx = axis(190, 610, h - 40, xlo, xhi, ticks, xlabel)
    body = [ax, '<line class="zero" x1="%g" y1="30" x2="%g" y2="%g"/><text class="tk" x="%g" y="24" text-anchor="middle">%s</text>' % (sx(0), sx(0), h - 40, sx(0), esc(zero_label))]
    for k, (name, est, lo_, hi_) in enumerate(items):
        y = 54 + 34 * k
        body.append('<text class="nm" x="182" y="%g" text-anchor="end">%s</text>' % (y + 4, esc(name)))
        body.append('<line class="ciline" x1="%g" y1="%g" x2="%g" y2="%g"/><circle class="pt" cx="%g" cy="%g" r="5"/>' % (sx(lo_), y, sx(hi_), y, sx(est), y))
        body.append('<text class="tk" x="%g" y="%g" text-anchor="middle">%s</text>' % (sx(est), y - 9, esc(fmt % est)))
    return svg(h, "\n".join(body), "Estimates with their 95% intervals and the line of no effect")


def fig_slopes(rows):
    """each region's weekly sales before and after, one line per region"""
    lo, hi = 85, 170
    y0, y1 = 24, 214
    sy = lambda v: y1 - (v - lo) / (hi - lo) * (y1 - y0)
    body = ['<line class="ax" x1="200" y1="%g" x2="200" y2="%g"/><line class="ax" x1="480" y1="%g" x2="480" y2="%g"/>' % (y0, y1, y0, y1)]
    for t in (90, 110, 130, 150, 170):
        body.append('<text class="tk" x="192" y="%g" text-anchor="end">%d</text>' % (sy(t) + 4, t))
    for region, arm, b, a, _ in rows:
        cls = "pb" if arm == "tv" else "pa"
        body.append('<line class="wh %s" x1="200" y1="%g" x2="480" y2="%g"/>' % (cls, sy(b), sy(a)))
    body.append('<text class="lb" x="200" y="238" text-anchor="middle">month before</text><text class="lb" x="480" y="238" text-anchor="middle">month after</text>')
    body.append('<text class="nm" x="500" y="40" fill="currentColor"><tspan class="tk">green: TV campaign</tspan></text><text class="tk" x="500" y="58">blue: no campaign</text>')
    body.append('<text class="tk" x="90" y="120" text-anchor="middle" transform="rotate(-90 90 120)">weekly sales, k euros</text>')
    return svg(250, "\n".join(body), "Weekly sales of twelve regions before and after the campaign month")


# ------------------------------------------------------------------------------------------------ the stories
def stories(R, folder, files):
    sections = []
    # 12 -------------------------------------------------------------------------------------------------------
    holdout = make_holdout()
    BIG["holdout.csv"] = csv_rows("user,arm,converted,revenue", holdout)
    (folder / "holdout.csv").write_text(BIG["holdout.csv"])
    h1 = ["ci", "prop", "--cols", "converted", "--success", "1", "--group", "arm"]
    h2 = ["test", "prop", "--cols", "converted", "--success", "1", "--group", "arm"]
    h3 = ["ci", "diff", "--cols", "revenue", "--group", "arm"]
    h4 = ["resample", "perm", "--cols", "revenue", "--group", "arm", "--seed", "3", "--permutations", "9999"]
    d1, _ = R.run(h1, "holdout.csv")
    d2, _ = R.run(h2, "holdout.csv")
    d3, _ = R.run(h3, "holdout.csv")
    d4, _ = R.run(h4, "holdout.csv")
    i1, t2, i3, t4 = d1["intervals"][0], d2["tests"][0], d3["intervals"][0], d4["tests"][0]
    lift = i1["interval"]
    n_ads = i1["n"]["a"]["n"]
    claim(i1["levels"] == {"a": "ads", "b": "holdout"} and lift["low"] > 0 and t2["p_value"] < 0.01, "story 12: the conversion lift is positive, its interval excludes zero and the test is below 0.01")
    claim(i3["interval"]["low"] < 0 < i3["interval"]["high"] and 0.05 < t4["p"] < 0.1, "story 12: the revenue lift's interval contains zero and the permutation p is between 0.05 and 0.1")
    spend = 0.30 * n_ads
    inc_conv = i1["estimate"] * n_ads
    inc_rev = i3["estimate"] * n_ads
    claim(inc_rev > spend > i3["interval"]["low"] * n_ads, "story 12: the incremental revenue point estimate pays the spend, the low end of the interval does not")
    sections.append(story(12, "Did the ads cause the sales?", "marketing: the randomized holdout", ["ci prop --group", "test prop", "ci diff", "resample perm"],
        "&#8220;The dashboard says %s%% of the people we showed the ad to bought. Is that the ad?&#8221; &#8212; the head of growth, who has a budget meeting on Thursday" % num(100 * i1["n"]["a"]["successes"] / n_ads, 1),
        ["<p>The percentage on the dashboard cannot say. Some of those buyers would have bought anyway; what the ad is worth is the <em>extra</em> buyers, and the only way to see people who were never shown the ad is to withhold it from some of them at random. The team did that: 8000 users, alternately assigned (the assignment is the marketer&#8217;s coin, not the user&#8217;s choice) to see the ads or to be held out. A user row is: arm, bought (0 or 1), and revenue (0 for a user who did not buy).</p>",
         "<p>First the conversion rate in each arm and the difference between them, with the Newcombe interval for a difference of two proportions (the one that stays sensible when the rates are small):</p>",
         term(shown(h1, "holdout.csv"), pretty({"intervals": [{"estimate": i1["estimate"], "interval": {k: lift[k] for k in ("method", "level", "low", "high")}, "levels": i1["levels"], "n": i1["n"]}]})),
         term(shown(h2, "holdout.csv"), pretty({"tests": [{k: t2[k] for k in ("statistic", "p_value", "p_method", "effect_size")}]})),
         fig_rows([("conversion lift, points", i1["estimate"] * 100, lift["low"] * 100, lift["high"] * 100)], -1, 3, [-1, 0, 1, 2, 3], "extra buyers per 100 users shown the ad (95% interval)", "no effect", "%.2f"),
         '<div class="read"><b>Reading it.</b> %d of %d users shown the ad bought, against %d of %d held out: the ad adds %s points (the interval runs from %s to %s, and the test gives p = %s). The ad works. And the dashboard&#8217;s %s%% was misleading: %s of every 100 buyers would have bought anyway, and the incremental share is the %s points.</div>' % (
             i1["n"]["a"]["successes"], n_ads, i1["n"]["b"]["successes"], i1["n"]["b"]["n"], num(i1["estimate"] * 100, 2), num(lift["low"] * 100, 2), num(lift["high"] * 100, 2), num(t2["p_value"], 3),
             num(100 * i1["n"]["a"]["successes"] / n_ads, 1), num(100 * (i1["n"]["b"]["successes"] / i1["n"]["b"]["n"]) / (i1["n"]["a"]["successes"] / n_ads), 0), num(i1["estimate"] * 100, 2)),
         "<p>Now money. The budget meeting wants revenue, not conversions, and revenue per user is a nasty column: most users spend nothing and a few spend a lot. The difference of the mean revenue per user, with an interval; then the same question answered by shuffling the arm labels 9999 times, which assumes nothing about the shape of the column:</p>",
         term(shown(h3, "holdout.csv"), pretty({"intervals": [{"estimate": i3["estimate"], "interval": {k: i3["interval"][k] for k in ("low", "high")}, "n": i3["n"]}]})),
         term(shown(h4, "holdout.csv"), pretty({"tests": [{k: t4[k] for k in ("statistic", "p", "extreme", "permutations", "mc_error")}]})),
         fig_rows([("revenue per user, euros", i3["estimate"], i3["interval"]["low"], i3["interval"]["high"])], -0.5, 1.5, [-0.5, 0, 0.5, 1, 1.5], "extra revenue per user shown the ad, euros (95% interval)", "no effect", "%.2f"),
         '<div class="read"><b>Reading it.</b> The point estimate is %s euros extra per user: for the %d users shown the ad, %s euros of incremental revenue, against %s euros of media spend at an assumed 0.30 euros a user. That is a return of %s on each euro, if the estimate is right. But the interval for revenue includes zero (its low end is %s) and the shuffle test agrees (p = %s, give or take %s from the shuffling): at this sample size the lift in <em>money</em> is not distinguished from none, even though the lift in <em>buyers</em> is. The honest sentence for Thursday is &#8220;the ad brings buyers; whether it pays for itself needs more users or a longer test&#8221;. (The arithmetic on the printed numbers is the page&#8217;s, not the tool&#8217;s: the media cost is invented.)</div>' % (
             num(i3["estimate"], 2), n_ads, num(inc_rev, 0), num(spend, 0), num(inc_rev / spend, 2), num(i3["interval"]["low"], 2), num(t4["p"], 3), num(t4["mc_error"], 3))]))
    # 13 -------------------------------------------------------------------------------------------------------
    observed = make_observed()
    BIG["observed.csv"] = csv_rows("user,arm,sessions,converted", observed)
    (folder / "observed.csv").write_text(BIG["observed.csv"])
    for name, keep in (("quiet", lambda r: r[2] <= 3), ("busy", lambda r: r[2] > 3)):
        BIG["observed_%s.csv" % name] = csv_rows("user,arm,sessions,converted", [r for r in observed if keep(r)])
        (folder / ("observed_%s.csv" % name)).write_text(BIG["observed_%s.csv" % name])
    o1 = ["ci", "prop", "--cols", "converted", "--success", "1", "--group", "arm"]
    o2 = ["describe", "--by", "arm", "--cols", "sessions"]
    d1, _ = R.run(o1, "observed.csv")
    d2, _ = R.run(o2, "observed.csv")
    dq, _ = R.run(o1, "observed_quiet.csv")
    db, _ = R.run(o1, "observed_busy.csv")
    n1, nq, nb = d1["intervals"][0], dq["intervals"][0], db["intervals"][0]
    g = {grp["value"]: grp["columns"][0] for grp in d2["groups"]}
    claim(n1["levels"] == {"a": "exposed", "b": "none"} and n1["interval"]["low"] > lift["high"] and g["exposed"]["mean"] > 1.2 * g["none"]["mean"], "story 13: the naive lift is above the whole interval of the experiment and the exposed users are more active")
    claim(3.5 < n1["estimate"] / i1["estimate"] < 4.2, "story 13: the naive lift is nearly four times the experimental one")
    claim(lift["high"] < nq["estimate"] < n1["estimate"] and lift["high"] < nb["estimate"] < n1["estimate"], "story 13: both strata estimates are below the naive one and still above the experiment's whole interval")
    sections.append(story(13, "The dashboard that was too good to be true", "marketing: why exposed-versus-not is not a lift", ["ci prop --group", "describe --by"],
        "&#8220;We did not run a holdout for the autumn campaign, but we know who saw the ad and who did not. Why not just compare them?&#8221; &#8212; a marketer, with a spreadsheet",
        ["<p>The autumn campaign was targeted: the platform shows the ad to users who are active, because active users click. %d users, %d of them reached. Comparing buyers among the reached with buyers among the rest is the comparison every ad dashboard makes:</p>" % (n1["n"]["a"]["n"] + n1["n"]["b"]["n"], n1["n"]["a"]["n"]),
         term(shown(o1, "observed.csv"), pretty({"intervals": [{"estimate": n1["estimate"], "interval": {k: n1["interval"][k] for k in ("low", "high")}, "levels": n1["levels"], "n": n1["n"]}]})),
         "<p>A lift of %s points: nearly four times the %s of the holdout, and its whole interval lies above the whole interval of the experiment. One column explains why. How active were the two groups?</p>" % (num(n1["estimate"] * 100, 1), num(i1["estimate"] * 100, 1)),
         term(shown(o2, "observed.csv"), pretty({"by": d2["by"], "groups": [{"value": k, "sessions": {x: v[x] for x in ("count", "mean", "sd", "quantiles")}} for k, v in g.items()]})),
         '<div class="read"><b>Reading it.</b> The reached users had a mean of %s sessions and the others %s. Busy users are reached more often <em>and</em> buy more often with or without the ad, so the comparison counts their habit as the ad&#8217;s doing. This is selection, and no interval, however narrow, cures it: the interval above only measures how well the biased difference is known.</div>' % (num(g["exposed"]["mean"], 2), num(g["none"]["mean"], 2)),
         "<p>What the data can do, with no holdout, is compare like with like: split the users into quiet (3 sessions or fewer) and busy (4 or more), and compare reached with not reached inside each band (here the split files are given; with cancho-table it is a filter):</p>",
         term(shown(o1, "observed_quiet.csv"), pretty({"intervals": [{"estimate": nq["estimate"], "interval": {k: nq["interval"][k] for k in ("low", "high")}, "n": nq["n"]}]})),
         term(shown(o1, "observed_busy.csv"), pretty({"intervals": [{"estimate": nb["estimate"], "interval": {k: nb["interval"][k] for k in ("low", "high")}, "n": nb["n"]}]})),
         fig_rows([("everyone, exposed vs not", n1["estimate"] * 100, n1["interval"]["low"] * 100, n1["interval"]["high"] * 100), ("quiet users only", nq["estimate"] * 100, nq["interval"]["low"] * 100, nq["interval"]["high"] * 100),
                   ("busy users only", nb["estimate"] * 100, nb["interval"]["low"] * 100, nb["interval"]["high"] * 100), ("the holdout, above", i1["estimate"] * 100, lift["low"] * 100, lift["high"] * 100)],
                  -3, 9, [-2, 0, 2, 4, 6, 8], "difference in conversion, points (95% interval)", "no effect", "%.1f"),
         '<div class="read"><b>Reading it.</b> Inside the bands the lift falls from %s points to %s (quiet users) and %s (busy users). Comparing like with like moved the answer toward the experiment, and both estimates are still above the whole interval of the holdout (%s to %s): inside a band the busier users are still the ones reached, so some of the bias stays. That is the lesson. Adjustment narrows the gap and does not close it; only the experiment closes it. A campaign worth a large budget deserves a holdout, and a small one (even 5%% of the audience) costs far less than being wrong about the return.</div>' % (
             num(n1["estimate"] * 100, 1), num(nq["estimate"] * 100, 1), num(nb["estimate"] * 100, 1), num(lift["low"] * 100, 1), num(lift["high"] * 100, 1))]))
    # 14 -------------------------------------------------------------------------------------------------------
    geo = make_geo()
    files["geo.csv"] = csv_rows("region,arm,before_keur,after_keur,change_keur", geo)
    (folder / "geo.csv").write_text(files["geo.csv"])
    e1 = ["test", "t2", "--cols", "change_keur", "--group", "arm"]
    e2 = ["resample", "perm", "--cols", "change_keur", "--group", "arm", "--permutations", "all"]
    e3 = ["ci", "diff", "--cols", "change_keur", "--group", "arm"]
    d1, _ = R.run(e1, "geo.csv")
    d2, _ = R.run(e2, "geo.csv")
    d3, _ = R.run(e3, "geo.csv")
    w1, p2, c3 = d1["tests"][0], d2["tests"][0], d3["intervals"][0]
    claim(d2["levels"] == {"a": "none", "b": "tv"} and p2["statistic"] < 0 and p2["permutations"] == 924 and 0.02 < p2["p"] < 0.05 and w1["p_value"] < 0.05 and c3["interval"]["high"] < 0, "story 14: the TV regions grew more; exact p between 0.02 and 0.05 over 924 arrangements; Welch below 0.05; the interval excludes zero")
    sections.append(story(14, "Twelve regions and a TV campaign", "marketing: a geo experiment, difference in differences", ["test t2", "resample perm --permutations all", "ci diff"],
        "&#8220;We cannot hold out individual viewers of a TV spot. We can hold out whole regions.&#8221; &#8212; a media planner",
        ["<p>Television reaches a region, not a user, so the unit of the experiment is the region. Twelve sales regions; six were drawn at random to get the campaign for a month. Comparing the TV regions with the others <em>after</em> the month would mix the campaign with the fact that some regions are simply bigger. What does not depend on size is the <em>change</em> in each region: sales after minus sales before, in thousands of euros a week. Comparing changes between the two groups is the difference in differences, and here it is one column, <code>change_keur</code>, compared between two groups of six.</p>",
         fig_slopes(geo),
         "<p>Six against six is a small sample, so two answers are given, the t test and the exact one: the permutation test with <code>--permutations all</code> lists every one of the 924 ways to choose which six regions got the campaign, and counts those at least as extreme as the one drawn.</p>",
         term(shown(e1, "geo.csv"), pretty({"tests": [{k: w1[k] for k in ("statistic", "df", "p_value", "interval", "n")}]})),
         term(shown(e2, "geo.csv"), pretty({"p_method": d2["p_method"], "tests": [{k: p2[k] for k in ("statistic", "p", "extreme", "permutations", "n")}]})),
         term(shown(e3, "geo.csv"), pretty({"intervals": [{"estimate": c3["estimate"], "interval": {k: c3["interval"][k] for k in ("low", "high")}, "n": c3["n"]}]})),
         fig_rows([("TV effect, k euros", -c3["estimate"], -c3["interval"]["high"], -c3["interval"]["low"])], -1, 10, [0, 2, 4, 6, 8, 10], "extra weekly sales in a TV region, thousands of euros (95% interval)", "no effect", "%.1f"),
         '<div class="read"><b>Reading it.</b> The TV regions grew %s thousand euros a week more than the others (the tool prints the difference as none minus tv, hence the sign). The exact test says %d of the 924 arrangements are at least this extreme: p = %s, the t test p = %s. The interval runs from %s to %s thousand euros a week: the campaign probably worked, and could have worked a lot or barely. With six regions a side a planner should say exactly that, and the exact count is the safer of the two p-values to quote.</div>' % (
             num(-c3["estimate"], 1), p2["extreme"], num(p2["p"], 3), num(w1["p_value"], 3), num(-c3["interval"]["high"], 1), num(-c3["interval"]["low"], 1))]))
    return sections


def shown(args, file):
    import re
    return "stats " + " ".join(a if re.fullmatch(r"[\w.,:=/-]+", a) else "'%s'" % a for a in args) + " " + file
