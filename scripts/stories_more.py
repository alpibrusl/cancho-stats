"""Two more topics of docs/stories.html: medicine (a randomized trial and its side effects) and finance (what returns look like, and what a backtest can and cannot show). Invented seeded data; every number is drawn from a real run and `claim` stops the build when the narration stops matching them."""
import math
import random

from storylib import claim, pretty, esc, axis, svg, term, story
from stories_marketing import fig_rows, num, shown


# ------------------------------------------------------------------------------------------------ the invented data
def make_trial():
    """120 patients, alternately assigned to the drug or a placebo: change of systolic and diastolic pressure (mmHg) and of heart rate (beats a minute) after eight weeks"""
    r = random.Random(33)
    rows = []
    for i in range(1, 121):
        d = i % 2 == 0
        rows.append((i, "drug" if d else "placebo", round(r.gauss(-9 if d else -3, 10), 1), round(r.gauss(-4.5 if d else -2.5, 7), 1), round(r.gauss(-0.8 if d else 0, 6), 1)))
    return rows


def make_events():
    """the same 120 patients: did they report nausea or dizziness, and was there a serious adverse event (none, in this trial)"""
    r = random.Random(13)
    rows = []
    for i in range(1, 121):
        d = i % 2 == 0
        rows.append((i, "drug" if d else "placebo", "y" if r.random() < (0.15 if d else 0.04) else "n", "y" if r.random() < 0.07 else "n", "n"))
    return rows


def make_returns():
    """750 trading days of a share, in percent: heavy tails (a Student t with three degrees of freedom) and five crash days"""
    r = random.Random(5)
    out = []
    for i in range(750):
        t = r.gauss(0, 1) / math.sqrt(r.gammavariate(1.5, 2 / 3))
        out.append(round(0.09 + max(-3.0, min(3.0, t * 0.55)), 3))
    for k, v in zip((60, 190, 330, 470, 610), (-4.1, -3.6, -4.6, -3.9, -3.3)):
        out[k] = v
    return out


def make_strategies():
    """250 trading days: the index and six strategies that follow it closely and add noise; none has any edge"""
    r = random.Random(35)
    rows = []
    for _ in range(250):
        idx = r.gauss(0.04, 1.0)
        rows.append([round(idx, 3)] + [round(0.9 * idx + r.gauss(0, 0.45), 3) for _ in range(6)])
    return rows


def csv_rows(header, rows):
    return header + "\n" + "\n".join(",".join(str(x) for x in r) for r in rows) + "\n"


def fig_histogram(values, lo, hi, step, label, mark=None):
    """a histogram with a bar per bin, and a dashed line at `mark` (a value) if given"""
    bins = int(round((hi - lo) / step))
    counts = [0] * bins
    for v in values:
        k = min(bins - 1, max(0, int((v - lo) / step)))
        counts[k] += 1
    top = max(counts)
    ax, sx = axis(60, 610, 190, lo, hi, [t for t in range(int(lo), int(hi) + 1, 2)], label)
    body = [ax]
    for k, c in enumerate(counts):
        x0, x1 = sx(lo + k * step), sx(lo + (k + 1) * step)
        h = 150 * c / top
        body.append('<rect class="pbar pa" x="%g" y="%g" width="%g" height="%g"/>' % (x0 + 0.5, 190 - h, max(x1 - x0 - 1, 1), h))
    if mark is not None:
        body.append('<line class="zero" x1="%g" y1="20" x2="%g" y2="190"/>' % (sx(mark), sx(mark)))
    return svg(240, "\n".join(body), "Histogram of " + label)


def fig_pvalues(items):
    """p-values on a line from 0 to 1, with the 0.05 level marked; items: (label, raw p, adjusted p)"""
    ax, sx = axis(190, 610, 40 + 34 * len(items), 0, 1, [0, 0.05, 0.25, 0.5, 0.75, 1], "p-value (raw: hollow, after Holm: filled)")
    body = [ax, '<line class="zero" x1="%g" y1="20" x2="%g" y2="%g"/>' % (sx(0.05), sx(0.05), 40 + 34 * len(items))]
    for k, (name, raw, adj) in enumerate(items):
        y = 40 + 34 * k - 14
        body.append('<text class="nm" x="182" y="%g" text-anchor="end">%s</text>' % (y + 4, esc(name)))
        body.append('<line class="band" x1="%g" y1="%g" x2="%g" y2="%g"/>' % (sx(min(raw, adj)), y, sx(max(raw, adj)), y))
        body.append('<circle cx="%g" cy="%g" r="5" fill="none" stroke="currentColor" style="stroke:var(--accent);stroke-width:2"/>' % (sx(raw), y))
        body.append('<circle class="pt" cx="%g" cy="%g" r="5"/>' % (sx(adj), y))
    return svg(80 + 34 * len(items), "\n".join(body), "Raw and Holm-adjusted p-values against the 0.05 level")


# ------------------------------------------------------------------------------------------------ the stories
def stories(R, folder, files):
    sections = []
    # 15 ----------------------------------------------------------------------------------------------------
    trial = make_trial()
    files["trial.csv"] = csv_rows("patient,arm,sbp_change,dbp_change,hr_change", trial)
    (folder / "trial.csv").write_text(files["trial.csv"])
    a1 = ["test", "t2", "--cols", "sbp_change", "--group", "arm"]
    a2 = ["ci", "diff", "--cols", "sbp_change", "--group", "arm"]
    a3 = ["test", "t2", "--cols", "sbp_change,dbp_change,hr_change", "--group", "arm", "--correct", "holm"]
    d1, _ = R.run(a1, "trial.csv")
    d2, _ = R.run(a2, "trial.csv")
    d3, _ = R.run(a3, "trial.csv")
    t1, c2 = d1["tests"][0], d2["intervals"][0]
    ep = [(x["column"], x["p_value"], adj) for x, adj in zip(d3["tests"], d3["family"]["p_adjusted"])]
    claim(d1["levels"] == {"a": "drug", "b": "placebo"} and t1["p_value"] < 0.01 and c2["interval"]["high"] < 0, "story 15: the drug lowers the systolic pressure more than the placebo, with a clear interval")
    claim(ep[0][2] < 0.01 and ep[1][1] < 0.05 and ep[1][2] >= 0.05 and ep[2][1] > 0.3, "story 15: the diastolic p is below 0.05 alone and above it after Holm; the heart rate shows nothing")
    sections.append(story(15, "Did the drug lower the blood pressure?", "medicine: a randomized trial, three endpoints", ["test t2", "ci diff", "--correct holm"],
        "&#8220;Eight weeks, 120 patients, half on the drug and half on a placebo, assigned by a coin. The primary question is the systolic pressure. We also looked at the diastolic pressure and the heart rate.&#8221; &#8212; the trial statistician, before unblinding the report",
        ["<p>One question was written down before the trial: does the drug change the systolic pressure more than the placebo does? One test answers it, and an interval says by how much.</p>",
         term(shown(a1, "trial.csv"), pretty({"levels": d1["levels"], "tests": [{k: t1[k] for k in ("statistic", "df", "p_value", "effect_size", "n")}]})),
         term(shown(a2, "trial.csv"), pretty({"intervals": [{"estimate": c2["estimate"], "interval": {k: c2["interval"][k] for k in ("low", "high")}}]})),
         fig_rows([("drug minus placebo, mmHg", c2["estimate"], c2["interval"]["low"], c2["interval"]["high"])], -12, 2, [-12, -9, -6, -3, 0], "change in systolic pressure, mmHg (95% interval)", "no effect", "%.1f"),
         '<div class="read"><b>Reading it.</b> The drug lowers the systolic pressure by %s mmHg more than the placebo (95%% interval %s to %s; Welch&#8217;s t test p = %s). That is the pre-specified answer, and it does not need a correction: it is one question.</div>' % (num(-c2["estimate"], 1), num(-c2["interval"]["high"], 1), num(-c2["interval"]["low"], 1), num(t1["p_value"], 4)),
         "<p>The two other measurements are another matter. Three questions asked of the same patients, at the 5% level each, and the chance that at least one comes out significant by luck alone is about one in seven. Holm&#8217;s correction keeps that chance at 5% for the family: the smallest p-value is multiplied by the number of tests, the next by one fewer, and so on.</p>",
         term(shown(a3, "trial.csv"), pretty({"tests": [{"column": x["column"], "p_value": x["p_value"]} for x in d3["tests"]], "family": d3["family"]})),
         fig_pvalues([("systolic", ep[0][1], ep[0][2]), ("diastolic", ep[1][1], ep[1][2]), ("heart rate", ep[2][1], ep[2][2])]),
         '<div class="read"><b>Reading it.</b> The systolic result stands after the correction (p = %s). The diastolic pressure has p = %s on its own, which would pass 0.05, and %s after Holm, which does not: with three endpoints looked at, that result is what a family of tests would produce by chance one time in several, and the report should call it suggestive, not shown. The heart rate shows nothing (p = %s). Which questions were pre-specified, and how many were asked, is part of the result.</div>' % (num(ep[0][2], 4), num(ep[1][1], 3), num(ep[1][2], 3), num(ep[2][1], 2))]))
    # 16 ----------------------------------------------------------------------------------------------------
    events = make_events()
    files["events.csv"] = csv_rows("patient,arm,nausea,dizziness,serious", events)
    (folder / "events.csv").write_text(files["events.csv"])
    b1 = ["test", "fisher", "--cols", "nausea", "--vs", "arm"]
    b2 = ["test", "fisher", "--cols", "dizziness", "--vs", "arm"]
    b3 = ["ci", "prop", "--cols", "serious", "--success", "y", "--method", "clopper-pearson"]
    e1, _ = R.run(b1, "events.csv")
    e2, _ = R.run(b2, "events.csv")
    e3, _ = R.run(b3, "events.csv")
    f1, f2, p3 = e1["tests"][0], e2["tests"][0], e3["intervals"][0]
    n_drug = sum(1 for r_ in events if r_[1] == "drug")
    nau = {a: sum(1 for r_ in events if r_[1] == a and r_[2] == "y") for a in ("drug", "placebo")}
    diz = {a: sum(1 for r_ in events if r_[1] == a and r_[3] == "y") for a in ("drug", "placebo")}
    claim(f1["p_value"] < 0.05 and f2["p_value"] > 0.3 and nau["drug"] > nau["placebo"], "story 16: nausea is more frequent on the drug and Fisher's p is below 0.05; dizziness does not differ")
    claim(p3["n"]["successes"] == 0 and 0.02 < p3["interval"]["high"] < 0.04, "story 16: no serious event, and the exact interval reaches about 3%")
    sections.append(story(16, "No serious events: how safe is that?", "medicine: side effects and the zero that proves little", ["test fisher", "ci prop --method clopper-pearson"],
        "&#8220;Nausea was more common on the drug. And nobody had a serious adverse event, so the drug is safe.&#8221; &#8212; a sentence from a first draft of the safety section",
        ["<p>Side effects are counts, often small ones, and the arithmetic of small counts is exactly what Fisher&#8217;s exact test is for: it lists every table with the same margins and adds up those at least as unlikely as the observed one.</p>",
         term(shown(b1, "events.csv"), pretty({"tests": [{k: f1[k] for k in ("p_value", "p_method", "effect_size", "interval")}]})),
         term(shown(b2, "events.csv"), pretty({"tests": [{k: f2[k] for k in ("p_value", "p_method")}]})),
         '<div class="read"><b>Reading it.</b> Nausea: %d of %d on the drug against %d of %d on the placebo, Fisher p = %s. Dizziness: %d against %d, p = %s, no difference that 120 patients can show. The nausea result is a real signal to put in the label; the dizziness one is not a signal in either direction.</div>' % (nau["drug"], n_drug, nau["placebo"], n_drug, num(f1["p_value"], 4), diz["drug"], diz["placebo"], num(f2["p_value"], 3)),
         "<p>Then the sentence from the draft. Zero events in 120 patients. The exact (Clopper&#8211;Pearson) interval for a proportion says what rates are still compatible with seeing none:</p>",
         term(shown(b3, "events.csv"), pretty({"intervals": [{"estimate": p3["estimate"], "interval": {k: p3["interval"][k] for k in ("method", "low", "high")}, "n": p3["n"]}]})),
         fig_rows([("serious events, % of patients", 0.0, 0.0, p3["interval"]["high"] * 100)], -1, 5, [0, 1, 2, 3, 4, 5], "percent of patients with a serious adverse event (95% interval)", "none", "%.1f"),
         '<div class="read"><b>Reading it.</b> No serious event in %d patients is compatible with a true rate of up to %s%%, about three in a hundred (the &#8220;rule of three&#8221;: three divided by the number of patients). A serious event that happens to one patient in a hundred could easily go unseen in a trial this size. The honest sentence is &#8220;no serious adverse event was observed; the data cannot exclude a rate up to %s%%&#8221;, and a drug taken by millions needs a monitoring plan for what 120 patients cannot see.</div>' % (p3["n"]["n"], num(p3["interval"]["high"] * 100, 1), num(p3["interval"]["high"] * 100, 1))]))
    # 17 ----------------------------------------------------------------------------------------------------
    rets = make_returns()
    files["returns.csv"] = csv_rows("ret_pct", [(x,) for x in rets])
    (folder / "returns.csv").write_text(files["returns.csv"])
    r1 = ["describe", "--cols", "ret_pct", "--bins", "12"]
    r2 = ["ci", "mean", "--cols", "ret_pct"]
    r3 = ["ci", "boot", "--cols", "ret_pct", "--statistic", "quantile", "--prob", "0.01", "--method", "bca", "--seed", "7", "--replicates", "2000"]
    g1, _ = R.run(r1, "returns.csv")
    g2, _ = R.run(r2, "returns.csv")
    g3, _ = R.run(r3, "returns.csv")
    dc = g1["columns"][0]
    m2 = g2["intervals"][0]["interval"]
    q3 = g3["intervals"][0]
    normal1 = dc["mean"] - 2.326 * dc["sd"]
    claim(dc["kurtosis"] > 3 and dc["skewness"] < 0 and m2["low"] < 0 < m2["high"] and q3["interval"]["high"] < normal1, "story 17: heavy tails and a negative skew; the mean's interval contains zero; the whole interval of the 1% quantile lies below the normal model's")
    sections.append(story(17, "How bad is a bad day?", "finance: returns are not normal", ["describe", "ci boot --statistic quantile", "ci mean"],
        "&#8220;The risk report says the daily volatility is 0.9%, so a really bad day, one in a hundred, is a loss of about 2.0%. And our average day is positive.&#8221; &#8212; a risk analyst, rereading the report with some suspicion",
        ["<p>750 trading days of a share, in percent. Before any model, the shape: the skewness and the excess kurtosis (zero for a normal distribution) say how far from a bell this is.</p>",
         term(shown(r1, "returns.csv"), pretty({"columns": [{k: dc[k] for k in ("count", "mean", "sd", "min", "max", "skewness", "kurtosis")}]})),
         fig_histogram(rets, -6, 8, 1, "daily return, percent", 0),
         '<div class="read"><b>Reading it.</b> The excess kurtosis is %s (a normal distribution has 0), and the skewness is %s: more days far into the tails than a bell has, and the long tail is on the loss side (the worst day is %s%%, %s standard deviations from the mean, a day that a bell with a volatility of %s%% would make rarer than one in a million). It happened.</div>' % (num(dc["kurtosis"], 1), num(dc["skewness"], 2), num(dc["min"], 1), num(abs(dc["min"] - dc["mean"]) / dc["sd"], 1), num(dc["sd"], 2)),
         "<p>The average day first, because &#8220;our average day is positive&#8221; is the sentence people repeat. The interval for the mean, then the question that matters for risk: how bad is a really bad day, the loss that is exceeded one trading day in a hundred (the 1% quantile)? The bootstrap needs no assumption about the shape, and BCa corrects its bias and skew; the seed makes it repeatable.</p>",
         term(shown(r2, "returns.csv"), pretty({"intervals": [{"estimate": g2["intervals"][0]["estimate"], "interval": {k: m2[k] for k in ("low", "high")}}]})),
         term(shown(r3, "returns.csv"), pretty({"intervals": [{"statistic": g3["statistic"], "p": g3["p"], "estimate": q3["estimate"], "interval": {k: q3["interval"][k] for k in ("method", "low", "high")}}]})),
         fig_rows([("mean daily return, %", g2["intervals"][0]["estimate"], m2["low"], m2["high"]), ("1% quantile, % (a bad day)", q3["estimate"], q3["interval"]["low"], q3["interval"]["high"])], -4, 1, [-4, -3, -2, -1, 0, 1], "daily return, percent (95% intervals)", "zero", "%.2f"),
         '<div class="read"><b>Reading it.</b> The mean daily return is %s%% and its interval, %s to %s, contains zero: after 750 days the data do not show that the average day is positive. The one-in-a-hundred bad day is %s%%, with a 95%% interval from %s to %s: one day in a hundred loses more than about %s%%, and the whole interval lies beyond the %s%% that the bell of the report gives. The mean is hard to know and the tail is better known (and worse than the model said): a risk report should lead with the second.</div>' % (num(g2["intervals"][0]["estimate"], 3), num(m2["low"], 3), num(m2["high"], 3), num(q3["estimate"], 2), num(q3["interval"]["low"], 2), num(q3["interval"]["high"], 2), num(-q3["estimate"], 1), num(-normal1, 1))]))
    # 18 ----------------------------------------------------------------------------------------------------
    strat = make_strategies()
    files["strategies.csv"] = csv_rows("index,s1,s2,s3,s4,s5,s6", strat)
    (folder / "strategies.csv").write_text(files["strategies.csv"])
    pairs = ",".join("s%d:index" % k for k in range(1, 7))
    s1 = ["test", "paired", "--pairs", pairs]
    s2 = ["test", "paired", "--pairs", pairs, "--correct", "holm"]
    h1, _ = R.run(s1, "strategies.csv")
    h2, _ = R.run(s2, "strategies.csv")
    raw = [x["p_value"] for x in h1["tests"]]
    adj = h2["family"]["p_adjusted"]
    best = min(range(6), key=lambda k: raw[k])
    claim(sum(1 for p in raw if p < 0.05) == 1 and raw[best] < 0.03 and min(adj) > 0.1, "story 18: exactly one of six backtests passes 0.05 on its own, and none after Holm")
    bt = h1["tests"][best]
    sections.append(story(18, "Six backtests, one winner", "finance: the strategy that looked good by luck", ["test paired", "--correct holm"],
        "&#8220;We tried six variants of the strategy against the index over the last 250 days. Number six beats it, p = %s. We should trade it.&#8221; &#8212; a portfolio manager, with a chart" % num(raw[best], 3),
        ["<p>Each strategy and the index have a return on each day, so the comparison is paired: the same days, the difference of the two, and a test of whether its mean is zero. The data are invented and built so that <em>none</em> of the six has any edge: each follows the index and adds noise.</p>",
         term(shown(s1, "strategies.csv"), pretty({"tests": [{"pair": x["pair"], "p_value": x["p_value"], "mean_difference": x["interval"]["estimate"]} for x in h1["tests"]]})),
         '<div class="read"><b>Reading it.</b> One of the six, <code>s%d</code>, passes 0.05 on its own (p = %s, a mean advantage of %s percentage points a day). The chance that the best of six unrelated strategies passes 0.05 is about one in four, so this result is what a search over variants produces, not evidence about <code>s%d</code>.</div>' % (best + 1, num(raw[best], 3), num(bt["interval"]["estimate"], 3), best + 1),
         "<p>The remedy is to count every strategy that was tried, not only the one that was kept: Holm&#8217;s correction over the six.</p>",
         term(shown(s2, "strategies.csv"), pretty({"tests": [{"pair": x["pair"], "p_value": x["p_value"]} for x in h2["tests"]], "family": h2["family"]})),
         fig_pvalues([("s%d vs index" % (k + 1), raw[k], adj[k]) for k in range(6)]),
         '<div class="read"><b>Reading it.</b> After the correction the smallest p-value is %s and nothing is left: six backtests, no evidence of an edge in any. A backtest searched over variants has to carry the number of variants tried; &#8220;we tried six and kept the best&#8221; is the sentence that makes the p-value mean what it says. (Out-of-sample days, kept aside before the search, are the other cure, and cost data.)</div>' % num(min(adj), 3)]))
    return sections
