"""Where a mutant's pattern lives (scripts/mutants_*.py): the reader of the files moved from tools/stats/main.cho to lib/csvx/scan.cho and rule.cho (2,000 lines is the limit), and a mutant pattern is
found in whichever of the files that may hold it contains it. `where(root, pattern)` answers the first of SOURCES that has the text; the script then counts the occurrences in that file."""
SOURCES = ["tools/stats/main.cho", "lib/csvx/scan.cho", "lib/csvx/rule.cho", "lib/csvx/pick.cho", "lib/stest/answer.cho", "lib/stest/answer2.cho", "lib/stest/ttest.cho", "lib/stest/welch.cho", "lib/csvx/cross.cho", "lib/stest/chisq.cho", "lib/stest/answer3.cho", "lib/stest/prop.cho", "lib/stest/answer4.cho", "lib/stest/correct.cho", "lib/stest/anova.cho", "lib/stest/kruskal.cho", "lib/stest/wilcoxon.cho", "lib/stest/answer8.cho", "lib/stest/mannwhitney.cho", "lib/stest/exact.cho", "lib/stest/answer7.cho", "lib/stest/answer6.cho", "lib/stest/welchf.cho", "lib/stest/answer5.cho", "lib/csvx/run.cho"]


def where(root, pattern):
    for f in SOURCES:
        if pattern in (root / f).read_text():
            return f
    return SOURCES[0]


def relocate(root, M):
    """the 4-tuples (name, file, pattern, replacement) with the file changed to the one of SOURCES that has the pattern, where the named file does not"""
    out = []
    for n, f, o, w in M:
        if o not in (root / f).read_text():
            f = where(root, o)
        out.append((n, f, o, w))
    return out
