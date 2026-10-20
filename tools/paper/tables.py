"""The paper's tables, figure data and every number its text quotes, from results.json (written by the
development repository's scripts/dev/paper_results.sh and paper_assemble.py). main_text.py and supplement.py
import TABLES, STABLES, V and POINTS from here, so no number in the paper is typed by hand.

    python3 tools/paper/tables.py    # prints V, the numbers the text quotes
"""
import json
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RES = json.load(open(os.path.join(HERE, "results.json")))

TOOLS = [("pgphase", "pgphase"), ("hiphase", "HiPhase"), ("longphase", "LongPhase"), ("whatshap", "WhatsHap"),
         ("whatshap_opt", "WhatsHap-dg"), ("longcalld", "longcallD")]
NAME = dict(TOOLS)
RIVALS = [t for t, _ in TOOLS if t != "pgphase"]
DEV = ["chr12", "chr17", "chr18", "chr19", "chr20"]
HELD = ["chr6", "chr16", "chr21"]
CHROMS = DEV + HELD
MHC = "chr6_mhc"
P = RES["profiles"]


def reads(p, tool, cls):
    return P[p][tool]["reads"][cls]


def total(tool, cls, chroms=CHROMS):
    out = {"correct": 0, "discordant": 0, "unphased": 0}
    for p in chroms:
        for k in out:
            out[k] += reads(p, tool, cls)[k]
    return out


def pct(part, whole, digits=2):
    return f"{100 * part / whole:.{digits}f}"


def n(x):
    return f"{x:,}"


def correct_pct(p, tool, cls="arms_phaseable"):
    r = reads(p, tool, cls)
    return 100 * r["correct"] / sum(r.values())


def bold_best(cells, values, lower=False):
    """Bold the best of `values` (several when tied)."""
    best = (min if lower else max)(values)
    return [f"**{c}**" if v == best else c for c, v in zip(cells, values)]


def census(tool, chroms=CHROMS):
    out = {"switches": 0, "real": 0, "representation": 0, "undecided": 0}
    for row in RES["census"]:
        if row["tool"] == tool and row["profile"] in chroms:
            for k in out:
                out[k] += row[k]
    return out


def var(tool, key, chroms=CHROMS):
    return sum(P[p][tool]["variants"][key] for p in chroms)


def label(p):
    return "MHC (chr6:28.5–34.0 Mb)" if p == MHC else p


# ------------------------------------------------------------------------------- Table: reads per chromosome
rows = [["Used in development", "", "", "", "", "", ""]]
for group, chroms in (("Used in development", DEV), ("Not used for tuning", HELD)):
    if group == "Not used for tuning":
        rows.append(["Not used for tuning", "", "", "", "", "", ""])
    for p in chroms:
        vals = [correct_pct(p, t) for t, _ in TOOLS]
        rows.append([p] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
rows.append(["Region of chr6", "", "", "", "", "", ""])
vals = [correct_pct(MHC, t) for t, _ in TOOLS]
rows.append([label(MHC)] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
rows.append(["All eight chromosomes", "", "", "", "", "", ""])
tot = {t: total(t, "arms_phaseable") for t, _ in TOOLS}
vals = [100 * tot[t]["correct"] / sum(tot[t].values()) for t, _ in TOOLS]
rows.append(["Correct (%)"] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
for k, lab in (("discordant", "Discordant (reads)"), ("unphased", "Unphased (reads)")):
    vals = [tot[t][k] for t, _ in TOOLS]
    rows.append([lab] + bold_best([n(v) for v in vals], vals, lower=True))
conf = {t: total(t, "arms_confident") for t, _ in TOOLS}
vals = [conf[t]["discordant"] for t, _ in TOOLS]
rows.append(["Discordant (reads), HapQ ≥ 30"] + bold_best([n(v) for v in vals], vals, lower=True))

TABLES = {}
TABLES["reads"] = {
    "caption": "Phaseable reads on the chromosome arms that each tool phases correctly. A read is phaseable when "
               "its alignment spans a heterozygous truth variant; it is correct when its haplotype tag agrees with "
               "the parental orientation its phase set takes over its truth-labelled reads, discordant when it "
               "does not, unphased when it has none. Per chromosome the share of phaseable reads phased correctly "
               "(%); below, the totals over the eight chromosomes, and the discordant reads whose parental assignment "
               "is confident (HapQ ≥ 30). The MHC lies inside chr6 and is not added to the totals. Best in bold. "
               "WhatsHap-dg: WhatsHap with `--distrust-genotypes`.",
    "header": ["Chromosome"] + [lab for _, lab in TOOLS],
    "align": "lrrrrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- Table: switches and variants
rows = []
for t, lab in TOOLS:
    c = census(t)
    rows.append([lab, n(c["switches"]), n(c["real"]), n(c["representation"]), n(c["undecided"]),
                 n(var(t, "switches_cen")), n(var(t, "flips_arms")), n(var(t, "truth_matched_phased_hets")),
                 n(sum(P[p][t]["blocks"]["phase_sets"] for p in CHROMS)),
                 f"{statistics.median(P[p][t]['blocks']['n50'] for p in CHROMS) / 1e3:,.0f}"])
for col, lower in ((1, True), (2, True), (5, True), (6, True), (7, False), (8, True), (9, False)):
    vals = [int(r[col].replace(",", "")) for r in rows]
    best = (min if lower else max)(vals)
    for r, v in zip(rows, vals):
        if v == best:
            r[col] = f"**{r[col]}**"
TABLES["switches"] = {
    "caption": "Switch and flip errors on the chromosome arms and switch errors in the centromeric and "
               "pericentromeric satellites, phased truth heterozygotes and blocks, over the eight chromosomes. A switch is a lasting change of a phase set's orientation against the truth between "
               "consecutive truth-matched heterozygotes, a flip a single heterozygote out of line with its "
               "neighbours. Each switch is classified by comparing the tool's two haplotype sequences around its two "
               "heterozygotes with the truth's (Section 2.8): *real*, a genuine phasing error, when the two windows "
               "pair the tool's haplotypes with the truth's in opposite ways; *representation* when both pair them "
               "the same way (the tool writes the same haplotypes as the truth with different records); "
               "*undecided* when either window's two pairings are within 30% "
               "of each other (or 2 edits). Phased hets: truth-matched heterozygous variants written phased. Phase "
               "sets summed over the chromosomes, block N50 the median of the per-chromosome N50s.",
    "header": ["Tool", "All", "Real", "Repr.", "Undec.", "Switches", "Arm flips", "Phased hets", "Phase sets", "N50 (kb)"],
    "groups": [("Arm switches", 1, 4), ("Satellites", 5, 5), ("Variants", 6, 7), ("Blocks", 8, 9)],
    "align": "lrrrrrrrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- Table: centromeres
rows = []
for t, lab in TOOLS:
    c = total(t, "cen_phaseable")
    rows.append([lab, n(c["correct"]), n(c["discordant"]), n(c["unphased"]), f"{c['correct'] / max(1, c['discordant']):.1f}"])
for col, lower in ((1, False), (2, True), (3, True), (4, False)):
    vals = [float(r[col].replace(",", "")) for r in rows]
    best = (min if lower else max)(vals)
    for r, v in zip(rows, vals):
        if v == best:
            r[col] = f"**{r[col]}**"
TABLES["cen"] = {
    "caption": "Phaseable reads in centromeric and pericentromeric satellites (one interval per chromosome, "
               "Supplementary Table S4) over the eight chromosomes, and the ratio of correct to discordant reads.",
    "header": ["Tool", "Correct", "Discordant", "Unphased", "Correct : discordant"],
    "align": "lrrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- Table: resources
RS = RES["resources"]
# The chromosomes every multi-threaded tool ran on with 6 threads (the rivals' chr12/18/20 runs used 20).
WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight"}
SAME_THREADS = [p for p in CHROMS if all(RS[p][t]["threads"] in (1, 6) for t in RS[p])]
# The phasers given DeepVariant's calls need DeepVariant run first; pgphase and longcallD find their own variants.
CALLS_OWN = {"pgphase", "longcalld"}
DV_TOTAL = sum(RS[p]["deepvariant"]["wall_s"] for p in SAME_THREADS)
rows = []
for t, lab in TOOLS:
    wall = sum(RS[p][t]["wall_s"] for p in SAME_THREADS)
    rss = max(RS[p][t]["rss_gb"] for p in SAME_THREADS)
    threads = sorted({RS[p][t]["threads"] for p in SAME_THREADS})
    with_calls = wall + (0 if t in CALLS_OWN else DV_TOTAL)
    rows.append([lab, "/".join(map(str, threads)), f"{RS['chr6'][t]['wall_s']:.0f}", f"{wall / 60:.1f}",
                 f"{with_calls / 60:.1f}", f"{rss:.1f}"])
vals = [float(r[4]) for r in rows]
for r, v in zip(rows, vals):
    if v == min(vals):
        r[4] = f"**{r[4]}**"
TABLES["resources"] = {
    "caption": "Wall-clock time and peak memory (GNU time) on one machine with 10 cores, on the "
               f"{WORDS[len(SAME_THREADS)]} chromosomes where every multi-threaded tool ran with 6 threads "
               f"({', '.join(SAME_THREADS)}): chr6, the largest, then the total time and the largest peak memory. "
               "Total: the tool alone. With calling: adding DeepVariant (6 threads), whose calls HiPhase, "
               "LongPhase and WhatsHap phase; pgphase and longcallD find the sample's variants themselves, so their "
               "two totals are the same. pgphase's time includes the haplotagged reads; no time includes read "
               "alignment, which every tool needs.",
    "header": ["Tool", "Threads", "chr6 (s)", "Total (min)", "With calling (min)", "Memory (GB)"],
    "align": "lrrrrr",
    "fit": True,
    "rows": rows,
}

# ------------------------------------------------------------------------------- supplementary tables
STABLES = {}
rows = []
for p in CHROMS + [MHC]:
    rows.append([label(p), "", "", "", "", "", "", ""])
    group = []
    for t, lab in TOOLS:
        r = reads(p, t, "arms_phaseable")
        c = reads(p, t, "arms_confident")
        v = P[p][t]["variants"]
        group.append([lab, r["correct"], r["discordant"], r["unphased"], c["discordant"],
                      v["switches_arms"], v["flips_arms"], v["truth_matched_phased_hets"]])
    # The best of each column in bold: most correct and phased hets, fewest of the rest.
    for col, lower in ((1, False), (2, True), (3, True), (4, True), (5, True), (6, True), (7, False)):
        best = (min if lower else max)(row[col] for row in group)
        for row in group:
            row[col] = f"**{n(row[col])}**" if row[col] == best else n(row[col])
    rows += group
STABLES["perchrom"] = {
    "label": "S1", "long": True,
    "caption": "Phaseable reads on the arms (correct, discordant, unphased), discordant reads whose parental assignment "
               "is confident (HapQ ≥ 30), arm switches and flips, and phased truth heterozygotes, per chromosome "
               "and tool.",
    "header": ["Tool", "Correct", "Discordant", "Unphased", "Disc. HapQ≥30", "Switches", "Flips", "Phased hets"],
    "align": "lrrrrrrr",
    "rows": rows,
}
rows = []
for p in CHROMS:
    rows.append([p, "", "", "", ""])
    for t, lab in TOOLS:
        c = census(t, [p])
        rows.append([lab, n(c["switches"]), n(c["real"]), n(c["representation"]), n(c["undecided"])])
STABLES["census"] = {
    "label": "S2", "long": True,
    "caption": "Arm switches per chromosome and tool, classified by the haplotype sequences around them (Section 2.8).",
    "header": ["Tool", "Switches", "Real", "Representation", "Undecided"],
    "align": "lrrrr",
    "rows": rows,
}
rows = []
for p in CHROMS:
    rows.append([p, "", "", ""])
    for t, lab in TOOLS:
        if t in RS.get(p, {}):
            r = RS[p][t]
            rows.append([lab, str(r["threads"]), f"{r['wall_s']:.0f}", f"{r['rss_gb']:.2f}"])
    if "deepvariant" in RS.get(p, {}):
        r = RS[p]["deepvariant"]
        rows.append(["DeepVariant (calls for HiPhase, LongPhase, WhatsHap)", str(r["threads"]), f"{r['wall_s']:.0f}", "–"])
STABLES["resources"] = {
    "label": "S3", "long": True,
    "caption": "Wall-clock time and peak memory per chromosome and tool, and DeepVariant's time (its memory was not "
               "measured: it ran in a container).",
    "header": ["Tool", "Threads", "Time (s)", "Peak memory (GB)"],
    "align": "lrrr",
    "rows": rows,
}

SC = RES["scoring"]
rows = []
for p in CHROMS + [MHC]:
    beg, end = SC[p]["cen"]
    excl = "; ".join(f"{c.split('#')[-1]}:{a / 1e6:.2f}–{b / 1e6:.2f} Mb" for c, a, b in SC[p]["exclude"]) or "none"
    rows.append([label(p), f"{beg / 1e6:.2f}–{end / 1e6:.2f}", excl])
STABLES["scoring"] = {
    "label": "S4",
    "caption": "Scoring settings per chromosome: the interval counted as centromeric and pericentromeric "
               "(CHM13 cenSat v2.1 satellite span on chr6, chr16, chr17, chr19 and chr21, the acrocentric short arm "
               "included on chr21; rounded spans on chr12, chr18 and chr20, set earlier in development), and the "
               "regions not scored because the truth contradicts itself there.",
    "header": ["Chromosome", "Centromeric interval (Mb)", "Not scored"],
    "align": "lll",
    "rows": rows,
}

# ------------------------------------------------------------------------------- parameters, from the source
import glob
import re as _re

_SRC = "".join(open(f).read() for f in sorted(glob.glob(os.path.join(HERE, "..", "..", "src", "*.[ch]pp"))))


def const(name):
    """The value of a constexpr in src/, which must be defined once (or several times to one value)."""
    vals = set(_re.findall(r"constexpr[\w:\s]+?\b" + name + r"\s*=\s*([0-9.e+-]+)", _SRC))
    assert len(vals) == 1, (name, vals)
    v = vals.pop()
    return f"{float(v):g}" if "." in v else f"{int(v):,}"


PARAMS = [
 ("Graph evidence and catalog sites", [
  ("kDefaultGraphMinMapq", "minimum mapping quality of graph alignments without a BAM", "2.2"),
  ("kUnionAlignmentMinMapq", "minimum mapping quality of graph and linear alignments with a BAM", "2.2"),
  ("kDefaultChunkSize", "chunk length without a BAM (bp)", "2.1"),
  ("kDefaultGraphBamChunkSize", "chunk length with a BAM (bp)", "2.1"),
  ("kDefaultMinDepth", "minimum reads at a site", "2.2"),
  ("kDefaultMinAltDepth", "minimum reads carrying the alternative allele", "2.2"),
  ("kDefaultMinAf", "minimum alternative allele fraction of a heterozygous site", "2.2"),
  ("kDefaultMaxAf", "maximum alternative allele fraction of a heterozygous site", "2.2"),
  ("kDefaultMinBaseq", "minimum base quality of an aligned base", "2.3"),
  ("kFamilywiseError", "family-wise error of the test for an absent reference base at catalog SNPs", "2.3"),
  ("kMinDeletionFraction", "deletion share of reads that keeps a catalog SNP from being reclassified", "2.3")]),
 ("The sample's own heterozygous variants", [
  ("kDefaultMinMapq", "mapping quality of a confident read", "2.4"),
  ("kUnionMinConfidentAlleleReads", "confident reads needed on each allele of an MSA-verified variant", "2.4"),
  ("kUnionMaxDepthRatio", "largest depth, relative to the chunk median, of an added site", "2.4"),
  ("kUnionMsaMaxRegionLength", "longest noisy region resolved by MSA (bp)", "2.4"),
  ("kUnionMaxComplementRowSpan", "longest catalog site whose calls only fill reads without a graph allele (bp)", "2.4"),
  ("kAddedReadMinMapq", "minimum mapping quality of reads only the linear alignment holds", "2.4"),
  ("kMarginErrorLambda", "slope \\(\\lambda\\) of the per-read error \\(1/(1+e^{\\lambda m})\\) of a realigned call", "2.4"),
  ("kMsaUndecidedError", "error of an MSA call that realignment cannot decide", "2.4"),
  ("kMsaContradictedError", "error of an MSA call that realignment contradicts", "2.4"),
  ("kGraphOnlyCallError", "error of a graph call at a catalog SNP whose base the alignment does not give", "2.4")]),
 ("Expectation-maximization", [
  ("kEmCleanStartError", "starting error of clean sites", "2.5"),
  ("kEmOtherStartError", "starting error of other sites", "2.5"),
  ("kEmMinError", "smallest site error", "2.5"),
  ("kEmMaxError", "largest site and call error", "2.5"),
  ("kEmClassPriorReads", "pseudo-reads of the class prior", "2.5"),
  ("kEmPoolMinWeight", "teaching weight (reads) for a site to enter its class pool", "2.5"),
  ("kEmPoolMinSites", "sites needed for the pool median (else the weighted mean)", "2.5"),
  ("kEmLearnMinMapq", "minimum mapping quality of a read that updates the model", "2.5"),
  ("kEmMinMinorAlleleFraction", "minimum share of the rarer allele at a site", "2.5"),
  ("kEmIterations", "iterations before switch moves", "2.5"),
  ("kEmSwitchRounds", "most rounds of switch moves", "2.5"),
  ("kEmUnreliableError", "learned error at which a site is unphased and labels no read", "2.5"),
  ("kEmReadLabelPosterior", "posterior for a read label", "2.5")]),
 ("Joins", [
  ("kEmBlockCutLogLikelihood", "join margin below which a block is cut (log units)", "2.6"),
  ("kEmMinJoinReads", "reads with evidence on both sides of a join", "2.6"),
  ("kEmWeakBridgeCutLogLikelihood", "margin needed by a join through weak sites (log units)", "2.6"),
  ("kEmMinStrongJoinReads", "reads at reliable, non-weak sites for a strong margin", "2.6"),
  ("kEmJoinMaxDepthRatio", "callers, relative to the chunk median, above which a site joins nothing", "2.6"),
  ("kEmDupDepthRatio", "callers, relative to the median, of a duplicated-haplotype site", "2.6"),
  ("kEmDupMinorFraction", "rarer-allele share below which such a site joins nothing", "2.6"),
  ("kEmJoinMinMeanMapq", "mean mapping quality of a site's callers below which it joins nothing", "2.6"),
  ("kEmEchoMaxGap", "distance within which a site can repeat another (bp)", "2.6"),
  ("kEmEchoMinShare", "share of reads two such sites must have in common", "2.6"),
  ("kEmEchoMinConcordance", "share of shared reads that must split the same way", "2.6"),
  ("kEmEchoMinReads", "minimum shared reads", "2.6")]),
 ("Labels for the remaining reads, and output", [
  ("kLocusMergeGap", "noisy calls merged into one local window within (bp)", "2.7"),
  ("kLocusPad", "padding of a local window (bp)", "2.7"),
  ("kLocusMaxPeriod", "longest tandem repeat period a window extends through (bp)", "2.7"),
  ("kLocusMaxRepeatScan", "longest repeat scan (bp)", "2.7"),
  ("kLocusMaxLength", "longest local window (bp)", "2.7"),
  ("kLocusMaxReads", "most reads per local window", "2.7"),
  ("kLocusMinMapq", "minimum mapping quality of a window's reads", "2.7"),
  ("kSeedLabelMinAgreement", "agreement of a read's phased sites for it to seed a haplotype", "2.7"),
  ("kLocusSeedSample", "seeded reads per haplotype medoid", "2.7"),
  ("kLocusSplitRounds", "rounds of the unseeded split", "2.7"),
  ("kLocusMinMargin", "edits by which a read must be closer to one haplotype or allele", "2.7"),
  ("kLocusMinSide", "reads needed on the rarer side of a window", "2.7"),
  ("kLocusMinMinorFraction", "share of reads needed on the rarer side", "2.7"),
  ("kLocusMaxObsError", "largest per-read error of a window call", "2.7"),
  ("kAlleleWindowFlank", "flank of an allele window beyond its repeat (bp)", "2.7"),
  ("kAlleleWindowMaxReads", "most reads per allele window", "2.7"),
  ("kPairMaxRefFraction", "largest reference share of an indel in a two-allele window", "2.7"),
  ("kLengthMinTract", "shortest tandem repeat tract of a length window (bp)", "2.7"),
  ("kLengthMinReads", "spanning reads needed to genotype a length window", "2.7"),
  ("kLengthMinHetLLR", "two lengths over one, for a heterozygous length window (log units)", "2.7"),
  ("kLengthMinReadLLR", "a read's least likelihood ratio between the two lengths", "2.7"),
  ("kStitchMinNetReads", "reads by which one orientation must lead to join chunks", "2.7"),
  ("kStitchMinShare", "share of reads that orientation needs", "2.7"),
  ("kRescueSiteOrientationPValue", "p-value of the association that orients an excluded site", "2.7"),
  ("kRescueSingletonMaxDiscordance", "upper 95% bound on discordance for a site oriented by one locus", "2.7"),
  ("kEmDeferredLabelPosterior", "posterior for a deferred label", "2.7"),
  ("kEmUnverifiedReportError", "largest learned error of an unverified call that labels reads", "2.7"),
  ("kEmUnverifiedLabelPosterior", "posterior for a label from unverified calls", "2.7"),
  ("kIslandOutsideMaxError", "largest disagreement of an island site with labelled reads", "2.7"),
  ("kIslandMinReads", "unlabelled reads an island site must call", "2.7"),
  ("kOneLocusSpan", "span of one locus for an unphased phase set (bp)", "2.7"),
  ("kEmVcfMaxSnpError", "learned error from which a SNP is written unphased", "2.7"),
  ("kEmVcfMaxIndelError", "learned error from which an indel is written unphased", "2.7")]),
]
rows = []
for group, items in PARAMS:
    rows.append([group, "", "", ""])
    for name, what, sec in items:
        rows.append([what[0].upper() + what[1:], const(name), sec, f"`{name}`"])
STABLES["params"] = {
    "label": "S5", "long": True,
    "caption": "Every threshold of pgphase, with its value as compiled (read from the source when this table is "
               "generated), the Methods section that uses it, and its name in the source.",
    "header": ["Parameter", "Value", "Section", "Name"],
    "align": "lrrl",
    "wrap": [0],
    "rows": rows,
}

# ------------------------------------------------------------------------------- figure data
# Per tool and chromosome, the discordant and unphased shares of the phaseable arm reads (%).
POINTS = {t: [(p, 100 * reads(p, t, "arms_phaseable")["discordant"] / sum(reads(p, t, "arms_phaseable").values()),
               100 * reads(p, t, "arms_phaseable")["unphased"] / sum(reads(p, t, "arms_phaseable").values()))
              for p in CHROMS] for t, _ in TOOLS}

# Per tool, over the eight chromosomes, the discordant and unphased shares of the phaseable arm reads (%).
TOTALS = {t: (100 * tot[t]["discordant"] / sum(tot[t].values()), 100 * tot[t]["unphased"] / sum(tot[t].values()))
          for t, _ in TOOLS}

# ------------------------------------------------------------------------------- numbers the text quotes
V = {}
pg = tot["pgphase"]
V["phaseable"] = n(sum(pg.values()))
V["pg_correct"] = n(pg["correct"])
V["pg_correct_pct"] = pct(pg["correct"], sum(pg.values()))
V["pg_correct_pct1"] = pct(pg["correct"], sum(pg.values()), 1)
V["phaseable_m"] = f"{sum(pg.values()) / 1e6:.1f}"
V["pg_disc_pct"] = pct(pg["discordant"], sum(pg.values()))
V["pg_unph_pct"] = pct(pg["unphased"], sum(pg.values()))
best_rival = max(RIVALS, key=lambda t: tot[t]["correct"])
V["best_rival"] = NAME[best_rival]
V["best_rival_correct_pct"] = pct(tot[best_rival]["correct"], sum(tot[best_rival].values()))
V["best_rival_correct_pct1"] = pct(tot[best_rival]["correct"], sum(tot[best_rival].values()), 1)
V["more_correct_than_best"] = n(pg["correct"] - tot[best_rival]["correct"])
fewest_err_rival = min(RIVALS, key=lambda t: tot[t]["discordant"] + tot[t]["unphased"])
V["fewest_err_rival"] = NAME[fewest_err_rival]
V["err_ratio"] = f"{(tot[fewest_err_rival]['discordant'] + tot[fewest_err_rival]['unphased']) / (pg['discordant'] + pg['unphased']):.1f}"
V["pg_conf_disc"] = n(conf["pgphase"]["discordant"])
low_conf_rival = min(RIVALS, key=lambda t: conf[t]["discordant"])
V["low_conf_rival"] = NAME[low_conf_rival]
V["low_conf_rival_disc"] = n(conf[low_conf_rival]["discordant"])
V["conf_ratio"] = f"{conf[low_conf_rival]['discordant'] / conf['pgphase']['discordant']:.1f}"
lowest_disc_rival = min(RIVALS, key=lambda t: tot[t]["discordant"])
V["lowest_disc_rival"] = NAME[lowest_disc_rival]
V["lowest_disc_rival_n"] = n(tot[lowest_disc_rival]["discordant"])
V["pg_disc"] = n(pg["discordant"])
V["pg_unph"] = n(pg["unphased"])
V["rows_best_correct"] = sum(1 for p in CHROMS if correct_pct(p, "pgphase") == max(correct_pct(p, t) for t, _ in TOOLS))
V["held_margin_min"] = f"{min(correct_pct(p, 'pgphase') - max(correct_pct(p, t) for t in RIVALS) for p in HELD):.1f}"
V["held_margin_max"] = f"{max(correct_pct(p, 'pgphase') - max(correct_pct(p, t) for t in RIVALS) for p in HELD):.1f}"
V["dev_margin_min"] = f"{min(correct_pct(p, 'pgphase') - max(correct_pct(p, t) for t in RIVALS) for p in DEV):.1f}"
V["dev_margin_max"] = f"{max(correct_pct(p, 'pgphase') - max(correct_pct(p, t) for t in RIVALS) for p in DEV):.1f}"
disc_lost = [p for p in CHROMS if reads(p, "pgphase", "arms_phaseable")["discordant"] > min(reads(p, t, "arms_phaseable")["discordant"] for t in RIVALS)]
by_tool = {}
for p in disc_lost:
    by_tool.setdefault(NAME[min(RIVALS, key=lambda t: reads(p, t, "arms_phaseable")["discordant"])], []).append(p)
V["disc_lost"] = "; ".join(f"on {' and '.join(ps)}, {tool} left fewer" for tool, ps in by_tool.items())
V["n_disc_lost"] = len(disc_lost)
V["mhc_pg_pct"] = f"{correct_pct(MHC, 'pgphase'):.2f}"
V["mhc_best_rival"] = NAME[max(RIVALS, key=lambda t: correct_pct(MHC, t))]
V["mhc_best_rival_pct"] = f"{max(correct_pct(MHC, t) for t in RIVALS):.2f}"
cs = {t: census(t) for t, _ in TOOLS}
V["pg_switches"] = cs["pgphase"]["switches"]
V["pg_real"] = cs["pgphase"]["real"]
V["pg_repr"] = cs["pgphase"]["representation"]
V["pg_undec"] = cs["pgphase"]["undecided"]
fewest_sw = min(RIVALS, key=lambda t: cs[t]["switches"])
V["fewest_sw_rival"] = NAME[fewest_sw]
V["fewest_sw_rival_n"] = cs[fewest_sw]["switches"]
V["fewest_sw_rival_real"] = cs[fewest_sw]["real"]
V["hiphase_switches"] = cs["hiphase"]["switches"]
V["hiphase_real"] = cs["hiphase"]["real"]
V["rival_real_min"] = min(cs[t]["real"] for t in RIVALS)
V["rival_real_max"] = max(cs[t]["real"] for t in RIVALS)
V["pg_cen_switches"] = n(var("pgphase", "switches_cen"))
V["rival_cen_switches_min"] = n(min(var(t, "switches_cen") for t in RIVALS))
V["rival_cen_switches_max"] = n(max(var(t, "switches_cen") for t in RIVALS))
assert var("pgphase", "switches_cen") < min(var(t, "switches_cen") for t in RIVALS)
V["pg_flips"] = n(var("pgphase", "flips_arms"))
fewest_fl = min(RIVALS, key=lambda t: var(t, "flips_arms"))
V["fewest_flip_rival"] = NAME[fewest_fl]
V["fewest_flip_rival_n"] = n(var(fewest_fl, "flips_arms"))
V["pg_hets"] = n(var("pgphase", "truth_matched_phased_hets"))
most_hets_rival = max(RIVALS, key=lambda t: var(t, "truth_matched_phased_hets"))
V["most_hets_rival"] = NAME[most_hets_rival]
V["most_hets_rival_n"] = n(var(most_hets_rival, "truth_matched_phased_hets"))
V["pg_ps"] = n(sum(P[p]["pgphase"]["blocks"]["phase_sets"] for p in CHROMS))
V["hiphase_ps"] = n(sum(P[p]["hiphase"]["blocks"]["phase_sets"] for p in CHROMS))
V["pg_n50_kb"] = f"{statistics.median(P[p]['pgphase']['blocks']['n50'] for p in CHROMS) / 1e3:,.0f}"
V["hiphase_n50_kb"] = f"{statistics.median(P[p]['hiphase']['blocks']['n50'] for p in CHROMS) / 1e3:,.0f}"
cen = {t: total(t, "cen_phaseable") for t, _ in TOOLS}
V["pg_cen_correct"] = n(cen["pgphase"]["correct"])
V["pg_cen_ratio"] = f"{cen['pgphase']['correct'] / cen['pgphase']['discordant']:.1f}"
most_cen_rival = max(RIVALS, key=lambda t: cen[t]["correct"])
V["most_cen_rival"] = NAME[most_cen_rival]
V["most_cen_rival_correct"] = n(cen[most_cen_rival]["correct"])
best_ratio = max(RIVALS, key=lambda t: cen[t]["correct"] / max(1, cen[t]["discordant"]))
V["cen_best_ratio_rival"] = NAME[best_ratio]
V["cen_best_ratio"] = f"{cen[best_ratio]['correct'] / cen[best_ratio]['discordant']:.1f}"
V["cen_best_ratio_correct"] = n(cen[best_ratio]["correct"])
V["cen_best_ratio_unph"] = n(cen[best_ratio]["unphased"])
V["most_cen_rival_ratio"] = f"{cen[most_cen_rival]['correct'] / cen[most_cen_rival]['discordant']:.1f}"
V["n_same_threads"] = WORDS[len(SAME_THREADS)]
V["pg_wall_total_min"] = f"{sum(RS[p]['pgphase']['wall_s'] for p in SAME_THREADS) / 60:.1f}"
V["hiphase_wall_total_min"] = f"{sum(RS[p]['hiphase']['wall_s'] for p in SAME_THREADS) / 60:.1f}"
V["longphase_wall_total_min"] = f"{sum(RS[p]['longphase']['wall_s'] for p in SAME_THREADS) / 60:.1f}"
V["pg_wall_chr6"] = f"{RS['chr6']['pgphase']['wall_s']:.0f}"
V["pg_rss_max"] = f"{max(RS[p]['pgphase']['rss_gb'] for p in SAME_THREADS):.1f}"
V["hiphase_wall_chr6"] = f"{RS['chr6']['hiphase']['wall_s']:.0f}"
V["longphase_wall_chr6"] = f"{RS['chr6']['longphase']['wall_s']:.0f}"
nohet = {t: total(t, "arms_nohet") for t, _ in TOOLS}
V["pg_nohet_disc"] = n(nohet["pgphase"]["discordant"])
V["pg_nohet_phased"] = n(nohet["pgphase"]["correct"] + nohet["pgphase"]["discordant"])
V["pg_real_chr19"] = census("pgphase", ["chr19"])["real"]
V["pg_switch_ties_mhc"] = all(P[MHC][t]["variants"]["switches_arms"] >= P[MHC]["pgphase"]["variants"]["switches_arms"] for t in RIVALS)
V["fewest_err_every_chrom"] = all(
    sum(reads(p, "pgphase", "arms_phaseable")[k] for k in ("discordant", "unphased")) <
    min(sum(reads(p, t, "arms_phaseable")[k] for k in ("discordant", "unphased")) for t in RIVALS) for p in CHROMS)
V["fewest_switch_every_chrom"] = all(
    P[p]["pgphase"]["variants"]["switches_arms"] < min(P[p][t]["variants"]["switches_arms"] for t in RIVALS) for p in CHROMS)
V["n_disc_best_word"] = WORDS[8 - len(disc_lost)]
DV = ["hiphase", "longphase", "whatshap", "whatshap_opt"]
best_dv = max(DV, key=lambda t: var(t, "truth_matched_phased_hets"))
V["best_dv_hets_rival"] = NAME[best_dv]
V["best_dv_hets_n"] = n(var(best_dv, "truth_matched_phased_hets"))
assert var("pgphase", "truth_matched_phased_hets") > var(best_dv, "truth_matched_phased_hets")
assert STABLES["scoring"]["label"] == "S4"  # the centromere caption cites it
folds = []
for p in CHROMS:
    pg = POINTS["pgphase"][CHROMS.index(p)]
    best = min(sum(POINTS[t][CHROMS.index(p)][1:]) for t in RIVALS)
    folds.append(best / sum(pg[1:]))
V["fold_min"], V["fold_max"] = f"{min(folds):.1f}", f"{max(folds):.1f}"
V["dv_total_min"] = f"{DV_TOTAL / 60:.0f}"
_pg_total = sum(RS[p]["pgphase"]["wall_s"] for p in SAME_THREADS)
_dv_based = {t: sum(RS[p][t]["wall_s"] for p in SAME_THREADS) + DV_TOTAL for t in ("hiphase", "longphase")}
V["hiphase_with_calls_min"] = f"{_dv_based['hiphase'] / 60:.0f}"
V["longphase_with_calls_min"] = f"{_dv_based['longphase'] / 60:.0f}"
V["calls_speedup"] = f"{min(_dv_based.values()) / _pg_total:.0f}"
V["commit"] = RES["meta"].get("git", "").split(" ")[0]

if __name__ == "__main__":
    for k, v in V.items():
        print(f"{k:28} {v}")
