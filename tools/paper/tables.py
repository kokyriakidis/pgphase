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
for group, chroms in (("Used in development", DEV), ("Held out", HELD)):
    if group == "Held out":
        rows.append(["Held out", "", "", "", "", "", ""])
    for p in chroms:
        vals = [correct_pct(p, t) for t, _ in TOOLS]
        rows.append([p] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
rows.append(["Region of chr6, held out", "", "", "", "", "", ""])
vals = [correct_pct(MHC, t) for t, _ in TOOLS]
rows.append([label(MHC)] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
rows.append(["All eight chromosomes", "", "", "", "", "", ""])
tot = {t: total(t, "arms_phaseable") for t, _ in TOOLS}
vals = [100 * tot[t]["correct"] / sum(tot[t].values()) for t, _ in TOOLS]
rows.append(["Correct (%)"] + bold_best([f"{v:.2f}" for v in vals], [round(v, 2) for v in vals]))
for k, lab in (("discordant", "Discordant"), ("unphased", "Unphased")):
    vals = [tot[t][k] for t, _ in TOOLS]
    rows.append([lab] + bold_best([n(v) for v in vals], vals, lower=True))
conf = {t: total(t, "arms_confident") for t, _ in TOOLS}
vals = [conf[t]["discordant"] for t, _ in TOOLS]
rows.append(["Discordant, HapQ ≥ 30"] + bold_best([n(v) for v in vals], vals, lower=True))

TABLES = {}
TABLES["reads"] = {
    "caption": "Phaseable reads on the chromosome arms that each tool phases correctly. A read is phaseable when "
               "its alignment spans a heterozygous truth variant; it is correct when its haplotype tag agrees with "
               "the parental orientation its phase set takes over its truth-labelled reads, discordant when it "
               "does not, unphased when it has none. Per chromosome the share of phaseable reads phased correctly "
               "(%); below, the totals over the eight chromosomes, and the discordant reads whose parental label "
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
                 n(var(t, "flips_arms")), n(var(t, "truth_matched_phased_hets")),
                 n(sum(P[p][t]["blocks"]["phase_sets"] for p in CHROMS)),
                 f"{statistics.median(P[p][t]['blocks']['n50'] for p in CHROMS) / 1e3:,.0f}"])
for col, lower in ((1, True), (2, True), (5, True), (6, False), (7, True), (8, False)):
    vals = [int(r[col].replace(",", "")) for r in rows]
    best = (min if lower else max)(vals)
    for r, v in zip(rows, vals):
        if v == best:
            r[col] = f"**{r[col]}**"
TABLES["switches"] = {
    "caption": "Switch and flip errors on the chromosome arms, phased truth heterozygotes and blocks, over the eight "
               "chromosomes. A switch is a lasting change of a phase set's orientation against the truth between "
               "consecutive truth-matched heterozygotes, a flip a single heterozygote out of line with its "
               "neighbours. Each switch is classed by comparing the tool's two haplotype sequences around its two "
               "heterozygotes with the truth's (Section 2.8): *real* when the two windows pair the haplotypes in "
               "opposite ways, *representation* when both pair them the same way (the tool writes the same "
               "haplotypes with different records), *undecided* when either window's two pairings are within 30% "
               "of each other. Phase sets summed, block N50 the median over the chromosomes.",
    "header": ["Tool", "Switches", "Real", "Repr.", "Undec.", "Flips", "Phased truth hets", "Phase sets", "N50 (kb)"],
    "align": "lrrrrrrrr",
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
    "caption": "Phaseable reads in centromeric and pericentromeric satellites (CHM13 cenSat v2.1 span; the "
               "acrocentric short arm of chr21 included) over the eight chromosomes, and the ratio of correct to "
               "discordant reads.",
    "header": ["Tool", "Correct", "Discordant", "Unphased", "Correct : discordant"],
    "align": "lrrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- Table: resources
RS = RES["resources"]
# The chromosomes every multi-threaded tool ran on with 6 threads (the rivals' chr12/18/20 runs used 20).
SAME_THREADS = [p for p in CHROMS if all(RS[p][t]["threads"] in (1, 6) for t in RS[p])]
rows = []
for t, lab in TOOLS:
    wall = sum(RS[p][t]["wall_s"] for p in SAME_THREADS)
    rss = max(RS[p][t]["rss_gb"] for p in SAME_THREADS)
    threads = sorted({RS[p][t]["threads"] for p in SAME_THREADS})
    rows.append([lab, "/".join(map(str, threads)), f"{RS['chr6'][t]['wall_s']:.0f}", f"{wall / 60:.1f}", f"{rss:.1f}"])
TABLES["resources"] = {
    "caption": "Wall-clock time and peak memory (GNU time) on one machine with 10 cores, on the "
               f"{len(SAME_THREADS)} chromosomes where every multi-threaded tool ran with 6 threads "
               f"({', '.join(SAME_THREADS)}): chr6, the largest, then the total time and the largest peak memory. "
               "The DeepVariant-based phasers' times exclude DeepVariant; pgphase's include its own variant "
               "discovery and the haplotagged reads, and exclude graph alignment.",
    "header": ["Tool", "Threads", "chr6 (s)", "Total (min)", "Memory (GB)"],
    "align": "lrrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- supplementary tables
STABLES = {}
rows = []
for p in CHROMS + [MHC]:
    rows.append([label(p), "", "", "", "", "", "", ""])
    for t, lab in TOOLS:
        r = reads(p, t, "arms_phaseable")
        c = reads(p, t, "arms_confident")
        v = P[p][t]["variants"]
        rows.append([lab, n(r["correct"]), n(r["discordant"]), n(r["unphased"]), n(c["discordant"]),
                     n(v["switches_arms"]), n(v["flips_arms"]), n(v["truth_matched_phased_hets"])])
STABLES["perchrom"] = {
    "label": "S1", "long": True,
    "caption": "Phaseable reads on the arms (correct, discordant, unphased), discordant reads whose parental label "
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
    "caption": "Arm switches per chromosome and tool, classed by the haplotype sequences around them (Section 2.8).",
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
STABLES["resources"] = {
    "label": "S3", "long": True,
    "caption": "Wall-clock time and peak memory per chromosome and tool.",
    "header": ["Tool", "Threads", "Time (s)", "Peak memory (GB)"],
    "align": "lrrr",
    "rows": rows,
}

# ------------------------------------------------------------------------------- figure data
# Per tool and chromosome, the discordant and unphased shares of the phaseable arm reads (%).
POINTS = {t: [(p, 100 * reads(p, t, "arms_phaseable")["discordant"] / sum(reads(p, t, "arms_phaseable").values()),
               100 * reads(p, t, "arms_phaseable")["unphased"] / sum(reads(p, t, "arms_phaseable").values()))
              for p in CHROMS] for t, _ in TOOLS}

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
V["n_same_threads"] = len(SAME_THREADS)
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
V["commit"] = RES["meta"].get("git", "").split(" ")[0]

if __name__ == "__main__":
    for k, v in V.items():
        print(f"{k:28} {v}")
