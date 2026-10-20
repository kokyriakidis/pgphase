"""The paper's references, written once: render.py makes paper/references.bib and the HTML reference list from
them. Every DOI is checked against Crossref by check_refs.py; entries without one are software or data.

Each entry: key -> dict(type, authors [family, given], title, venue fields, year, doi, cite (text citation)).
"""

def A(*names):
    """Authors as (family, given) pairs from 'Family, Given' strings; 'others' ends a truncated list."""
    return [tuple(n.split(", ", 1)) if ", " in n else (n, "") for n in names]


REFS = {
 # Read-based phasers compared here, and others
 "Patterson2015": dict(type="article", authors=A("Patterson, Murray", "Marschall, Tobias", "Pisanti, Nadia", "van Iersel, Leo",
     "Stougie, Leen", "Klau, Gunnar W.", "Schönhuth, Alexander"),
     title="WhatsHap: weighted haplotype assembly for future-generation sequencing reads", journal="Journal of Computational Biology",
     year="2015", volume="22", number="6", pages="498--509", doi="10.1089/cmb.2014.0157"),
 "Martin2016": dict(type="misc", authors=A("Martin, Marcel", "Patterson, Murray", "Garg, Shilpa", "Fischer, Sarah O.", "Pisanti, Nadia",
     "Klau, Gunnar W.", "Schönhuth, Alexander", "Marschall, Tobias"),
     title="WhatsHap: fast and accurate read-based phasing", journal="bioRxiv", year="2016", doi="10.1101/085050"),
 "Edge2017": dict(type="article", authors=A("Edge, Peter", "Bafna, Vineet", "Bansal, Vikas"),
     title="HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies", journal="Genome Research",
     year="2017", volume="27", number="5", pages="801--812", doi="10.1101/gr.213462.116"),
 "Lin2022": dict(type="article", authors=A("Lin, Jyun-Hong", "Chen, Liang-Chi", "Yu, Shu-Qi", "Huang, Yao-Ting"),
     title="LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants", journal="Bioinformatics",
     year="2022", volume="38", number="7", pages="1816--1822", doi="10.1093/bioinformatics/btac058"),
 "Holt2024": dict(type="article", authors=A("Holt, James M.", "Saunders, Christopher T.", "Rowell, William J.", "Kronenberg, Zev",
     "Wenger, Aaron M.", "Eberle, Michael"),
     title="HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing", journal="Bioinformatics",
     year="2024", volume="40", number="2", pages="btae042", doi="10.1093/bioinformatics/btae042"),
 "Gao2026": dict(type="misc", authors=A("Gao, Yan", "Liao, Wen-Wei", "Qin, Qian", "Hall, Ira M.", "Li, Heng"),
     title="LongcallD: joint calling and phasing of small, structural and mosaic variants from long reads", journal="bioRxiv",
     year="2026", doi="10.64898/2026.03.20.713111"),
 "Shafin2021": dict(type="article", authors=A("Shafin, Kishwar", "Pesout, Trevor", "Chang, Pi-Chuan", "others"),
     title="Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads",
     journal="Nature Methods", year="2021", volume="18", number="11", pages="1322--1332", doi="10.1038/s41592-021-01299-w"),
 "Poplin2018": dict(type="article", authors=A("Poplin, Ryan", "Chang, Pi-Chuan", "Alexander, David", "others"),
     title="A universal SNP and small-indel variant caller using deep neural networks", journal="Nature Biotechnology",
     year="2018", volume="36", number="10", pages="983--987", doi="10.1038/nbt.4235"),
 "Dempster1977": dict(type="article", authors=A("Dempster, Arthur P.", "Laird, Nan M.", "Rubin, Donald B."),
     title="Maximum likelihood from incomplete data via the EM algorithm",
     journal="Journal of the Royal Statistical Society: Series B (Methodological)", year="1977", volume="39", number="1",
     pages="1--22", doi="10.1111/j.2517-6161.1977.tb01600.x"),
 # Sequencing, references, pangenomes, graph alignment
 "Wenger2019": dict(type="article", authors=A("Wenger, Aaron M.", "Peluso, Paul", "Rowell, William J.", "others"),
     title="Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome",
     journal="Nature Biotechnology", year="2019", volume="37", number="10", pages="1155--1162", doi="10.1038/s41587-019-0217-9"),
 "Nurk2022": dict(type="article", authors=A("Nurk, Sergey", "Koren, Sergey", "Rhie, Arang", "others"),
     title="The complete sequence of a human genome", journal="Science", year="2022", volume="376", number="6588",
     pages="44--53", doi="10.1126/science.abj6987"),
 "Altemose2022": dict(type="article", authors=A("Altemose, Nicolas", "Logsdon, Glennis A.", "Bzikadze, Andrey V.", "others"),
     title="Complete genomic and epigenetic maps of human centromeres", journal="Science", year="2022", volume="376",
     number="6588", pages="eabl4178", doi="10.1126/science.abl4178"),
 "Liao2023": dict(type="article", authors=A("Liao, Wen-Wei", "Asri, Mobin", "Ebler, Jana", "others"),
     title="A draft human pangenome reference", journal="Nature", year="2023", volume="617", number="7960", pages="312--324",
     doi="10.1038/s41586-023-05896-x"),
 "Hickey2024": dict(type="article", authors=A("Hickey, Glenn", "Monlong, Jean", "Ebler, Jana", "others"),
     title="Pangenome graph construction from genome alignments with Minigraph-Cactus", journal="Nature Biotechnology",
     year="2024", volume="42", number="4", pages="663--673", doi="10.1038/s41587-023-01793-w"),
 "Garrison2018": dict(type="article", authors=A("Garrison, Erik", "Sirén, Jouni", "Novak, Adam M.", "others"),
     title="Variation graph toolkit improves read mapping by representing genetic variation in the reference",
     journal="Nature Biotechnology", year="2018", volume="36", number="9", pages="875--879", doi="10.1038/nbt.4227"),
 "Siren2021": dict(type="article", authors=A("Sirén, Jouni", "Monlong, Jean", "Chang, Xian", "others"),
     title="Pangenomics enables genotyping of known structural variants in 5202 diverse genomes", journal="Science",
     year="2021", volume="374", number="6574", pages="abg8871", doi="10.1126/science.abg8871"),
 "Siren2024": dict(type="article", authors=A("Sirén, Jouni", "Eskandar, Parsa", "Ungaro, Matteo Tommaso", "others"),
     title="Personalized pangenome references", journal="Nature Methods", year="2024", volume="21", number="11",
     pages="2017--2023", doi="10.1038/s41592-024-02407-2"),
 "Chang2025": dict(type="misc", authors=A("Chang, Xian", "Novak, Adam M.", "Eizenga, Jordan M.", "others"),
     title="Rapid, accurate long- and short-read mapping to large pangenome graphs with vg Giraffe", journal="bioRxiv",
     year="2025", doi="10.1101/2025.09.29.678807"),
 "Siren2022": dict(type="article", authors=A("Sirén, Jouni", "Paten, Benedict"),
     title="GBZ file format for pangenome graphs", journal="Bioinformatics", year="2022", volume="38", number="22",
     pages="5012--5018", doi="10.1093/bioinformatics/btac656"),
 "Paten2018": dict(type="article", authors=A("Paten, Benedict", "Eizenga, Jordan M.", "Rosen, Yohei M.", "Novak, Adam M.",
     "Garrison, Erik", "Hickey, Glenn"),
     title="Superbubbles, ultrabubbles, and cacti", journal="Journal of Computational Biology", year="2018", volume="25",
     number="7", pages="649--663", doi="10.1089/cmb.2017.0251"),
 "Li2020": dict(type="article", authors=A("Li, Heng", "Feng, Xiaowen", "Chu, Chong"),
     title="The design and construction of reference pangenome graphs with minigraph", journal="Genome Biology",
     year="2020", volume="21", pages="265", doi="10.1186/s13059-020-02168-z"),
 "Ebler2022": dict(type="article", authors=A("Ebler, Jana", "Ebert, Peter", "Clarke, Wayne E.", "others"),
     title="Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes",
     journal="Nature Genetics", year="2022", volume="54", number="4", pages="518--525", doi="10.1038/s41588-022-01043-w"),
 # Shared with the dinara-align paper (DOIs checked there)
 "Li2018": dict(type="article", authors=A("Li, Heng"), title="Minimap2: pairwise alignment for nucleotide sequences",
     journal="Bioinformatics", year="2018", volume="34", number="18", pages="3094--3100", doi="10.1093/bioinformatics/bty191"),
 "Gao2021": dict(type="article", authors=A("Gao, Yan", "Liu, Yongzhuang", "Ma, Yanmei", "Liu, Bo", "Wang, Yadong", "Xing, Yi"),
     title="abPOA: an SIMD-based C library for fast partial order alignment using adaptive band", journal="Bioinformatics",
     year="2021", volume="37", number="15", pages="2209--2211", doi="10.1093/bioinformatics/btaa963"),
 "Zook2016": dict(type="article", authors=A("Zook, Justin M.", "Catoe, David", "McDaniel, Jennifer", "others"),
     title="Extensive sequencing of seven human genomes to characterize benchmark reference materials",
     journal="Scientific Data", year="2016", volume="3", pages="160025", doi="10.1038/sdata.2016.25"),
 "Hansen2025": dict(type="article", authors=A("Hansen, Nancy F.", "Dwarshuis, Nathan", "Ji, Hyun Joo", "Rhie, Arang", "others"),
     title="A complete diploid human genome benchmark for personalized genomics", journal="Cell", year="2026",
     volume="189", number="16", pages="4857--4875.e31", doi="10.1016/j.cell.2026.06.016"),
 # Software and data without a paper
 "HipHap": dict(type="misc", authors=A("HipHap developers"),
     title="HipHap: the best alignment of each read to each haplotype of a diploid assembly", year="2026",
     url="https://github.com/jheinz27/hiphap"),
 "HPRC2": dict(type="misc", authors=A("Human Pangenome Reference Consortium"),
     title="Human pangenome reference, release 2 (Minigraph-Cactus graph, CHM13 reference)", year="2026",
     url="https://humanpangenome.org"),
}


def cite(key):
    """(names, year) as an author-year citation shows them: one or two families, else 'First et al.'."""
    r = REFS[key]
    fam = [f for f, _ in r["authors"] if f != "others"]
    many = len(r["authors"]) > 2 or any(f == "others" for f, _ in r["authors"])
    names = fam[0] if len(fam) == 1 and not many else (f"{fam[0]} and {fam[1]}" if not many else f"{fam[0]} et al.")
    return names, r["year"]


def bibtex():
    out = ["% Generated by tools/paper/render.py from refs.py; edit the source there, not this file.",
           "% Every DOI is checked against Crossref by tools/paper/check_refs.py.", ""]
    for key, r in REFS.items():
        def esc(s):
            return s.replace("&", "\\&")
        authors = " and ".join("others" if f == "others" else (f"{{{f}}}" if not g else f"{f}, {g}") for f, g in r["authors"])
        fields = [("author", authors), ("title", "{" + esc(r["title"]) + "}")]
        if r["type"] == "article":
            fields.append(("journal", esc(r["journal"])))
        elif r.get("journal"):
            fields.append(("howpublished", esc(r["journal"])))
        for f in ("year", "volume", "number", "pages", "doi", "url"):
            if r.get(f):
                fields.append((f, r[f]))
        body = ",\n".join(f"  {k:<12} = {{{v}}}" for k, v in fields)
        out.append(f"@{r['type']}{{{key},\n{body}\n}}\n")
    return "\n".join(out)


def html_entry(key):
    """Authors (year) Title. Venue, volume, pages. doi."""
    r = REFS[key]
    names = [("et al." if f == "others" else (f"{f}," + "".join(p[0] + "." for p in g.replace("-", " ").split()) if g else f))
             for f, g in r["authors"]]
    if len(names) > 1 and names[-1] == "et al.":
        authors = ", ".join(names[:-1]) + " et al."
    else:
        authors = ", ".join(names[:-1]) + (" and " if len(names) > 1 else "") + names[-1]
    venue = r.get("journal", "")
    for f, sep in (("volume", ", "), ("pages", ", ")):
        if r.get(f):
            venue += sep + r[f].replace("--", "–")
    return authors, r["title"] + ".", venue, r.get("doi"), r.get("url")
