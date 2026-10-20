"""The main paper, written once; render.py turns it into paper/ (LaTeX) and build/paper/paper.html.

Inline markup: **bold**, *italic*, `code`, [@Key; @Key] parenthetical citation, [@!Key] narrative citation,
\\( \\) inline maths, \\[ \\] display maths, [#tab:x] / [#fig:x] cross-references.
Blocks: ("sec", title), ("sub", title), ("p", text), ("list", [items]) numbered, ("table", key), ("figure", key).
Every number comes from tables.V (results.json), never typed here.
"""
from tables import TABLES, V

TITLE = "pgphase: phasing long reads at a pangenome's known alleles and the sample's own heterozygous variants"
# Authors in order, each with an email. DRAFT: co-authors to be added.
AUTHORS = [("Konstantinos Kyriakidis", "kkyriaki@ucsc.edu")]
AFFILIATION = "UC Santa Cruz Genomics Institute, Santa Cruz, CA, USA"
KEYWORDS = ["phasing", "haplotagging", "long reads", "pangenome graph", "graph alignment", "expectation maximization"]

ABSTRACT = [
 ("Motivation",
  "Read-based phasers assign the heterozygous variants of a diploid genome, and the long reads that carry them, to "
  "the two parental haplotypes. They phase the variants a caller found on a linear reference: where the caller "
  "misses the heterozygous variants of a region its reads are left unphased, and where it calls the differences "
  "between copies of a collapsed duplication as heterozygous, they are phased against one another and join blocks "
  "in the wrong orientation. Applications that use the reads themselves, such as haplotype-resolved assembly and "
  "allele-specific analyses, depend on how many reads are phased, and phased correctly."),
 ("Results",
  "We present pgphase, which phases reads at the known alleles of a pangenome graph, observed through the reads' "
  "own graph alignments without calling them, together with the heterozygous variants the sample carries outside "
  "the graph, found in its linear alignments. One expectation-maximization over both learns each site's error "
  "rate, and blocks are joined only where several reads verify the join. On HG002 PacBio HiFi reads aligned to the "
  "Human Pangenome Reference Consortium graph, with the parent of origin of every read as truth, on eight "
  f"chromosomes of which three were held out from development, pgphase phased {V['pg_correct_pct']}% of "
  f"{V['phaseable']} phaseable reads correctly, against {V['best_rival_correct_pct']}% for the best of HiPhase, "
  f"LongPhase, WhatsHap and longcallD. It left {V['err_ratio']} times fewer reads misphased or unphased than "
  f"{V['fewest_err_rival']}, made {V['pg_switches']} switch errors against {V['fewest_sw_rival_n']} for the "
  f"best of them, and phased more truth heterozygous variants than any phaser given DeepVariant's calls."),
 ("Availability and implementation",
  "pgphase is free software under the MIT license, available at https://github.com/kokyriakidis/pgphase."),
]

BODY = [
('sec', 'Introduction'),
('p',
 'A long read spans several heterozygous variants, so the reads of a diploid sample and its heterozygous variants '
 'can be phased together: each read is assigned to the parental haplotype it was sequenced from, and each variant '
 'allele to the haplotype whose reads carry it [@Patterson2015; @Edge2017]. Phased variants give haplotypes; '
 'phased reads, tagged with their haplotype, are themselves the input of haplotype-resolved analyses, from '
 'allele-specific methylation and expression to the local assembly and structural variant calling of each '
 'haplotype. Accurate long reads such as PacBio HiFi [@Wenger2019] make this possible on a single sample.'),
('p',
 'Read-based phasers differ in how they solve the problem. WhatsHap minimises the number of read alleles that must '
 'be corrected for the reads to split into two haplotypes [@Patterson2015; @Martin2016], HapCUT2 maximises a '
 'likelihood by iterated cuts [@Edge2017], LongPhase finds the two heaviest disjoint paths through a graph of '
 'variants linked by reads, and phases structural variants with them [@Lin2022], and HiPhase searches with A* and '
 'phases small, structural and tandem repeat variants jointly [@Holt2024]. Each phases a set of variants a caller '
 'found on a linear reference, typically DeepVariant [@Poplin2018; @Shafin2021]; longcallD calls and phases '
 'variants together, resolving noisy regions with local multiple sequence alignment [@Gao2026].'),
('p',
 'Phasing from a linear reference has two blind spots that matter for the reads more than for the variants. '
 'Where the caller finds no heterozygous variant, in tandem repeats and other regions where read alignments '
 'disagree, the reads have nothing to be phased by and are left unphased. Where a region is duplicated in the '
 "sample but not in the reference, the reads of both copies align to one place, and the copies' differences look "
 'like heterozygous variants: they co-segregate as perfectly as real ones, and the reads they label, and the '
 'blocks they join, can be phased in the wrong orientation by the thousand. Phasers are usually compared by their '
 'phased variants, through switch errors and block lengths [@Lin2022; @Holt2024], which count neither the reads '
 'left unphased nor a block of reads joined wrongly by a few variants.'),
('p',
 'A pangenome graph describes the known alleles of many genomes at once [@Liao2023; @Hickey2024], and graph '
 'aligners place reads on it [@Garrison2018; @Siren2021]. Where a read crosses a site of the graph, a snarl '
 '[@Paten2018], the walk it takes names the allele it carries: an observation of a known allele, made without '
 'calling it, including alleles a linear alignment smears across a repeat. Known alleles have been genotyped this '
 'way [@Siren2021; @Ebler2022]. They are not enough to phase a sample, though, since every genome carries '
 'heterozygous variants that no other genome of the graph does.'),
('p', 'pgphase phases reads with both. Our contributions are:'),
('list',
 ['**Phasing at known alleles from graph alignments** (Sections 2.2–2.3): a catalog of the sites of a pangenome '
  "graph, and each read's allele at them read from its graph alignment.",
  "**One solve over the catalog's sites and the sample's own** (Sections 2.4–2.5): heterozygous variants found in "
  "the linear alignments, including those a multiple sequence alignment resolves in noisy regions, join the "
  "catalog's sites before phasing, and an expectation-maximization learns every site's error rate, shrunk toward "
  "that of sites of its kind.",
  '**Joins that the reads verify** (Section 2.6): a block is joined to the next only when several reads carry '
  'evidence on both sides at reliable sites; sites at collapsed duplications, sites read only at low mapping '
  'quality and one-base length differences cannot join blocks, though they still label reads.',
  '**Read-level evaluation** (Section 2.8): every read labelled with its parent of origin, chromosomes held out '
  "from development, and each switch error classed by the tool's haplotype sequences against the truth's."]),
('p',
 'On eight chromosomes of HG002, pgphase phased more reads correctly, and left fewer misphased or unphased, than '
 'HiPhase, LongPhase, WhatsHap and longcallD on every chromosome, held-out ones included, with the fewest switch '
 'errors ([#tab:reads]; [#tab:switches]).'),

('sec', 'Materials and methods'),
('sub', 'Overview'),
('p',
 "pgphase takes the reads' alignments to a pangenome graph, as a GAF file [@Li2020], the graph's site catalog, and "
 'optionally the same reads aligned to the linear reference, as a BAM file ([#fig:overview]). It works on chunks '
 'of the reference (1 Mb by default with a BAM), in parallel, and joins the chunks through the reads they share. '
 'Its two parts are the *graph arm*, which phases the catalog\'s sites from the graph alignments, and *union gap '
 "phasing*, which adds the sample's own heterozygous variants from the linear alignments and phases everything in "
 'one solve. Its outputs are a phased VCF and the reads\' haplotype (HP) and phase set (PS) tags.'),
('figure', 'overview'),
('sub', 'Site catalog and read evidence from a pangenome graph'),
('p',
 "The catalog is built once per graph (`build-snarl-catalog`, from a GBZ graph [@Siren2022]): each snarl "
 '[@Paten2018] of the graph is a site, its reference allele the walk of the reference path through it and its '
 'alternative alleles the distinct walks of the haplotype paths, ordered by how many haplotypes take them. The '
 "reads' graph alignments are indexed once by reference position (`index-gaf`). For each chunk, every "
 "read crossing a site contributes the allele whose walk its own alignment follows; the identity of an allele is "
 'its walk, so two reads carry the same allele exactly when they take the same path through the site. A '
 'multi-allelic site is split into one reference-against-alternative contrast per alternative allele, and sites '
 'are kept by depth and allele fraction. Indels in homopolymers, tandem repeats and low-complexity sequence are '
 'set apart as noisy candidates, which do not take part in the first solve.'),
('sub', 'Phasing the catalog sites'),
('p',
 "The catalog's clean heterozygous sites are phased first by k-means clustering of the reads on their alleles, "
 "the algorithm of longcallD [@Gao2026] applied to walk observations. With a BAM, each graph SNP is checked against "
 'the bases of well-mapped reads: one whose reference base no such read shows, and no deletion explains, is '
 'reclassified homozygous. Where a read takes a walk that differs from a selected allele only away from a SNP, the '
 'walk is read as that SNP\'s allele, so that reads whose alignment takes a rarer path through a site still vote. '
 'The clean solve\'s labels start the joint solve below.'),
('sub', "The sample's own heterozygous variants"),
('p',
 'With a BAM, each chunk is also solved from its linear alignments, by the alignment pipeline pgphase ports from '
 'longcallD [@Gao2026]: candidate variants from the reads\' differences to the reference, classified by depth, '
 'allele fraction and sequence context, with noisy regions up to 2 kb resolved by aligning each haplotype\'s reads '
 'with partial order alignment [@Gao2021]. Its clean heterozygous variants, and those verified by the multiple '
 'sequence alignment, join the chunk\'s sites with their per-read alleles. A verified noisy variant must be seen '
 'by two reads of mapping quality 30 on each of its two alleles, and no injected site may be read by more than '
 '1.5 times the chunk\'s median depth: at a collapsed duplication the two "alleles" are the copies, and a run of '
 'such sites agrees with itself perfectly. Variants that fail these rules still label reads, but only as islands '
 '(Section 2.7). The same indel placed at several copies of a repeat becomes one site; calls inside a catalog site '
 'of up to 2 kb only give alleles to reads the graph left without one there; reads only the linear alignment '
 'holds are added. At injected indels, every read crossing the site without a call is realigned to the two '
 'alleles in reference context and takes the closer one, with a per-read error of '
 '\\(1/(1+e^{3m})\\) for a margin of \\(m\\) edits. At catalog SNPs a read\'s graph allele is replaced by its '
 'aligned base where the two disagree: measured against the parents, the alignment\'s base was right in about 90% '
 'of such conflicts.'),
('sub', 'One expectation-maximization over all sites'),
('p',
 'Every site with two observed alleles, from the catalog or the sample, gets a phase \\(h_s\\in\\{0,1\\}\\) and '
 'an error rate \\(\\varepsilon_s\\), and every read \\(r\\) a posterior probability of each haplotype. With '
 '\\(x_{rs}\\) the allele read \\(r\\) shows at site \\(s\\) and \\(e_{rs}\\) the per-read error of that call where '
 'one is known, a read on haplotype \\(a\\) has likelihood'),
('p',
 '\\[ L_r(a) = \\prod_{s} \\bigl(1-\\varepsilon_{rs}\\bigr)^{[x_{rs} = h_s \\oplus a]}\\, \\varepsilon_{rs}^{[x_{rs} \\neq h_s \\oplus a]}, \\]'),
('p',
 'where \\(\\varepsilon_{rs} = 1-(1-\\varepsilon_s)(1-e_{rs})\\) combines the site\'s error with the call\'s, '
 'and expectation-maximization [@Dempster1977] alternates the reads\' posteriors with the sites\' phases and '
 'errors. A site\'s error is the weighted disagreement of its reads, shrunk toward the median error of the sites '
 'of its class in the chunk (clean catalog heterozygotes, other catalog sites, injected sites and local windows, '
 'each split by whether the site is a one-base length difference) by five pseudo-reads, an empirical Bayes prior '
 'that keeps a noisy repeat call from overruling a clean SNP that shares its reads. Only reads of mapping quality '
 '20 or more teach the model; the others are scored against it and labelled. A site whose rarer allele is held by '
 'fewer than a quarter of its reads is left out: a skewed split is the signature of a paralog or an error. After '
 'the first iterations, the solve tries switch moves, flipping every site past the boundary whose flip most '
 'increases the likelihood. A site whose learned error is 0.3 or more is reported unphased and labels no read; a '
 'read is labelled in the block where its own evidence is strongest when its posterior there is at least 0.8.'),
('sub', 'Joining blocks only where reads verify the join'),
('p',
 'A block ends wherever flipping everything after a boundary would cost less than 4 log-likelihood units, or '
 'fewer than two reads carry evidence on both sides: one molecule can be chimeric or misplaced. Beyond that, '
 'pgphase restricts which sites can carry a join at all, because a wrong join flips every read past it:'),
('list',
 ['**One-base length differences are weak.** A homopolymer length call in a homozygous stretch looks like a '
  'heterozygous site; a join that holds only through such sites needs a margin of 10 log units instead of 4.',
  '**Strong joins need reads.** A join counts as strong only when at least four reads carry evidence at '
  'reliable, non-weak sites on both sides of the boundary.',
  '**Collapsed duplications and poorly mapped sites carry no joins.** A site called by more than twice the '
  "chunk's median number of reads, or by more than 1.3 times with its rarer allele under 40% of them, or whose "
  "reads' mean mapping quality is under 40, labels reads but joins nothing.",
  '**One locus, one vote.** A local window and a call inside it, or several records of one event within 30 bp '
  'whose reads split the same way, are counted once.',
  '**No hinges.** A block that hangs on a single site with weak joins on both sides, crossed by no read, is cut '
  'at the weaker join.']),
('p',
 'These rules leave more, shorter blocks where the reads cannot verify a join, and the reads inside each block '
 'phased correctly.'),
('sub', 'Labels for the remaining reads, and output'),
('p',
 'Reads the solve leaves unlabelled are labelled from evidence the solve does not vote with. Local haplotype '
 'windows (up to 2 kb) contrast the two haplotype sequences of a noisy locus, each built as the medoid of reads '
 'labelled on either side; allele windows contrast the two exact alleles of a heterozygous indel the solve does '
 'not observe, including two-allele sites with no reference read; and length windows genotype tandem repeat '
 'tracts under a stutter model learned from the chunk. A second pass of the solve uses the first two kinds as '
 'sites. A read whose calls in a block all agree takes the block\'s label from posterior 0.7, and a read still '
 'unlabelled is phased by its own call at one site that cannot vote, in a block of that site alone, an *island* '
 'that never joins anything, so that a wrong island costs only its own reads. Chunks are joined by the reads '
 'they share, only when one orientation leads by at least three reads and 90%; a read without a label takes one '
 'from the phase set its alleles at excluded sites associate with (one-sided binomial test, \\(p \\le 0.01\\)). '
 'In the VCF, a phase set whose heterozygotes all lie at one locus is written unphased, as are SNPs of learned '
 'error 0.1 or more and indels of 0.25 or more.'),
('sub', 'Evaluation'),
('p',
 '**Data.** HG002 [@Zook2016] PacBio HiFi reads (DRAFT: source and chemistry) were aligned with vg Giraffe '
 'v1.73.0 [@Siren2021] (HiFi preset, minimizers \\(k=31\\), \\(w=50\\)) to a 16-haplotype sample [@Siren2024] '
 'of the Human Pangenome Reference Consortium release 2 Minigraph-Cactus graph on CHM13 [@Liao2023; @Hickey2024; '
 '@HPRC2], from which HG002 and its parents had been removed (229 assembled samples, CHM13 and GRCh38), and '
 'surjected onto CHM13 [@Nurk2022]; they covered chr19 at 59-fold. The catalog was built from the whole graph. Every tool was given the same reads: pgphase their graph '
 'alignments and the surjected BAM, the other tools the BAM. HiPhase 1.6.0, LongPhase 2.0.2 and WhatsHap 2.8, by '
 'default and with `--distrust-genotypes` (WhatsHap-dg), phased DeepVariant 1.10.0 calls (PACBIO model) made on '
 'the same BAM; longcallD (v0.0.11) called and phased its own. pgphase was developed on chr12, chr17, chr18, chr19 '
 'and chr20, whose truth guided its rules and constants; chr6, chr16 and chr21, and the MHC of chr6, were held out '
 'until the comparison.'),
('p',
 '**Truth.** Each read\'s parent of origin was found by aligning it with minimap2 [@Li2018] to both haplotypes of '
 'the HG002 v1.1 assembly [@Hansen2025] and taking the better alignment with HipHap [@HipHap], which also scores '
 'the confidence of the assignment (HapQ). Truth heterozygous variants are the Q100 HG002 small-variant benchmark '
 'on CHM13 [@Hansen2025].'),
('p',
 '**Read metrics.** A primary read is *phaseable* when its alignment spans a heterozygous truth variant. Each '
 'phase set takes the parental orientation of the majority of its truth-labelled reads; a read is *correct* when '
 'its tag agrees with it, *discordant* when it does not, and *unphased* when it has no tag. Reads are on the '
 '*arms* or in the centromeric and pericentromeric satellites of CHM13\'s cenSat annotation [@Altemose2022] (for '
 'chr21, its acrocentric short arm too). We also report the discordant reads whose parental label is confident '
 '(HapQ ≥ 30): below HapQ 10, every tool is 25–37% discordant, and most such reads carry the other parent\'s bases '
 'at their own clean SNVs, so their label rather than their phase is wrong.'),
('p',
 '**Variant metrics.** A tool\'s phased heterozygous biallelic records are matched to the truth by position, '
 'reference and alternative allele. Within a phase set, a single matched heterozygote whose orientation '
 'disagrees with its agreeing neighbours is a *flip* and a lasting change a *switch*. In repeats, the same '
 'haplotypes can be written as different records, and a record-level scorer then counts switches that are not. '
 'We therefore class each switch by sequence: around each of the two heterozygotes bounding it, a window of 75 bp '
 'on either side, grown until no record of either VCF crosses its edges (at most 20 kb), the tool\'s two '
 'haplotype sequences are compared with the truth\'s by edit distance, and each window pairs them in the cheaper '
 'way. A switch is *real* when the two windows pair them in opposite ways, *representation* when both pair them '
 'the same way, and *undecided* when either window\'s two pairings are within 30% (or 2 edits) of each other.'),

('sec', 'Results'),
('sub', 'More reads phased correctly, fewer misphased'),
('p',
 f'Over the eight chromosomes, {V["phaseable"]} reads on the arms span a heterozygous truth variant. pgphase '
 f'phased {V["pg_correct_pct"]}% of them correctly, {V["pg_disc_pct"]}% discordantly and left '
 f'{V["pg_unph_pct"]}% unphased ([#tab:reads]). The best of the other tools, {V["best_rival"]}, phased '
 f'{V["best_rival_correct_pct"]}% correctly: pgphase phased {V["more_correct_than_best"]} more reads correctly, '
 f'and left {V["err_ratio"]} times fewer reads discordant or unphased than the tool with the fewest such reads, '
 f'{V["fewest_err_rival"]}. It phased the most reads correctly on every chromosome, and on '
 f'{8 - V["n_disc_lost"]} of the eight left the fewest discordant reads; {V["disc_lost"]}, while leaving many '
 f'more reads unphased ([#fig:accuracy]). Among reads whose parental label is '
 f'confident, pgphase phased {V["pg_conf_disc"]} discordantly, against {V["low_conf_rival_disc"]} for '
 f'{V["low_conf_rival"]}, the lowest of the others.'),
('table', 'reads'),
('figure', 'accuracy'),
('p',
 f'The margin held on chromosomes the method was not developed on. On the held-out chr6, chr16 and chr21, pgphase '
 f'phased {V["held_margin_min"]} to {V["held_margin_max"]} percentage points more phaseable reads correctly than '
 f'the best other tool on each, against {V["dev_margin_min"]} to {V["dev_margin_max"]} points on the development '
 f'chromosomes. In the MHC, whose HLA genes are among the most polymorphic of the genome, it phased '
 f'{V["mhc_pg_pct"]}% correctly, against {V["mhc_best_rival_pct"]}% for {V["mhc_best_rival"]}.'),
('sub', 'Switch errors'),
('p',
 f'pgphase made {V["pg_switches"]} switch errors on the arms of the eight chromosomes, against '
 f'{V["fewest_sw_rival_n"]} for {V["fewest_sw_rival"]}, the fewest among the others, and {V["hiphase_switches"]} '
 f'for HiPhase ([#tab:switches]). Classed by haplotype sequence, {V["pg_real"]} of pgphase\'s switches were real, '
 f'{V["pg_repr"]} wrote the truth\'s haplotypes with different records, and {V["pg_undec"]} were undecided, '
 f'against {V["rival_real_min"]} to {V["rival_real_max"]} real switches for the other tools. Undecided switches '
 'are mostly tandem repeats where neither of a tool\'s haplotypes matches either of the truth\'s, so that no '
 f'pairing is better than the other. pgphase made {V["pg_flips"]} flips, against {V["fewest_flip_rival_n"]} for '
 f'{V["fewest_flip_rival"]}.'),
('table', 'switches'),
('sub', 'Phased variants and blocks'),
('p',
 f'pgphase phased {V["pg_hets"]} truth heterozygous variants, against {V["most_hets_rival_n"]} for '
 f'{V["most_hets_rival"]}, the most of the others ([#tab:switches]). Its blocks are shorter than HiPhase\'s: '
 f'{V["pg_ps"]} phase sets against {V["hiphase_ps"]}, a median block N50 of {V["pg_n50_kb"]} kb against '
 f'{V["hiphase_n50_kb"]} kb. The difference is where the join rules refuse to join: across collapsed '
 'duplications and stretches read only at low mapping quality, where each site is its own block. The reads '
 'there are phased within their blocks, and a join that the reads cannot verify is exactly what flips them.'),
('sub', 'Centromeres'),
('p',
 f'In centromeric and pericentromeric satellites, where reads map ambiguously and copies differ by few bases, '
 f'pgphase phased {V["pg_cen_correct"]} reads correctly at {V["pg_cen_ratio"]} correct reads per discordant one '
 f'([#tab:cen]). {V["most_cen_rival"]} phased {V["most_cen_rival_correct"]} correctly at '
 f'{V["most_cen_rival_ratio"]} per discordant read, phasing more of them but with more of them wrong. '
 f'{V["cen_best_ratio_rival"]} reached {V["cen_best_ratio"]} correct reads per discordant one by phasing few: '
 f'{V["cen_best_ratio_correct"]} correctly, with {V["cen_best_ratio_unph"]} left unphased.'),
('table', 'cen'),
('sub', 'Running time and memory'),
('p',
 f'On the {V["n_same_threads"]} chromosomes where every tool ran with the same threads, pgphase took '
 f'{V["pg_wall_total_min"]} minutes with six, against {V["hiphase_wall_total_min"]} for HiPhase and '
 f'{V["longphase_wall_total_min"]} for LongPhase; on chr6, {V["pg_wall_chr6"]} s against {V["hiphase_wall_chr6"]} s '
 f'and {V["longphase_wall_chr6"]} s ([#tab:resources]). '
 'Its time includes finding the sample\'s variants and writing the haplotagged reads, while the others\' excludes '
 f'DeepVariant. It needed up to {V["pg_rss_max"]} GB of memory, more than the phasers given DeepVariant\'s calls.'),
('table', 'resources'),

('sec', 'Discussion'),
('p',
 'A phaser is judged by its phased variants when its output is a phased VCF, but much of the use of long-read '
 'phasing is in the reads: a read left unphased is lost to every haplotype-resolved analysis, and a block of '
 'reads joined in the wrong orientation is worse. pgphase gains reads in two ways. It has more sites to phase '
 'them by: the graph\'s known alleles, read from where each read\'s alignment goes rather than called, and the '
 'sample\'s own heterozygous variants, so that fewer reads lack a site. And it is slow to join: the sites that '
 'make phasers misphase blocks of reads, collapsed duplications, poorly mapped stretches and homopolymer length '
 'calls, still label the reads that cross them but cannot hold two blocks together. The cost is more and shorter '
 'blocks than HiPhase\'s, where the reads cannot tell how blocks connect.'),
('p',
 'pgphase has limits. Its catalog holds only the alleles of the genomes in the graph: at a tandem repeat whose '
 'expansion in the sample no genome shares, a read\'s graph alignment takes the closest known walk, and the '
 'most of the switch errors we could not class are of this kind. It needs graph alignments of '
 'the reads, whose cost we do not count here, and more memory than the phasers it was compared with. We evaluated '
 'one sample and one technology; ONT reads, other samples and other graphs remain to be tested, and the '
 'development chromosomes shaped its rules even though the held-out ones were phased as well. Finally, a '
 'record-level comparison counts as switches differences of representation in repeats; classing switches by '
 'haplotype sequence, as here, or comparing haplotypes rather than records, would make phasers easier to compare.'),
]

BACK = [
 ("Author contributions", "DRAFT: to be completed."),
 ("Supplementary data", "Supplementary data are available at *Bioinformatics* online."),
 ("Conflict of interest", "None declared."),
 ("Funding", "DRAFT: to be completed."),
 ("Data availability",
  "pgphase is available at https://github.com/kokyriakidis/pgphase. DRAFT: the evaluation scripts, the exact "
  "commands of every tool (Supplementary Section S4) and the scores behind every table are to be released with "
  f"the paper; the results reported here are of pgphase {V['commit']}."),
]

FIGURES = {
 "overview": "pgphase's pipeline. The graph arm phases the sites of a pangenome graph's catalog from the walks the "
             "reads' graph alignments take through them; with a linear alignment, union gap phasing adds the "
             "sample's own heterozygous variants, and one expectation-maximization phases all sites, joining blocks "
             "only where reads verify the join. Reads still unlabelled are phased from local windows, rescue and "
             "island blocks.",
 "accuracy": "Discordant against unphased phaseable reads on the arms, per tool and chromosome (eight points per "
             "tool; both axes logarithmic). Lower left is better.",
}

WIDE_TABLES = ("reads", "switches")
