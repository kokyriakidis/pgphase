"""The main paper, written once; render.py turns it into paper/ (LaTeX) and build/paper/paper.html.

Inline markup: **bold**, *italic*, `code`, [@Key; @Key] parenthetical citation, [@!Key] narrative citation,
\\( \\) inline maths, \\[ \\] display maths, [#tab:x] / [#fig:x] cross-references.
Blocks: ("sec", title), ("sub", title), ("p", text), ("list", [items]) numbered, ("table", key), ("figure", key).
Every number comes from tables.V (results.json), never typed here; the claims the text makes about every
chromosome are asserted against the same data, so the paper does not render once the data stops supporting them.
"""
from tables import TABLES, V
import supplement as S


def ssec(title):
    """A supplementary section's number, by its title."""
    return "S%d" % ([b[1] for b in S.BODY if b[0] == "sec"].index(title) + 1)


def stab(key):
    """A supplementary table's label."""
    return S.TABLES[key]["label"]

assert V["fewest_err_every_chrom"] and V["fewest_switch_every_chrom"] and V["rows_best_correct"] == 8

TITLE = "pgphase: pangenome graph-aware phasing of long reads with known and sample-specific variants"
# Authors in order, each with an email; both are corresponding authors (as in the dinara-align paper).
AUTHORS = [("Konstantinos Kyriakidis", "kkyriaki@ucsc.edu"), ("Benedict Paten", "bpaten@ucsc.edu")]
AFFILIATION = "UC Santa Cruz Genomics Institute, Santa Cruz, CA, USA"
KEYWORDS = ["phasing", "haplotagging", "long reads", "pangenome graph", "graph alignment", "expectation maximization"]

ABSTRACT = [
 ("Motivation",
  "Long reads can be assigned to the parental haplotype they were sequenced from, but current phasers do so "
  "through heterozygous variants called against a single linear reference. Reads in regions where the caller "
  "finds no heterozygous variant stay unphased, and in collapsed duplications the differences between paralogous "
  "copies mimic heterozygous variants and cause reads, and whole blocks, to be phased in the wrong orientation. Haplotype-"
  "resolved analyses that use the reads themselves need as many reads as possible phased, and phased correctly."),
 ("Results",
  "We present pgphase, which phases reads at the known alleles of a pangenome graph, observed directly in the reads' "
  "graph alignments, together with the sample's own heterozygous variants found in its linear alignments, so that "
  "no separate variant caller is needed. A single "
  "expectation-maximization learns how reliable each site is, and phase blocks are joined only where several reads "
  "support the join. On HG002 PacBio HiFi reads, with each read's parent of origin as truth on eight chromosomes, "
  f"three of them not used to tune the method, pgphase phased {V['pg_correct_pct1']}% of {V['phaseable_m']} "
  f"million phaseable reads correctly, compared with {V['best_rival_correct_pct1']}% for the best of HiPhase, "
  f"LongPhase, WhatsHap and longcallD. It left {V['err_ratio']} times fewer reads misphased or unphased than the "
  f"best alternative, and made {V['pg_switches']} switch errors against {V['fewest_sw_rival_n']} for the next best "
  f"tool. In centromeric satellites it phased the most reads correctly and made {V['pg_cen_switches']} switch errors "
  f"against {V['rival_cen_switches_min']} to {V['rival_cen_switches_max']}."),
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
 'haplotype. Long, accurate reads such as PacBio HiFi [@Wenger2019] make both possible from a single sample, '
 'without parental data.'),
('p',
 'Read-based phasers differ in how they solve the problem. WhatsHap minimises the number of read alleles that must '
 'be corrected for the reads to split into two haplotypes [@Patterson2015; @Martin2016], HapCUT2 maximises a '
 'likelihood by iterated cuts [@Edge2017], LongPhase finds the two heaviest disjoint paths through a graph of '
 'variants linked by reads, and phases structural variants with them [@Lin2022], and HiPhase searches with A* and '
 'phases small, structural and tandem repeat variants jointly [@Holt2024]. Each phases a set of variants a caller '
 'found on a linear reference, typically DeepVariant [@Poplin2018; @Shafin2021]; longcallD calls and phases '
 'variants together, resolving noisy regions with local multiple sequence alignment [@Gao2026].'),
('p',
 'Phasing through variants called on a linear reference has two weaknesses that affect reads more than variants. '
 'Where the caller finds no heterozygous variant, as in many tandem repeats and other regions where read '
 'alignments disagree, the reads there carry nothing to be phased by and stay unphased. Where a region is '
 "duplicated in the sample but not in the reference, the reads of both copies align to one place, and the copies' "
 'differences look like heterozygous variants. They are as consistent across the reads as real heterozygous '
 'variants, so a phaser can join two blocks through them in the wrong orientation and misphase every read beyond '
 'the join. Phasers are usually compared by their phased variants, through switch '
 'errors and block lengths [@Lin2022; @Holt2024]. These metrics count a wrong join as one switch however many '
 'reads it misphases, and do not count the reads left unphased.'),
('p',
 'A pangenome graph represents the known alleles of many genomes at once [@Liao2023; @Hickey2024], and graph '
 'aligners place reads on it [@Garrison2018; @Siren2021; @Chang2025]. Where a read crosses a site of the graph, a snarl '
 '[@Paten2018], the walk its alignment takes names the allele it carries, so the allele is observed without being '
 'called, including in repeats where alignments to a single reference place the same variant inconsistently. '
 'Pangenome references have been used to genotype known variants [@Siren2021; @Ebler2022]. Known alleles alone '
 'cannot phase a sample, however, because every genome carries rare and private heterozygous variants that the '
 'graph does not contain.'),
('p', 'pgphase phases reads with both kinds of variant. Our contributions are:'),
('list',
 ['**Phasing at known alleles from graph alignments** (Sections 2.2–2.3): a catalog of the sites of a pangenome '
  "graph, and each read's allele at them taken from its graph alignment.",
  "**One solve over the catalog's sites and the sample's own** (Sections 2.4–2.5): heterozygous variants from the "
  "linear alignments, including those resolved in noisy regions, join the catalog's sites, and an "
  "expectation-maximization learns every site's error rate. pgphase therefore needs no separate variant "
  "caller: the phasers we compare it with, except longcallD, phase the calls of DeepVariant run first.",
  '**Joins that the reads verify** (Section 2.6): two blocks are joined only when several reads carry evidence on '
  'both sides at reliable sites; sites at collapsed duplications, sites read only at low mapping quality and '
  'one-base length differences cannot join blocks, although they still label reads.',
  '**Read-level evaluation** (Section 2.8): the parent of origin of every read as truth, chromosomes not used for '
  "tuning, and each switch error classified by comparing the tool's haplotype sequences with the truth's."]),
('p',
 'On each of eight chromosomes of HG002, including the three not used for tuning, pgphase phased more reads '
 'correctly, left fewer reads misphased or unphased, and made fewer switch errors than HiPhase, LongPhase, WhatsHap '
 'and longcallD ([#tab:reads]; [#tab:switches]). It did so in centromeric satellites too, which the graph '
 'represents, phasing the most reads there correctly with a fraction of the other tools\' switch errors. Neither '
 'the graph alignments nor the linear alignments alone come close: together they phase 99% of reads correctly, '
 'each alone about 90% ([#tab:ablation]).'),

('sec', 'Materials and methods'),
('sub', 'Overview'),
('p',
 "pgphase takes the reads' alignments to a pangenome graph, as a GAF file [@Li2020], the graph's site catalog and, "
 'optionally, the same reads aligned to the linear reference, as a BAM file ([#fig:overview]). It processes the '
 'reference in chunks (1 Mb by default when a BAM is given) in parallel and joins adjacent chunks through the '
 "reads they share. It has two parts: the *graph arm*, which phases the catalog's sites from the graph "
 "alignments, and *union gap phasing*, which adds the sample's own heterozygous variants from the linear "
 'alignments and phases all sites in one solve. It writes a phased VCF and the reads\' haplotype (HP) and phase set '
 '(PS) tags. All thresholds below were set on the development chromosomes (Section 2.8); Supplementary '
 f'Table {stab("params")} lists every one with its value.'),
('figure', 'overview'),
('sub', 'Site catalog and read evidence from a pangenome graph'),
('p',
 "The catalog is built once per graph, from a GBZ graph [@Siren2022] with `build-snarl-catalog`: each snarl "
 '[@Paten2018] of the graph is a site, whose reference allele is the walk of the reference path through it and '
 'whose alternative alleles are the distinct walks of the haplotype paths, ordered by the number of haplotypes '
 "that take them. The reads' graph alignments are indexed once by reference position (`index-gaf`). For each "
 'chunk, every read with mapping quality 5 or more (1 or more when a BAM is given) that crosses a site contributes '
 'the allele whose walk its alignment follows. An allele is identified by its walk, so two reads carry the same '
 'allele exactly when they take the same path through the site. A multi-allelic site is split into one '
 'reference-against-alternative contrast per alternative allele, and sites are filtered by depth and allele '
 'fraction. Indels in homopolymers, tandem repeats and low-complexity sequence are classified as noisy and left '
 'out of the first solve.'),
('sub', 'Phasing the catalog sites'),
('p',
 "The catalog's clean heterozygous sites are phased first by k-means clustering of the reads on their alleles, "
 "the algorithm of longcallD [@Gao2026] applied to walk observations. With a BAM, each catalog SNP is tested "
 'against the aligned bases of reads of mapping quality 30 or more: a SNP whose reference base is absent from '
 'these reads (a binomial test, Bonferroni-corrected at 1% over the chunk\'s SNPs) and whose absence no deletion '
 'explains is reclassified as homozygous. Where the two alleles of a SNP site differ at a single node of the graph, '
 "a read whose walk follows another allele of the snarl takes the SNP allele whose node that walk contains, when "
 'it contains exactly one of them; reads on rarer walks through the site thus still contribute. The labels of this '
 'clean solve initialise the joint solve below.'),
('sub', "The sample's own heterozygous variants"),
('p',
 'With a BAM, each chunk is also solved from its linear alignments by the alignment pipeline pgphase ports from '
 'longcallD [@Gao2026]: candidate variants are collected from the reads\' differences to the reference and '
 'classified by depth, allele fraction and sequence context, and noisy regions of up to 2 kb are resolved by '
 "aligning each haplotype's reads with partial order alignment [@Gao2021]. Its clean heterozygous variants, and the "
 'noisy-region variants the multiple sequence alignment verifies, join the chunk\'s sites with their per-read '
 'alleles. A verified noisy variant must be observed by at least two reads of mapping quality 30 or more on each of '
 'its two alleles, and no site observed by more than 1.5 times the chunk\'s median depth is added: at a collapsed '
 'duplication the two "alleles" are the copies, and a run of such sites is perfectly consistent with itself. '
 'Verified noisy variants that fail these rules are still used, but only to label reads as islands (Section 2.7). '
 'A clean call inside a noisy region longer than the MSA limit cannot be verified; it is kept as an *unverified* '
 'call, which takes no part in joins and labels reads only as a last resort (Section 2.7). '
 'Placements of one indel at different copies of a repeat are merged into one site; a call inside a catalog site '
 'of up to 2 kb only supplies alleles to reads the graph left without one there, and one inside a longer site is '
 'added as a site of its own; and reads with mapping quality 20 or more that only the linear alignment holds are '
 'added. At each added indel, every read that crosses the site '
 'without a call is realigned to the two alleles in reference context and takes the closer one, with a per-read '
 'error of \\(1/(1+e^{3m})\\) for a margin of \\(m\\) edits; where the MSA called the read, its call is kept, '
 'with error 0.25 if realignment cannot decide and 0.35 if it calls the other allele. Where a read\'s allele at a '
 'catalog SNP disagrees with its aligned base, the base is used: measured against the parents on the development '
 'chromosomes, the aligned base was correct in 88–93% of such disagreements. A graph allele at a SNP whose base '
 'the read\'s alignment does not give keeps an error of 0.25.'),
('sub', 'One expectation-maximization over all sites'),
('p',
 'Every site with two observed alleles, from the catalog or the sample, is given a phase \\(h_s\\in\\{0,1\\}\\) '
 'and an error rate \\(\\varepsilon_s\\), and every read \\(r\\) a posterior probability for each haplotype. With '
 '\\(x_{rs}\\in\\{0,1\\}\\) the allele read \\(r\\) shows at site \\(s\\) and \\(e_{rs}\\) the per-read error of '
 'that call (zero when none is known), a read from haplotype \\(a\\) has likelihood'),
('p',
 '\\[ L_r(a) = \\prod_{s} \\bigl(1-\\varepsilon_{rs}\\bigr)^{[x_{rs} = h_s \\oplus a]}\\, \\varepsilon_{rs}^{[x_{rs} \\neq h_s \\oplus a]}, \\]'),
('p',
 'where \\(\\varepsilon_{rs} = 1-(1-\\varepsilon_s)(1-e_{rs})\\) combines the error of the site with that of '
 'the call. Expectation-maximization [@Dempster1977] alternates between the reads\' posteriors and the sites\' '
 'phases and errors, starting from error 0.05 at clean sites and 0.2 at the others and from the clean solve\'s '
 'labels; site errors are kept between 0.01 and 0.45. The error of a site is the weighted disagreement of its '
 'reads, shrunk by five pseudo-reads toward the median error of the sites of its class in the chunk that carry '
 'at least ten reads\' teaching weight (the weighted mean when fewer than five such sites exist) (clean catalog heterozygotes, other catalog sites, added '
 'variants, local windows and unverified calls from noisy regions longer than 2 kb, each divided by whether the '
 'site is a one-base length difference). This empirical Bayes prior keeps, for example, a noisy homopolymer call '
 'from overruling a clean SNP that shares its reads. Only reads of mapping quality 20 or more update the model; '
 'the others are scored against it and labelled. A read whose calls all lie at one site does not teach that '
 'site\'s error, since its posterior comes from the site alone, except at clean catalog heterozygotes, whose '
 'starting error is the better prior. A site whose rarer allele is held by fewer than a quarter of its '
 'reads is left out, since such a skewed split is typical of a paralogous copy or of an error. After 15 '
 'iterations, the solve flips all sites beyond the boundary whose flip most increases the likelihood, for as long '
 'as some flip increases it (at most 30 rounds). A site whose learned error is 0.3 or more is reported unphased '
 'and labels no read, and a read is labelled in the block where its evidence is strongest if its posterior there '
 'is at least 0.8.'),
('sub', 'Joining blocks only where reads verify the join'),
('p',
 'A block ends wherever flipping all sites beyond a boundary would cost less than 4 log-likelihood units, or where '
 'fewer than two reads carry evidence on both sides, since a single molecule can be chimeric or misplaced. The '
 'cut is judged over all reads, including those below mapping quality 20: a read too ambiguously mapped to shape '
 'the model still shows that its molecule continues across the boundary. Because '
 'a wrong join misphases every read beyond it, pgphase further restricts which sites can support a join:'),
('list',
 ['**One-base length differences are weak.** A homopolymer length call in a homozygous stretch can look '
  'heterozygous. Where the join margin from the other sites alone is below 4, a join that relies on such sites '
  'needs a margin of 10 log-likelihood units.',
  '**Strong joins need several reads.** A margin counts as strong only when at least four reads carry evidence at '
  'reliable, non-weak sites on both sides of the boundary.',
  '**Collapsed duplications and poorly mapped sites support no joins.** A site called by more than twice the '
  "chunk's median number of reads, by more than 1.3 times that number with its rarer allele below 40%, or by reads "
  'whose mean mapping quality is below 40, labels reads but supports no join.',
  '**Redundant evidence counts once.** In the join test, a read\'s call at a site inside a local window is ignored '
  'when the read also has a call at the window. A site within 30 bp of another that shares at least 80% of its '
  'reads, 95% of which split the same way at both, repeats that site\'s evidence and does not vote on read labels.',
  '**No hinges.** A site that alone connects two blocks through weak joins on both sides, with no read spanning '
  'it, is cut at the weaker join.']),
('sub', 'Labels for the remaining reads, and output'),
('p',
 'Reads that the solve leaves unlabelled are labelled from evidence the solve does not vote with. Local haplotype '
 'windows, up to 2 kb long, contrast the two haplotype sequences of a noisy locus, each the medoid of reads '
 'labelled on one side; allele windows contrast the two exact alleles of a heterozygous indel the solve does not '
 'observe, including two-allele sites with no reference read; and length windows genotype tandem repeat tracts '
 'under a stutter model learned from the chunk. A second pass of the solve uses haplotype and allele windows as '
 'sites. After the solve, adjacent chunks are joined in the orientation supported by the reads they share, only '
 'when one orientation leads by at least three reads and 90%. Sites left out of the solve are then oriented from '
 'the labelled reads when their alleles associate with one phase set (exact one-sided binomial test, '
 '\\(p \\le 0.01\\)), and label the unlabelled reads that cross them, using a read\'s linear-alignment allele '
 'where its graph alignment gives none. Last, a read whose calls in a block all '
 "agree takes the block's label if its posterior there is at least 0.7, and a read still unlabelled is phased by its own call at one "
 'site that cannot vote, in a block of that site alone (an *island*), which joins nothing, so that a wrong island '
 'affects only its own reads. Island sites are the sites the solve rejected, repeat indels it never sees, verified '
 'noisy variants that failed the rules of Section 2.4, and length windows; a site must be called by two '
 'unlabelled reads and have both alleles in at least a quarter of its reads, and a read takes the one with the '
 'lowest error. Before islands, a '
 'read with only unverified calls takes their label at posterior 0.8 when their learned error is below 0.1. In the '
 'VCF, multi-base substitutions are written as their '
 'single-base changes, indels are left-aligned, and two records of one variant are collapsed, keeping the phased '
 'one; a phase set whose heterozygotes all lie at one locus (within 10 bp) is '
 'written unphased, as are SNPs of learned error 0.1 or more, indels of 0.25 or more and unverified calls.'),
('sub', 'Evaluation'),
('p',
 '**Data.** HG002 [@Zook2016] PacBio HiFi reads (DRAFT: source and chemistry) were aligned with vg Giraffe '
 'v1.73.0 [@Siren2021; @Chang2025], in its HiFi preset with minimizers \\(k=31\\) and \\(w=50\\), to a '
 '16-haplotype sample [@Siren2024] '
 'of the Human Pangenome Reference Consortium release 2 Minigraph-Cactus graph on CHM13 [@Liao2023; @Hickey2024; '
 '@HPRC2], from which HG002 and its parents had been removed (229 assembled samples, CHM13 and GRCh38), and the '
 'alignments were surjected onto CHM13 [@Nurk2022]; they covered chr19 at 59-fold. The site catalog was built from '
 'the full graph rather than from the 16-haplotype sample. All tools received the same reads: pgphase their graph alignments and the surjected BAM, the other '
 'tools the BAM. HiPhase 1.6.0, LongPhase 2.0.2 and WhatsHap 2.8, run by default and with `--distrust-genotypes` '
 '(WhatsHap-dg), phased DeepVariant 1.10.0 calls (PACBIO model) made on the same BAM; longcallD v0.0.11 called '
 'and phased its own. pgphase was developed on chr12, chr17, chr18, chr19 and chr20, whose truth guided its rules '
 'and thresholds. chr6, chr16 and chr21 were not used to choose thresholds; they were scored during development to '
 'check that changes generalised, and one rule (the hinge rule, Section 2.6) was found by examining an error on '
 'chr6. The MHC (chr6:28.5–34.0 Mb) is reported separately.'),
('p',
 '**Truth.** The reads were aligned separately to the maternal and the paternal haplotype of the HG002 v1.1 '
 'assembly [@Hansen2025] with minimap2 2.31 [@Li2018] (preset `lr:hqae`, no secondary alignments). HipHap '
 '[@HipHap] compared the two alignments of each read, assigned the read to the haplotype with the better '
 'alignment score and scored the assignment\'s confidence (HapQ), using minimap2\'s match score rather than an '
 'estimate from a random sample of reads, so that the labels are reproducible. Reads without a parental label are not scored. Truth heterozygous variants are those of the Q100 '
 'HG002 small-variant benchmark on CHM13 [@Hansen2025].'),
('p',
 '**Read metrics.** A primary read is *phaseable* when its alignment spans a heterozygous truth variant. Each '
 'phase set takes the parental orientation of the majority of its truth-labelled reads; a read is *correct* when '
 'its haplotype tag agrees with that orientation, *discordant* when it does not, and *unphased* when it has no '
 'tag. A read is centromeric when its alignment midpoint lies in one interval per chromosome spanning its '
 'centromeric and pericentromeric satellites, as annotated on CHM13 [@Altemose2022] and listed in '
 f'Supplementary Table {stab("scoring")}, and on the arms '
 'otherwise. The graph represents the satellites: its catalog held '
 f'{V["alpha_density_min"]} to {V["alpha_density_max"]} sites per Mb in the alpha-satellite arrays of these '
 f'chromosomes, against {V["arm_density_min"]} to {V["arm_density_max"]} on their arms, and '
 f'{V["alpha_gaf_min"]} to {V["alpha_gaf_max"]} graph alignments per chromosome lay in the arrays. We report the '
 'satellites apart from the arms because reads map there ambiguously and the parental assignments are least '
 'reliable. One region at the end of chr20 where the truth labels contradict the truth variants is not scored '
 f'(Supplementary Table {stab("scoring")}). We also report the discordant reads whose parental assignment is '
 'confident (HapQ ≥ 30): on the development chromosomes, every tool was 25–37% discordant on reads below HapQ 10, '
 'most of which carry the other parent\'s alleles at their own clean SNVs, which points to an error in the label '
 'rather than in the phase. All tools were scored by the same script.'),
('p',
 '**Variant metrics.** A tool\'s phased heterozygous biallelic records are matched to the truth by position, '
 'reference and alternative allele. Within a phase set, a single matched heterozygote whose orientation '
 'disagrees with agreeing neighbours on both sides is a *flip*, and a lasting change of orientation is a '
 '*switch*. Switches and flips are counted separately on the arms and in the centromeric intervals. In repeats, '
 'the same haplotypes can be written as different records, and a record-level comparison '
 'then counts switches where the haplotypes agree. We therefore classify each switch by sequence: around each of '
 'the two heterozygotes bounding it, a window of 75 bp on either side, extended until no record of either VCF '
 "crosses its edges (at most 20 kb), the tool's two haplotype sequences are compared with the truth's by edit "
 'distance, and each window pairs them in the cheaper way. A switch is *real* when the two windows pair the '
 'haplotypes in opposite ways, *representation* when both pair them the same way, and *undecided* when the two '
 "pairings of either window differ by less than 30% (or by fewer than 2 edits). A phase set spans from its first "
 "to its last phased heterozygous record; we report the number of phase sets and the N50 (not NG50) of their "
 "spans. Running time and peak memory were measured "
 "with GNU time on one machine with an Intel Core i9-7900X processor (10 cores) and 91 GB of memory."),

('sub', 'Implementation'),
('p',
 'pgphase is written in C++17 on htslib, with abPOA [@Gao2021] for partial order alignment, WFA2-lib and edlib '
 'for pairwise alignment; its alignment pipeline is derived from longcallD [@Gao2026], whose MIT-licensed code it '
 'ports and tests against. Catalog building and GAF indexing are done by pgcat, a Rust program that reads the '
 'GBZ graph with the gbz library [@Siren2022]. Chunks are processed by a pool of worker threads, each with its own '
 'file handles, and joined in one sequential pass; output does not depend on the number of threads.'),

('sec', 'Results'),
('sub', 'More reads phased correctly, fewer misphased'),
('p',
 f'Over the eight chromosomes, {V["phaseable"]} reads on the arms span a heterozygous truth variant. pgphase '
 f'phased {V["pg_correct_pct"]}% of them correctly, {V["pg_disc_pct"]}% discordantly, and left '
 f'{V["pg_unph_pct"]}% unphased ([#tab:reads]). The best of the other tools, {V["best_rival"]}, phased '
 f'{V["best_rival_correct_pct"]}% correctly; pgphase phased {V["more_correct_than_best"]} more reads correctly, '
 f'and left {V["err_ratio"]} times fewer reads discordant or unphased than the tool with the fewest such reads, '
 f'{V["fewest_err_rival"]}. On each chromosome, it left {V["fold_min"]} to {V["fold_max"]} times fewer reads '
 f'misphased or unphased than the best other tool on that chromosome ([#fig:accuracy]b). It phased the most '
 f'reads correctly on every chromosome and left the fewest discordant '
 f'reads on {V["n_disc_best_word"]} of the eight; {V["disc_lost"]}, while leaving many more reads unphased '
 f'([#fig:accuracy]). Among reads whose parental assignment is confident, pgphase phased {V["pg_conf_disc"]} '
 f'discordantly, against {V["low_conf_rival_disc"]} for {V["low_conf_rival"]}, the fewest among the others.'),
('table', 'reads'),
('figure', 'accuracy'),
('p',
 f'The margin was similar on the chromosomes not used for tuning: on chr6, chr16 and chr21, pgphase phased '
 f'{V["held_margin_min"]} to {V["held_margin_max"]} percentage points more phaseable reads correctly than the best '
 f'other tool on each, against {V["dev_margin_min"]} to {V["dev_margin_max"]} points on the development '
 f'chromosomes. In the MHC (chr6:28.5–34.0 Mb), which holds the HLA genes, among the most polymorphic of the '
 f'genome, pgphase phased {V["mhc_pg_pct"]}% of {V["mhc_phaseable"]} phaseable reads correctly, against '
 f'{V["mhc_best_rival_pct"]}% for {V["mhc_best_rival"]}, and left {V["mhc_pg_err"]} misphased or unphased, '
 f'{V["mhc_err_ratio"]} times fewer than {V["mhc_next_rival"]} ({V["mhc_next_err"]}). {V["mhc_fewer_disc_tools"]} '
 f'misphased fewer reads there ({V["mhc_fewer_disc_min"]} to {V["mhc_fewer_disc_max"]}, against '
 f'{V["mhc_pg_disc"]}) but left at least {V["mhc_unph_min_rival"]} unphased, against {V["mhc_unph_pg"]}.'),
('sub', 'Graph and linear alignments complement each other'),
('p',
 f'Run on its graph alignments alone, pgphase phased {V["abl_graphonly_correct"]}% of phaseable reads correctly and '
 f'left {V["abl_graphonly_unph"]}% unphased; run on its linear alignments alone, through its port of longcallD\'s '
 f'pipeline, {V["abl_bamonly_correct"]}% and {V["abl_bamonly_unph"]}%, close to longcallD itself '
 f'({V["abl_lcd_correct"]}% and {V["abl_lcd_unph"]}%) ([#tab:ablation]). Together the two phased '
 f'{V["pg_correct_pct"]}% correctly and left {V["pg_unph_pct"]}% unphased: each kind of evidence phases reads the '
 'other leaves without a site, the graph\'s known alleles where linear alignments disagree and the sample\'s own '
 f'variants where the graph has none. The combination also made fewer arm switch errors, {V["pg_switches"]}, than '
 f'either part alone ({V["abl_graphonly_switches"]} and {V["abl_bamonly_switches"]}), as the joint solve sees more '
 'sites at every boundary.'),
('table', 'ablation'),
('sub', 'Switch errors'),
('p',
 f'pgphase made {V["pg_switches"]} switch errors on the arms of the eight chromosomes, fewer than any other tool on '
 f'each chromosome; {V["fewest_sw_rival"]} made the next fewest, {V["fewest_sw_rival_n"]}, and HiPhase '
 f'{V["hiphase_switches"]} ([#tab:switches]); that is {V["pg_switch_rate"]} switches per 10,000 phased truth '
 f'heterozygotes on the arms, against {V["rival_rate_min"]} to {V["rival_rate_max"]} for the other tools. A record-level count includes switches that are not phasing errors, '
 'so we checked each one against the haplotype sequences (Section 2.8). A switch is *real* when the tool\'s '
 'haplotypes around it pair with the truth\'s one way on one side and the other way on the other: the phase '
 'really changes there. It is a difference of *representation* when the haplotypes pair the same way on both '
 'sides, so the tool wrote the same haplotypes as the truth with different records, and *undecided* when the '
 'sequences cannot tell, because the tool\'s haplotypes differ from both of the truth\'s by similar amounts. '
 f'{V["pg_real"]} of pgphase\'s arm switches were real ({V["pg_real_chr19"]} of them on chr19), {V["pg_repr"]} '
 f'were differences of representation and {V["pg_undec"]} undecided, against {V["rival_real_min"]} to '
 f'{V["rival_real_max"]} real switches for the other tools. pgphase\'s undecided switches lie mostly in tandem '
 "repeats whose expansion in the sample no genome of the graph carries (Section 4). In the centromeric and "
 f'pericentromeric satellites, pgphase made {V["pg_cen_switches"]} switch errors, against '
 f'{V["rival_cen_switches_min"]} to {V["rival_cen_switches_max"]} for the other tools. On the arms it made '
 f'{V["pg_flips"]} flips, against {V["fewest_flip_rival_n"]} for {V["fewest_flip_rival"]}.'),
('table', 'switches'),
('sub', 'Phased variants and blocks'),
('p',
 f'{V["most_hets_rival"]}, which calls its own variants, phased the most truth heterozygous variants, '
 f'{V["most_hets_rival_n"]}; pgphase phased {V["pg_hets"]}, more than any of the phasers given DeepVariant\'s calls '
 f'({V["best_dv_hets_n"]} for {V["best_dv_hets_rival"]}; [#tab:switches]). Its blocks are shorter than HiPhase\'s: '
 f'{V["pg_ps"]} phase sets against {V["hiphase_ps"]}, and a median block N50 of {V["pg_n50_kb"]} kb against '
 f'{V["hiphase_n50_kb"]} kb. Many of the additional phase sets lie in centromeric satellites, which pgphase '
 f'phases and HiPhase largely leaves unphased ([#tab:cen]), and in collapsed duplications and stretches read only at low '
 'mapping quality, where the join rules leave each site in its own block.'),
('sub', 'Centromeres'),
('p',
 f'In centromeric and pericentromeric satellites, where reads map ambiguously and copies differ by few bases, '
 f'pgphase phased more reads correctly than any other tool, {V["pg_cen_correct"]}, with {V["pg_cen_ratio"]} correct '
 f'reads per discordant one, and made the fewest switch errors, {V["pg_cen_switches"]} against '
 f'{V["rival_cen_switches_min"]} to {V["rival_cen_switches_max"]} ([#tab:cen]; [#tab:switches]). {V["most_cen_rival"]} phased {V["most_cen_rival_correct"]} correctly with '
 f'{V["most_cen_rival_ratio"]} correct reads per discordant one. {V["cen_best_ratio_rival"]} had the highest ratio, '
 f'{V["cen_best_ratio"]}, but phased only {V["cen_best_ratio_correct"]} reads correctly and left '
 f'{V["cen_best_ratio_unph"]} unphased.'),
('table', 'cen'),
('sub', 'Running time and memory'),
('p',
 f'On the {V["n_same_threads"]} chromosomes where every multi-threaded tool ran with six threads, pgphase took '
 f'{V["pg_wall_total_min"]} minutes in total, against {V["hiphase_wall_total_min"]} for HiPhase and '
 f'{V["longphase_wall_total_min"]} for LongPhase; on chr6, {V["pg_wall_chr6"]} s against {V["hiphase_wall_chr6"]} s '
 f'and {V["longphase_wall_chr6"]} s (Supplementary Table {stab("resources")}). '
 'pgphase\'s time includes finding the sample\'s own variants and writing the haplotagged reads. HiPhase, '
 f'LongPhase and WhatsHap first need DeepVariant\'s calls, which took {V["dv_total_min"]} minutes on the same '
 f'chromosomes and threads: with them, HiPhase took {V["hiphase_with_calls_min"]} minutes and LongPhase '
 f'{V["longphase_with_calls_min"]}, {V["calls_speedup"]} times or more pgphase\'s time. No time includes read '
 f'alignment, which every tool needs. pgphase needed up to {V["pg_rss_max"]} GB of memory, more than the phasers '
 'given DeepVariant\'s calls.'),

('sec', 'Discussion'),
('p',
 'Phasers are usually compared through their phased variants, but many uses of long-read phasing depend on the '
 'reads: a read left unphased is lost to every haplotype-resolved analysis, and a block joined in the wrong '
 'orientation misassigns all of its reads. On these terms pgphase phased more reads correctly than the tools we '
 'compared it with, on every chromosome we tested, centromeric satellites included. Two design choices explain '
 'this. First, it has more sites to phase reads with: the graph\'s known alleles, observed in the reads\' graph '
 'alignments rather than called, and the sample\'s own heterozygous variants. Each source alone phased about 90% '
 'of phaseable reads correctly, no better than the other tools; together they phased 99% ([#tab:ablation]). '
 'Second, it is conservative in joining blocks: sites at collapsed duplications, poorly mapped sites and '
 'homopolymer length calls still label the reads that cross them but cannot join two blocks. The cost of the '
 'second is more and shorter blocks than HiPhase produces where the reads cannot verify how blocks connect. '
 'pgphase also needs no separate variant caller, which simplifies a pipeline and shortens it: the graph\'s '
 'known alleles are observed in the reads\' graph alignments and the sample\'s own heterozygous variants are '
 f'found by its alignment solve, whereas DeepVariant followed by HiPhase or LongPhase took {V["calls_speedup"]} '
 'times or more as long here. We '
 'did not compare its variant calls with DeepVariant\'s; the evaluation here concerns phasing.'),
('p',
 'Our evaluation has limits. All tools were given reads aligned to a pangenome graph and surjected onto CHM13, '
 'not the linear alignments (for example from minimap2 or pbmm2) the other tools are usually run on, and their '
 'accuracy with such alignments may differ. We evaluated one sample and one sequencing technology on eight '
 'chromosomes; other samples, ONT reads and other graphs remain to be tested. The development chromosomes shaped '
 "pgphase's rules and thresholds, and although the chromosomes not used for tuning showed a similar margin, they "
 'were scored during development. The read truth relies on alignment to the parental assemblies, and its labels '
 'are least reliable where the parental haplotypes are most similar.'),
('p',
 'pgphase also has limits of its own. Its catalog holds only the alleles of the genomes in the graph: at a tandem '
 "repeat whose expansion in the sample no genome shares, a read's graph alignment follows the closest known walk, "
 'and most of the switches we could not classify are of this kind. Representing such alleles from the reads '
 'themselves is a natural next step. pgphase also needs graph alignments of the reads and more memory than the '
 'phasers it was compared with. Finally, a record-level comparison counts differences of representation in '
 'repeats as switch errors; classifying switches by haplotype sequence, as here, or comparing haplotypes rather '
 'than records, would make phasers easier to compare.'),
]

BACK = [
 ("Author contributions",
  "Konstantinos Kyriakidis (Conceptualization [lead], Methodology [lead], Software [lead], Validation [lead], "
  "Formal analysis [lead], Investigation [lead], Data curation [lead], Visualization [lead], Writing—original "
  "draft [lead]), Benedict Paten (Supervision [lead], Writing—review & editing [equal])."),
 ("Supplementary data", "Supplementary data are available at *Bioinformatics* online."),
 ("Conflict of interest", "None declared."),
 ("Funding", "This work received no specific funding."),
 ("Data availability",
  f"pgphase is available at https://github.com/kokyriakidis/pgphase; the results reported here are of pgphase "
  f"commit {V['commit']}. The scores behind every table and figure are in the same repository (`tools/paper/results.json`), "
  "from which `tools/paper/` regenerates this paper. The evaluation scripts, which run every tool and score its "
  "output, are available at https://github.com/kokyriakidis/pgphase-dev, and the exact command of every tool is "
  f"given in Supplementary Section {ssec('Commands')}. The HG002 reads (DRAFT: source), the HPRC release 2 graph [@HPRC2], the HG002 "
  "v1.1 assembly and the Q100 benchmark [@Hansen2025] are publicly available."),
]

FIGURES = {
 "overview": "pgphase's pipeline. The reads are aligned by other tools (grey, dashed boxes); pgphase's steps lie inside the "
             "dashed outline. The graph arm (blue) builds the graph's site catalog and phases its sites from the walks "
             "the reads' graph alignments take through them; union gap phasing (orange), which runs only when a "
             "linear alignment is given, adds the sample's own heterozygous variants. One expectation-maximization "
             "phases all sites, joining blocks only where reads verify the join; local windows and a second pass "
             "follow, then chunk stitching, the rescue and island blocks label the remaining reads, and pgphase "
             "writes a phased VCF and the reads' haplotype tags.",
 "accuracy": "Phaseable reads on the chromosome arms that each tool misphases or leaves unphased. (a) Over the eight "
             "chromosomes, the share misphased (discordant) and the share left unphased, one row per tool; tools are "
             "ordered by the sum of the two. (b) Per chromosome, the share misphased or unphased: pgphase (vermillion) "
             "joined to the best other tool on that chromosome (dark), the other tools faint, and at the right how many "
             "times fewer reads pgphase leaves misphased or unphased. chr6, chr16 and chr21 were not used for tuning.",
}

WIDE_TABLES = ("reads", "switches")
