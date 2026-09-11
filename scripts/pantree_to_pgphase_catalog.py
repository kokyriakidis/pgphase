#!/usr/bin/env python3
"""Convert a pantree reference-tree variant catalog into a pgphase sites VCF.

pantree (Nowbandegani et al., Cell Genomics 2026; https://github.com/oclb/pantree)
defines a variant as a graph edge that is *not* in a spanning "reference tree" of
the pangenome graph.  That finds variants `vg deconstruct` cannot emit -- notably
SNPs whose REF and ALT are both non-GRCh38, which concentrate in segmental
duplications (+46.1% SNPs in SDs vs. the superbubble catalog).  Those are exactly
the paralog-distinguishing anchors pgphase lacks in the regions where its residual
switch errors live.

pantree's own VCF is not consumable by `collect-graph-variation`: it carries no
`AT` (allele traversal) field, which is the only thing pgphase matches reads
against (see src/graph_sites.cpp, src/graph_query.cpp).  This script loads the
graph once through pantree's Python API and emits a drop-in pgphase site catalog:

  * `AT=<ref_walk>,<alt_walk>` synthesized from the reference-tree path and the
    variant-edge path between the *same* pair of flanking handles, which is the
    boundary invariant `graph_site_validation_skip_reason()` enforces.
  * `ID=>u>v` -- the oriented variant edge, which pgphase uses as the site key.
  * `POS`/`REF`/`ALT`/`VT`/`AC`/`AN` computed by pantree's own writer internals,
    so this catalog is directly comparable to plain `pantree gfa2vcf` output.
  * No nesting fields.  Variant edges are independent by construction, so there
    is no `LV`/`PS`/`PA` parent gating -- pgphase treats every site as top level.

Non-reference variants (pantree `REF='.'`, real allele in `NR`) are emitted as
`REF=N` / `ALT=*` with `PTNONREF=1`.  They stay usable as phasing anchors -- the
allele walks are intact -- but `collect-graph-variation` skips `*` alleles when
deriving variant keys (src/graph_collect.cpp), so they never reach the output VCF.
Use `--nonref drop` for a catalog that is byte-comparable with a deconstruct one.

Usage:
    python3 pantree_to_pgphase_catalog.py graph.gfa \
        -o chr20.pantree.sites.vcf.gz \
        --chrom 'CHM13#0#chr20' --ref-name CHM13

Requires: pantree importable (`uv pip install -e /path/to/pantree`, or run under
`uv run` from a pantree checkout), plus `bgzip` and `tabix` on PATH.

The GFA must be the same graph pgphase phases against:
    vg chunk --gbz --contig chr20 -x full.gbz -o /tmp/chunk_chr20
    vg convert -f /tmp/chunk_chr20_graph_0_chr20.gbz > chr20.gfa

NOTE ON COST: pantree is a research-grade Python implementation.  The paper
reports 36 GB / 10 h for chr1.  Budget accordingly; this is one-time preprocessing
per graph version, like `build-snarl-catalog`.
"""

import argparse
import os
import subprocess
import sys
import tempfile

# pantree INFO tags that are re-emitted verbatim.  Descriptions are copied from
# pantree's writer so the header self-documents without importing private state.
PANTREE_INFO_HEADERS = {
    "NR": '##INFO=<ID=NR,Number=1,Type=String,Description="Non-reference allele (pantree)">',
    "VT": '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type (pantree)">',
    "TP": '##INFO=<ID=TP,Number=1,Type=Integer,Description="Tree position of the variant edge branch point (pantree)">',
    "AC": '##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">',
    "AN": '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    "HP": '##INFO=<ID=HP,Number=.,Type=String,Description="Haplotype positions at reference tree edge (pantree)">',
    "TR_MOTIF": '##INFO=<ID=TR_MOTIF,Number=1,Type=String,Description="Tandem repeat motif (pantree)">',
    "NIA": '##INFO=<ID=NIA,Number=1,Type=Integer,Description="Nearly identical alleles (pantree)">',
    "UIDX": '##INFO=<ID=UIDX,Number=1,Type=Integer,Description="Index of node u (pantree)">',
}

# pantree's RC is "REF allele count"; vg deconstruct's RC is "reference contig".
# Renamed so a pantree catalog can be concatenated with a deconstruct one without
# a conflicting header declaration or a silently mistyped field.
RENAMED_INFO = {"RC": "PTRC"}

EXTRA_INFO_HEADERS = [
    '##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency (AC/AN)">',
    '##INFO=<ID=PTRC,Number=1,Type=Integer,Description="REF allele count (pantree RC, renamed to avoid clashing with vg deconstruct RC)">',
    '##INFO=<ID=PTNONREF,Number=0,Type=Flag,Description="Variant whose REF allele is absent from the linear reference; REF/ALT masked to N/* so it anchors phasing but is not emitted as a call">',
]

VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def import_pantree():
    """Import pantree, including the private writer internals we reuse.

    Reusing `_VariantData` / `_VariantRecord` rather than reimplementing the
    allele and position rules keeps POS/REF/ALT identical to `pantree gfa2vcf`,
    which is what makes the differential comparison meaningful.
    """
    try:
        from pantree.graph import PangenomeGraph
        from pantree.dfs import dfs_methods
        from pantree.utils import node_recover, edge_complement
        from pantree.vcf import (
            _VariantData,
            _VariantRecord,
            _get_default_info_fields,
        )
    except ImportError as exc:
        sys.exit(
            f"error: cannot import pantree ({exc}).\n"
            "  Install it into this environment, e.g.:\n"
            "    git clone https://github.com/oclb/pantree.git && cd pantree && uv sync\n"
            "  then run this script with that environment's python, or `uv run`.\n"
            "  If the import failure names _VariantData/_VariantRecord, pantree's\n"
            "  private writer API has changed and this script needs updating."
        )
    return {
        "PangenomeGraph": PangenomeGraph,
        "dfs_methods": dfs_methods,
        "node_recover": node_recover,
        "edge_complement": edge_complement,
        "_VariantData": _VariantData,
        "_VariantRecord": _VariantRecord,
        "_get_default_info_fields": _get_default_info_fields,
    }


def allele_walk_nodes(graph, first, last, variant_edges):
    """Oriented node list from `first` to `last`, inclusive of both endpoints.

    `walk_with_variants` returns *interior* nodes only: it drops `first`, and
    drops `last` except when the walk flips strand (+ to -), where the '-' node's
    start is its complement's end and so it is already included.  pgphase needs
    both flanking handles present and identical across alleles, so re-add them.
    """
    interior = list(graph.walk_with_variants(first, last, variant_edges))
    nodes = [first] + interior
    last_already_included = graph.direction(first) == 1 and graph.direction(last) == -1
    if not last_already_included:
        nodes.append(last)
    return nodes


def build_allele_traversals(graph, edge, node_recover):
    """Return (ref_at, alt_at) as pgphase walk strings, or None if not derivable.

    The endpoints mirror `PangenomeGraph.ref_alt_alleles`: the walk runs from the
    branch point to `v` on the forward strand, or from `u` to the branch point on
    the reverse strand.  Both alleles therefore share `first` and `last`, which is
    the `incompatible_allele_boundaries` invariant in graph_sites.cpp.
    """
    u, v = edge
    branch_point = graph.edges[edge]["branch_point"]
    if graph.direction(u) == 1:
        first, last = branch_point, v
    else:
        first, last = u, branch_point

    ref_nodes = allele_walk_nodes(graph, first, last, [])
    alt_nodes = allele_walk_nodes(graph, first, last, [edge])
    if len(ref_nodes) < 2 or len(alt_nodes) < 2:
        return None

    ref_at = "".join(node_recover(n) for n in ref_nodes)
    alt_at = "".join(node_recover(n) for n in alt_nodes)
    if ref_at == alt_at:
        return None
    return ref_at, alt_at


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        out[key] = value
    return out


def render_info(pairs):
    parts = []
    for key, value in pairs:
        parts.append(key if value is None else f"{key}={value}")
    return ";".join(parts)


def passes_filters(info, ref_allele, alt_allele, args):
    """Catalog-level pruning.  Returns (ok, reason)."""
    vt = info.get("VT", ".")
    if args.snps_only and vt != "SNP":
        return False, "not_snp"
    if args.exclude_vt and vt in args.exclude_vt:
        return False, f"vt_{vt}"

    if args.max_allele_len is not None:
        longest = max(len(ref_allele), len(alt_allele))
        if longest > args.max_allele_len:
            return False, "allele_too_long"

    # AC/AN are only meaningful when genotypes were computed.  Frequency gating
    # is the main reason to keep genotypes on: a catalog trimmed to panel-common
    # sites is both smaller and more likely to carry a real het in the sample.
    ac_raw = info.get("AC", ".")
    an_raw = info.get("AN", ".")
    try:
        ac = int(ac_raw.split(",")[0])
    except ValueError:
        ac = None
    try:
        an = int(an_raw)
    except ValueError:
        an = None

    if args.min_ac is not None:
        if ac is None:
            return False, "no_ac"
        if ac < args.min_ac:
            return False, "ac_below_min"
    if (args.min_af is not None or args.max_af is not None):
        if ac is None or not an:
            return False, "no_af"
        af = ac / an
        if args.min_af is not None and af < args.min_af:
            return False, "af_below_min"
        if args.max_af is not None and af > args.max_af:
            return False, "af_above_max"
    return True, ""


def build_header(args, sources):
    lines = ["##fileformat=VCFv4.2"]
    lines.append(f"##source={sources}")
    lines.append(f"##pantree_dfs_method={args.dfs_method}")
    if args.priority_samples:
        lines.append(f"##pantree_priority_samples={args.priority_samples}")
    lines.append(f"##pantree_ref_name={args.ref_name}")
    lines.extend(EXTRA_INFO_HEADERS)
    lines.extend(PANTREE_INFO_HEADERS.values())
    if args.contig_length:
        lines.append(f"##contig=<ID={args.chrom},length={args.contig_length}>")
    else:
        lines.append(f"##contig=<ID={args.chrom}>")
    lines.append(VCF_COLUMNS)
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("gfa", help="Pangenome graph in GFA (vg convert -f from a GBZ chunk)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output sites VCF; .vcf.gz is bgzipped and tabix-indexed")
    parser.add_argument("--chrom", required=True,
                        help="CHROM value, matching the deconstruct catalog (e.g. 'CHM13#0#chr20')")
    parser.add_argument("--ref-name", default="CHM13",
                        help="Reference sample name in the graph [CHM13]")
    parser.add_argument("--contig-length", type=int, default=None,
                        help="Contig length for the ##contig header (optional)")
    parser.add_argument("--dfs-method", default="max_weight",
                        choices=["max_weight", "contiguous"],
                        help="Reference-tree DFS strategy [max_weight]")
    parser.add_argument("--priority-samples", default=None,
                        help="Comma-separated sample names prioritized by --dfs-method contiguous")
    parser.add_argument("--no-genotypes", action="store_true",
                        help="Skip genotype computation (faster; disables AC/AN/AF and all frequency filters)")
    parser.add_argument("--no-missingness", action="store_true",
                        help="Skip missingness computation (recommended for vg-built graphs)")
    parser.add_argument("--nonref", choices=["keep", "drop"], default="keep",
                        help="Variants whose REF is absent from the linear reference: "
                             "keep as N/* anchors, or drop [keep]")
    parser.add_argument("--snps-only", action="store_true",
                        help="Emit only VT=SNP sites")
    parser.add_argument("--exclude-vt", default=None,
                        help="Comma-separated VT values to drop (e.g. INV,DUP)")
    parser.add_argument("--max-allele-len", type=int, default=None,
                        help="Drop sites whose REF or ALT exceeds this length")
    parser.add_argument("--min-ac", type=int, default=None,
                        help="Minimum panel alt allele count")
    parser.add_argument("--min-af", type=float, default=None,
                        help="Minimum panel alt allele frequency (AC/AN)")
    parser.add_argument("--max-af", type=float, default=None,
                        help="Maximum panel alt allele frequency (AC/AN)")
    parser.add_argument("--stats", default=None,
                        help="Write a per-reason drop-count TSV here")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Log pantree progress to stderr")
    args = parser.parse_args()

    args.exclude_vt = set(args.exclude_vt.split(",")) if args.exclude_vt else set()
    if args.no_genotypes and (args.min_ac is not None or args.min_af is not None
                              or args.max_af is not None):
        sys.exit("error: --no-genotypes leaves AC/AN undefined; drop the frequency filters "
                 "or drop --no-genotypes")

    pt = import_pantree()
    PangenomeGraph = pt["PangenomeGraph"]
    node_recover = pt["node_recover"]
    edge_complement = pt["edge_complement"]
    _VariantData = pt["_VariantData"]
    _VariantRecord = pt["_VariantRecord"]

    logger = None
    if args.verbose:
        import logging
        logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                            format="[%(asctime)s] %(message)s")
        logger = logging.getLogger("pantree")

    priority_dict = None
    if args.priority_samples:
        names = [s.strip() for s in args.priority_samples.split(",")]
        priority_dict = {name: i for i, name in enumerate(names)}

    print(f"[pantree->pgphase] loading {args.gfa}", file=sys.stderr)
    graph = PangenomeGraph.from_gfa(
        args.gfa,
        ref_name=args.ref_name,
        logger=logger,
        dfs_method=pt["dfs_methods"][args.dfs_method],
        priority_dict=priority_dict,
    )

    if args.no_genotypes:
        sample_to_genotype = {}
    else:
        print("[pantree->pgphase] computing genotypes", file=sys.stderr)
        sample_to_genotype = graph.genotypes_from_gfa(
            args.gfa, True, skip_missing=args.no_missingness
        )

    info_fields = pt["_get_default_info_fields"]()
    reference_edges = graph.get_reference_edges()

    drops = {}
    def drop(reason):
        drops[reason] = drops.get(reason, 0) + 1

    n_written = 0
    n_nonref = 0

    # Stream to a temp file with an explicit sort prefix rather than buffering
    # every record: on a whole chromosome this is millions of lines on top of a
    # graph that already dominates RAM.  Tie-break on UIDX so ties resolve in
    # graph order, matching pantree's own sort.
    tmp_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    with tempfile.NamedTemporaryFile("w", dir=tmp_dir, prefix=".pantree_sites.",
                                     suffix=".unsorted", delete=False) as raw:
        raw_path = raw.name
        print("[pantree->pgphase] emitting variant edges", file=sys.stderr)
        for u, v in graph.sorted_variant_edges(exclude_terminus=True):
            if graph.direction(u) == -1 and graph.direction(v) == -1:
                u, v = edge_complement((u, v))
            edge = (u, v)

            try:
                reference_edge = reference_edges[edge]
            except KeyError:
                drop("no_reference_edge")
                continue

            try:
                variant_info = _VariantData.from_graph(
                    graph=graph,
                    edge=edge,
                    reference_edge=reference_edge,
                    chr_name=args.chrom,
                    sample_to_genotype=sample_to_genotype,
                )
                record = _VariantRecord.from_variant_data(
                    variant_info=variant_info,
                    edge=edge,
                    genotype_records=[],
                    info_fields=info_fields,
                )
            except Exception as exc:  # noqa: BLE001 - research code, keep going
                drop(f"pantree_record_error:{type(exc).__name__}")
                continue

            # A variant with no reference position cannot be retrieved by the
            # tabix region query that drives collect-graph-variation.
            if record.vcf_position is None or record.vcf_position < 1:
                drop("no_reference_position")
                continue

            walks = build_allele_traversals(graph, edge, node_recover)
            if walks is None:
                drop("no_allele_traversal")
                continue
            ref_at, alt_at = walks

            info = parse_info(record.info)
            ref_allele = record.ref_allele
            alt_allele = record.alt_allele

            is_nonref = ref_allele == "."
            if is_nonref:
                if args.nonref == "drop":
                    drop("nonref")
                    continue
                # Keep the site as a phasing anchor but make it un-emittable:
                # graph_collect.cpp skips ALT='*' when deriving variant keys, so
                # these never reach the phased VCF while their walks still carry
                # read-to-haplotype signal.  The real allele stays in NR.
                ref_allele, alt_allele = "N", "*"

            ok, reason = passes_filters(info, record.ref_allele, record.alt_allele, args)
            if not ok:
                drop(reason)
                continue

            uidx = info.get("UIDX", "-1")
            ordered = []
            for key, value in info.items():
                ordered.append((RENAMED_INFO.get(key, key), value))
            ordered.append(("AT", f"{ref_at},{alt_at}"))
            try:
                ac = int(info.get("AC", ".").split(",")[0])
                an = int(info.get("AN", "."))
                if an:
                    ordered.append(("AF", f"{ac / an:.6g}"))
            except ValueError:
                pass
            if is_nonref:
                ordered.append(("PTNONREF", None))
                n_nonref += 1

            line = "\t".join([
                args.chrom,
                str(record.vcf_position),
                record.variant_id,
                ref_allele,
                alt_allele,
                record.qual,
                record.filter_field,
                render_info(ordered),
            ])
            raw.write(f"{args.chrom}\t{record.vcf_position}\t{uidx}\t{line}\n")
            n_written += 1

    header = build_header(
        args,
        sources=f"pantree_to_pgphase_catalog.py (pantree reference-tree catalog, dfs={args.dfs_method})",
    )

    plain_out = args.output[:-3] if args.output.endswith(".gz") else args.output
    print(f"[pantree->pgphase] sorting {n_written} records", file=sys.stderr)
    with open(plain_out, "w") as out:
        out.write(header)
        out.flush()
        sort = subprocess.Popen(
            ["sort", "-k1,1", "-k2,2n", "-k3,3n", raw_path],
            stdout=subprocess.PIPE, env={**os.environ, "LC_ALL": "C"},
        )
        cut = subprocess.Popen(["cut", "-f4-"], stdin=sort.stdout, stdout=out)
        sort.stdout.close()
        cut.communicate()
        sort.wait()
        if sort.returncode or cut.returncode:
            sys.exit("error: sorting the catalog failed")
    os.unlink(raw_path)

    if args.output.endswith(".gz"):
        subprocess.run(["bgzip", "-f", plain_out], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", args.output], check=True)

    total_dropped = sum(drops.values())
    print(f"[pantree->pgphase] wrote {n_written} sites to {args.output} "
          f"({n_nonref} non-reference anchors, {total_dropped} dropped)", file=sys.stderr)
    for reason, count in sorted(drops.items(), key=lambda kv: -kv[1]):
        print(f"  dropped {count:>10}  {reason}", file=sys.stderr)

    if args.stats:
        with open(args.stats, "w") as fh:
            fh.write("metric\tcount\n")
            fh.write(f"written\t{n_written}\n")
            fh.write(f"nonref_anchors\t{n_nonref}\n")
            for reason, count in sorted(drops.items(), key=lambda kv: -kv[1]):
                fh.write(f"dropped_{reason}\t{count}\n")


if __name__ == "__main__":
    main()
