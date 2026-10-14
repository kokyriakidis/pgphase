"""A line-for-line port of HiPhase's allele assignment (evaluation only).

Source: PacificBiosciences/HiPhase, src/wfa_graph.rs
(WFAGraph::from_reference_variants_with_hom, edit_distance_with_pruning),
src/read_parsing.rs (global_realignment, local_realignment, the per-read
fallback in load_full_read_segments), src/data_types/variants.rs
(closest_allele_clip), src/block_gen.rs (filter_out_alignment_record,
get_variant_type), src/phaser.rs (variant loading, reference buffer) and
src/cli.rs defaults: min MAPQ 5, reference buffer 15, global max edit
distance 500, max edit distance ratio 0.1, pruning distance 500.

Alleles: 0 = the variant's allele0 (REF, or the lower ALT at a 1/2 site),
1 = allele1, 2 = ambiguous, 3 = no overlap (AlleleType order).
"""
import heapq
import math

REF, ALT, AMBIG, NOOVERLAP = 0, 1, 2, 3


class MaxEditDistance(Exception):
    pass


class Variant:
    """HiPhase's Variant: position = 0-based VCF POS, alleles as in the VCF,
    allele0 < allele1 by index, a hom-alt site loaded as (REF, ALT)."""
    __slots__ = ("position", "ref_len", "allele0", "allele1", "prefix_len", "postfix_len", "vtype", "key", "index0")

    def __init__(self, position, ref_len, allele0, allele1, vtype, key=None, index0=0):
        self.index0 = index0
        self.position = position
        self.ref_len = ref_len
        self.allele0 = allele0
        self.allele1 = allele1
        self.prefix_len = 0
        self.postfix_len = 0
        self.vtype = vtype
        self.key = key

    def truncated0(self):
        return self.allele0[self.prefix_len:len(self.allele0) - self.postfix_len]

    def truncated1(self):
        return self.allele1[self.prefix_len:len(self.allele1) - self.postfix_len]


def variant_type(ref, alts):
    """get_variant_type (block_gen.rs) for records without SVTYPE / TRID."""
    ref_len = len(ref)
    max_alt = max(len(a) for a in alts)
    if ref_len == 1:
        return "snv" if max_alt == 1 else "insertion"
    return "deletion" if max_alt == 1 else "indel"


def load_variants(records, fetch_ref, reference_buffer=15):
    """phaser.rs: het variants (sorted alleles, both loaded) with the reference
    prefix / postfix buffer, and hom-alt variants (allele0 forced to REF, no
    buffer). records: iterable of (pos0, ref, alts, gt) sorted by position."""
    hets, homs = [], []
    previous_het_end = 0
    for pos0, ref, alts, gt in records:
        if len(gt) != 2 or None in gt:
            continue
        a, b = sorted(gt)
        if a == b == 0:
            continue
        hom = a == b
        if hom:
            a = 0
        alleles = [ref] + list(alts)
        vt = variant_type(ref, alts)
        v = Variant(pos0, len(ref), alleles[a].upper(), alleles[b].upper(), vt, key=(pos0, ref, tuple(alts), gt), index0=a)
        if reference_buffer > 0 and not hom:
            ref_prefix_start = max(0, pos0 - reference_buffer)
            ref_postfix_start = pos0 + len(ref)
            if ref_prefix_start < previous_het_end and hets:
                last = hets[-1]
                current_end = last.position + last.ref_len + last.postfix_len
                truncate = min(current_end - pos0, last.postfix_len)
                if truncate > 0:
                    last.allele0 = last.allele0[:len(last.allele0) - truncate]
                    last.allele1 = last.allele1[:len(last.allele1) - truncate]
                    last.postfix_len -= truncate
                ref_prefix_start = min(previous_het_end, pos0)
            prefix = fetch_ref(ref_prefix_start, pos0)
            v.allele0 = prefix + v.allele0
            v.allele1 = prefix + v.allele1
            v.prefix_len += len(prefix)
            postfix = fetch_ref(ref_postfix_start, ref_postfix_start + reference_buffer)
            v.allele0 += postfix
            v.allele1 += postfix
            v.postfix_len += len(postfix)
            previous_het_end = pos0 + len(ref)
        (homs if hom else hets).append(v)
    return hets, homs


class WFAGraph:
    def __init__(self, max_edit_distance):
        self.seqs = []
        self.edges = []
        self.max_edit_distance = max_edit_distance

    def add_node(self, seq, parents):
        idx = len(self.seqs)
        for p in parents:
            self.edges[p].append(idx)
        self.seqs.append(seq)
        self.edges.append([])
        return idx


def build_graph(reference, ref_offset, hets, homs, ref_start, ref_end, max_ed):
    """from_reference_variants_with_hom. reference[i - ref_offset] is reference
    base i (0-based); ref_start inclusive, ref_end exclusive."""
    def ref(a, b):
        return reference[a - ref_offset:b - ref_offset]
    g = WFAGraph(max_ed)
    node_to_alleles = {}
    previous_end = ref_start
    reference_reconnect = []
    reference_alleles = []
    queue = []  # (reconnect position, insertion counter, alt node)
    counter = 0
    allv = [(v, i) for i, v in enumerate(hets)] + [(v, None) for v in homs]
    allv.sort(key=lambda x: x[0].position)  # stable: hets before homs at one position

    def pop_reconnects_upto(limit):
        nonlocal previous_end, reference_reconnect, reference_alleles
        while queue and queue[0][0] <= limit:
            alt_reconnect, _, alt_index = heapq.heappop(queue)
            assert alt_reconnect > previous_end
            ri = g.add_node(ref(previous_end, alt_reconnect), reference_reconnect)
            if reference_alleles:
                node_to_alleles[ri] = reference_alleles
                reference_alleles = []
            previous_end = alt_reconnect
            reference_reconnect = [ri, alt_index]
            while queue and queue[0][0] == alt_reconnect:
                _, _, ai2 = heapq.heappop(queue)
                reference_reconnect.append(ai2)

    for v, vi in allv:
        vp, rl = v.position, v.ref_len
        if vp < ref_start or vp + rl > ref_end:
            continue
        pop_reconnects_upto(vp)
        if previous_end < vp or len(g.seqs) == 0:
            ri = g.add_node(ref(previous_end, vp), reference_reconnect)
            if reference_alleles:
                node_to_alleles[ri] = reference_alleles
                reference_alleles = []
            reference_reconnect = [ri]
            previous_end = vp
        else:
            assert previous_end == vp
        # convert_index(AlleleType::Reference) != 0: allele0 is an ALT (a 1/2 site);
        # hom-alt variants are loaded with allele0 forced to REF (index 0)
        if v.index0 != 0:
            ai = g.add_node(v.truncated0(), list(reference_reconnect))
            if vi is not None:
                node_to_alleles[ai] = [(vi, 0)]
            counter += 1
            heapq.heappush(queue, (vp + rl, counter, ai))
        else:
            if vi is not None:
                reference_alleles.append((vi, 0))
        ai = g.add_node(v.truncated1(), list(reference_reconnect))
        if vi is not None:
            node_to_alleles[ai] = [(vi, 1)]
        counter += 1
        heapq.heappush(queue, (vp + rl, counter, ai))
    # reconnect everything downstream
    while queue:
        alt_reconnect, _, alt_index = heapq.heappop(queue)
        assert alt_reconnect > previous_end
        ri = g.add_node(ref(previous_end, alt_reconnect), reference_reconnect)
        if reference_alleles:
            node_to_alleles[ri] = reference_alleles
            reference_alleles = []
        previous_end = alt_reconnect
        reference_reconnect = [ri, alt_index]
        while queue and queue[0][0] == alt_reconnect:
            _, _, ai2 = heapq.heappop(queue)
            reference_reconnect.append(ai2)
    assert previous_end <= ref_end
    g.add_node(ref(previous_end, ref_end), reference_reconnect)
    assert not reference_alleles
    return g, node_to_alleles


def edit_distance_with_pruning(g, other, prune_distance):
    """edit_distance_with_pruning: returns (score, sorted traversed nodes)."""
    n_nodes = len(g.seqs)
    other_len = len(other)
    set_index = {frozenset([0]): 0}
    index_set = [frozenset([0])]

    def intern(s):
        i = set_index.get(s)
        if i is None:
            i = len(index_set)
            index_set.append(s)
            set_index[s] = i
        return i

    active = {0: {0: [[0, 0]]}}
    nxt = {}
    maxfronts = {}
    ed = 0
    farthest = 0
    min_progression = 0
    while True:
        for node_index in range(n_nodes):
            if node_index not in active:
                continue
            seq = g.seqs[node_index]
            node_len = len(seq)
            wavefront = active.pop(node_index)
            maxfront = maxfronts.setdefault(node_index, {})
            for other_start, vec_waves in wavefront.items():
                max_offset = 0
                for w in vec_waves:
                    offset = w[0]
                    op = other_start + offset
                    # extend while bases match (chunked compare, same result as base by base)
                    while offset + 32 <= node_len and op + 32 <= other_len and seq[offset:offset + 32] == other[op:op + 32]:
                        offset += 32
                        op += 32
                    while offset < node_len and op < other_len and seq[offset] == other[op]:
                        offset += 1
                        op += 1
                    w[0] = offset
                    if offset > max_offset:
                        max_offset = offset
                rec = maxfront.get(other_start, 0)
                if max_offset < rec or (other_start + max_offset) < min_progression:
                    maxfront.setdefault(other_start, 0)
                    continue
                maxfront[other_start] = max_offset
                farthest = max(farthest, other_start + max_offset)
                best_sets = sorted({s for o, s in vec_waves if o == max_offset})
                if len(best_sets) > 1:
                    u = frozenset().union(*(index_set[s] for s in best_sets))
                    best_set = intern(u)
                else:
                    best_set = best_sets[0]
                if max_offset == node_len:
                    if node_index == n_nodes - 1:
                        if other_start + max_offset < other_len:
                            nxt.setdefault(node_index, {}).setdefault(other_start + 1, []).append([max_offset, best_set])
                    else:
                        new_offset = other_start + max_offset
                        cur = index_set[best_set]
                        for succ in g.edges[node_index]:
                            ns = intern(cur | {succ})
                            active.setdefault(succ, {}).setdefault(new_offset, []).append([0, ns])
                else:
                    node_wf = nxt.setdefault(node_index, {})
                    node_wf.setdefault(other_start - 1, []).append([max_offset + 1, best_set])
                    if other_start + max_offset < other_len:
                        node_wf.setdefault(other_start, []).append([max_offset + 1, best_set])
                        node_wf.setdefault(other_start + 1, []).append([max_offset, best_set])
            if node_index == n_nodes - 1:
                finals = set()
                for other_start, vec_waves in wavefront.items():
                    for offset, s in vec_waves:
                        if offset == node_len and other_start + offset == other_len:
                            finals.add(s)
                if finals:
                    if len(finals) > 1:
                        u = frozenset().union(*(index_set[s] for s in finals))
                    else:
                        u = index_set[next(iter(finals))]
                    return ed, sorted(u)
        ed += 1
        active = nxt
        nxt = {}
        if farthest > prune_distance:
            min_progression = farthest - prune_distance
        if ed > g.max_edit_distance:
            raise MaxEditDistance(g.max_edit_distance)


def aligned_regions(read):
    """global_realignment: reference ranges covered by the read's matches,
    split at reference skips (N)."""
    regions = []
    ref_pos = read.reference_start
    start = None
    end = 0
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            if start is None:
                start = ref_pos
            ref_pos += ln
            end = ref_pos
        elif op == 2:
            ref_pos += ln
        elif op == 3:
            if start is not None:
                regions.append((start, end))
                start = None
            ref_pos += ln
    if start is not None:
        regions.append((start, end))
    return regions


def global_realignment(read, hets, homs, reference, ref_offset,
                       max_edit_distance=500, max_ed_ratio=0.1, prune_distance=500):
    """Alleles per het variant (index into hets) for one read; raises
    MaxEditDistance as HiPhase does (the caller falls back to local)."""
    lookup = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
    seq = read.query_sequence.upper()
    alleles = [NOOVERLAP] * len(hets)
    total_overlaps = 0
    positions = [v.position for v in hets]
    hpos = [v.position for v in homs]
    import bisect
    for rs, re_ in aligned_regions(read):
        min_position, max_position = rs, re_ - 1
        first = bisect.bisect_left(positions, min_position)
        last = bisect.bisect_right(positions, max_position)
        if last <= first:
            continue
        hf = bisect.bisect_left(hpos, min_position)
        hl = bisect.bisect_right(hpos, max_position)
        read_start, read_end = lookup[min_position], lookup[max_position]
        read_align = seq[read_start:read_end + 1]
        max_region_ed = math.ceil(len(read_align) * max_ed_ratio)
        max_wfa_ed = min(max_edit_distance, max_region_ed)
        g, node_to_alleles = build_graph(reference, ref_offset, hets[first:last], homs[hf:hl],
                                         min_position, max_position + 1, max_wfa_ed)
        score, nodes = edit_distance_with_pruning(g, read_align, prune_distance)
        for ni in nodes:
            for vi, a in node_to_alleles.get(ni, []):
                ci = first + vi
                if alleles[ci] == NOOVERLAP:
                    alleles[ci] = a
                elif alleles[ci] != a:
                    alleles[ci] = AMBIG
        total_overlaps += last - first
    return alleles, total_overlaps


def edit_distance(v1, v2):
    l1 = len(v1)
    prev = list(range(l1 + 1))
    for i, c2 in enumerate(v2):
        row = [i + 1] + [0] * l1
        for j, c1 in enumerate(v1):
            row[j + 1] = min(prev[j + 1] + 1, row[j] + 1, prev[j] + (0 if c1 == c2 else 1))
        prev = row
    return prev[l1]


def local_realignment(read, hets):
    """local_realignment for small variants (SNV / insertion / deletion / indel)."""
    lookup = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
    seq = read.query_sequence.upper()
    min_position = read.reference_start
    max_position = max(lookup) if lookup else min_position
    out = []
    for v in hets:
        vp = v.position
        first_start = vp - v.prefix_len
        last_start = vp + 1
        first_end = vp + v.ref_len
        last_end = vp + v.ref_len + v.postfix_len + 1
        closest_start = next((lookup[sc] for sc in range(last_start - 1, first_start - 1, -1) if sc in lookup), None)
        closest_end = next((lookup[ec] for ec in range(first_end, last_end) if ec in lookup), None)
        start_coord = end_coord = None
        start_clip = end_clip = 0
        if closest_start is not None and closest_end is not None:
            for sc in range(first_start, last_start):
                start_clip += 1
                if sc in lookup:
                    si = lookup[sc]
                    if closest_start - si > 2 * v.prefix_len:
                        continue
                    start_coord = si
                    for ec in range(last_end - 1, first_end - 1, -1):
                        end_clip += 1
                        if ec in lookup:
                            ni = lookup[ec]
                            if ni - closest_end > 2 * v.postfix_len:
                                continue
                            end_coord = ni
                            break
                    break
        if start_coord is not None and end_coord is not None:
            obs = seq[start_coord:end_coord]
            if obs == v.allele0:
                a = REF
            elif obs == v.allele1:
                a = ALT
            else:
                hc, tc = start_clip - 1, end_clip - 1
                d0 = edit_distance(obs, v.allele0[hc:len(v.allele0) - tc])
                d1 = edit_distance(obs, v.allele1[hc:len(v.allele1) - tc])
                a = REF if d0 < d1 else ALT if d1 < d0 else AMBIG
            out.append(a)
        elif start_coord is not None:
            out.append(AMBIG)
        else:
            out.append(AMBIG if min_position <= vp <= max_position else NOOVERLAP)
    return out


def filter_out(read, min_mapq=5):
    return read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate or read.mapping_quality < min_mapq
