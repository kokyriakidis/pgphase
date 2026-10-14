import itertools, random, sys
sys.path.insert(0, ".")
import hiphase_port as H


def ed(a, b):
    return H.edit_distance(a, b)


def paths(ref, start, end, variants):
    """every combination of allele choices whose spans do not overlap (HiPhase's
    graph cannot traverse overlapping variants), as (sequence, choices)."""
    out = []
    for choice in itertools.product((0, 1), repeat=len(variants)):
        used = []
        for v, c in zip(variants, choice):
            seqa = v.truncated0() if c == 0 else v.truncated1()
            is_ref = c == 0 and v.index0 == 0
            if not is_ref:
                used.append((v.position, v.position + v.ref_len, seqa))
        used.sort()
        if any(used[i][1] > used[i + 1][0] for i in range(len(used) - 1)):
            continue
        s, cur = [], start
        for a, b, seqa in used:
            s.append(ref[cur:a]); s.append(seqa); cur = b
        s.append(ref[cur:end])
        out.append(("".join(s), choice))
    return out


rng = random.Random(7)
bad = 0
for trial in range(300):
    L = 80
    ref = "".join(rng.choice("ACGT") for _ in range(L))
    vs = []
    used = set()
    for k in range(rng.randint(1, 4)):
        p = rng.randint(5, L - 15)
        if any(abs(p - q) < 4 for q in used):
            continue
        used.add(p)
        kind = rng.choice(["snv", "ins", "del", "multi"])
        r = ref[p]
        if kind == "snv":
            alt = rng.choice([b for b in "ACGT" if b != r]); refa = r; gt = (0, 1); alts = [alt]
        elif kind == "ins":
            refa = r; alts = [r + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4)))]; gt = (0, 1)
        elif kind == "del":
            n = rng.randint(2, 4); refa = ref[p:p + n]; alts = [r]; gt = (0, 1)
        else:
            refa = r; alts = [r + "A", r + "AA"]; gt = (1, 2)
        vs.append((p, refa, alts, gt))
    vs.sort()
    hets, homs = H.load_variants(vs, lambda a, b: ref[a:b], reference_buffer=15)
    # a read: one path with a couple of random errors
    allp = paths(ref, 0, L, hets)
    seq, choice = rng.choice(allp)
    seq = list(seq)
    for _ in range(rng.randint(0, 2)):
        i = rng.randrange(len(seq))
        seq[i] = rng.choice("ACGT")
    seq = "".join(seq)
    g, n2a = H.build_graph(ref, 0, hets, homs, 0, L, 1000)
    score, nodes = H.edit_distance_with_pruning(g, seq, 10 ** 9)
    brute = min(ed(seq, p) for p, c in allp)
    if score != brute:
        bad += 1
        print("MISMATCH", trial, score, brute, vs)
        continue
    # alleles: every variant whose optimal paths all agree must match the port's call
    best = [c for p, c in allp if ed(seq, p) == brute]
    alle = [H.NOOVERLAP] * len(hets)
    for ni in nodes:
        for vi, a in n2a.get(ni, []):
            alle[vi] = a if alle[vi] in (H.NOOVERLAP, a) else H.AMBIG
    for vi in range(len(hets)):
        opts = {c[vi] for c in best}
        if len(opts) == 1 and alle[vi] not in (next(iter(opts)),):
            bad += 1
            print("ALLELE", trial, vi, opts, alle[vi], vs)
print("trials 300, problems", bad)
