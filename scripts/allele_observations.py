"""Exact CIGAR-based sequence observations for multi-allelic VCF records."""


def observe_allele(alignment, start, reference, alleles):
    """Return an original VCF allele index, or None for partial/other sequence.

    Only the read sequence inside the reference replacement interval is used.
    Insertions after the VCF anchor are included. Equal-length alternate alleles
    are distinguished by sequence rather than by their shared length change.
    """
    end = start + len(reference)
    sequence = alignment.query_sequence
    if sequence is None or alignment.reference_start > start or \
            alignment.reference_end is None or alignment.reference_end < end:
        return None
    query_pos = 0
    ref_pos = alignment.reference_start
    observed = []
    for op, length in alignment.cigartuples or ():
        if op in (0, 7, 8):
            left, right = max(start, ref_pos), min(end, ref_pos + length)
            if left < right:
                observed.append(sequence[query_pos + left - ref_pos:
                                         query_pos + right - ref_pos])
            query_pos += length
            ref_pos += length
        elif op == 1:
            if start < ref_pos <= end:
                observed.append(sequence[query_pos:query_pos + length])
            query_pos += length
        elif op == 2:
            ref_pos += length
        elif op == 3:
            if ref_pos < end and ref_pos + length > start:
                return None
            ref_pos += length
        elif op == 4:
            query_pos += length
        if ref_pos > end:
            break
    observed = ''.join(observed).upper()
    matches = [i for i, allele in enumerate(alleles) if allele.upper() == observed]
    return matches[0] if len(matches) == 1 else None
