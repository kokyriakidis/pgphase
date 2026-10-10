def costs_and_counts(a, b):
    costs = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    ways = [[1] * (len(b) + 1) for _ in range(len(a) + 1)]
    costs[0] = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        costs[i][0] = i
        for j in range(1, len(b) + 1):
            edges = [(costs[i-1][j]+1, ways[i-1][j]),
                     (costs[i][j-1]+1, ways[i][j-1]),
                     (costs[i-1][j-1]+(a[i-1]!=b[j-1]), ways[i-1][j-1])]
            best = min(c for c, _ in edges)
            costs[i][j] = best
            ways[i][j] = sum(n for c, n in edges if c == best)
    return costs, ways

def oracle(ref, alt):
    costs, ways = costs_and_counts(ref, alt)
    backwards, suffixes = costs_and_counts(ref[::-1], alt[::-1])
    optimum, total = costs[-1][-1], ways[-1][-1]
    matches = []
    for i, base in enumerate(ref):
        for j, alt_base in enumerate(alt):
            if base == alt_base and costs[i][j] + backwards[len(ref)-i-1][len(alt)-j-1] == optimum:
                through = ways[i][j] * suffixes[len(ref)-i-1][len(alt)-j-1]
                if through == total:
                    matches.append((i, j))
    spans = []
    for i, j in matches:
        if spans and spans[-1][0] + spans[-1][2] == i and spans[-1][1] + spans[-1][2] == j:
            spans[-1][2] += 1
        else:
            spans.append([i, j, 1])
    return optimum, spans, total

