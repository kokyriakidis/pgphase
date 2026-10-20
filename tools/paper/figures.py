"""The paper's figures, each defined once and drawn twice: TikZ for paper/ and SVG for the HTML draft.

  overview   the pipeline, boxes and arrows
  accuracy   per tool and chromosome, discordant against unphased phaseable arm reads (log-log)
"""
import math

from tables import CHROMS, NAME, POINTS, TOOLS

WIDE = {"overview": True, "accuracy": False}

# One colour and marker per tool, shared by both renderings.
STYLE = {
    "pgphase": ("black", "000000", "circle"),
    "hiphase": ("blue!75!black", "1f3fbf", "square"),
    "longphase": ("orange!85!black", "d97a00", "diamond"),
    "whatshap": ("red!75!black", "bf1f1f", "triangle"),
    "whatshap_opt": ("violet", "8f00ff", "itriangle"),
    "longcalld": ("green!50!black", "1f7f1f", "cross"),
}
TIKZ_STYLES = "\\usetikzlibrary{plotmarks}\n\\tikzset{box/.style={draw, rounded corners=2pt, align=center, font=\\footnotesize, inner sep=3pt}," \
              " inbox/.style={box, fill=gray!10}, outbox/.style={box, fill=gray!10, very thick}, arr/.style={-latex, thick}}"

# ------------------------------------------------------------------------------- overview
# Boxes: key -> (x, y, width, height, kind, lines); coordinates in cm, y upwards.
BOXES = {
    "graph": (0.0, 4.2, 2.6, 0.9, "in", ["Pangenome graph", "(GBZ)"]),
    "reads": (0.0, 2.3, 2.6, 0.9, "in", ["Long reads"]),
    "catalog": (3.4, 4.2, 2.8, 0.9, "step", ["Site catalog", "(snarl alleles, pgcat)"]),
    "gaf": (3.4, 2.9, 2.8, 0.9, "step", ["Graph alignments", "(GAF, indexed)"]),
    "bam": (3.4, 1.6, 2.8, 0.9, "step", ["Linear alignments", "(BAM)"]),
    "clean": (7.0, 3.55, 3.2, 1.2, "step", ["Catalog sites phased", "from walk observations", "(k-means)"]),
    "solve": (7.0, 1.6, 3.2, 1.2, "step", ["Sample-specific hets", "(candidates, noisy-region", "MSA)"]),
    "em": (11.0, 2.6, 3.3, 1.6, "step", ["Global EM per chunk:", "site phases and errors,", "joins the reads verify"]),
    "windows": (15.1, 3.45, 3.0, 1.0, "step", ["Haplotype and allele", "windows, EM pass 2"]),
    "fill": (15.1, 1.75, 3.0, 1.0, "step", ["Stitching, rescue,", "island blocks"]),
    "out": (18.9, 2.6, 2.4, 1.3, "out", ["Phased VCF", "Haplotagged", "reads (HP/PS)"]),
}
ARROWS = [("graph", "catalog"), ("reads", "gaf"), ("reads", "bam"), ("graph", "gaf"), ("catalog", "clean"),
          ("gaf", "clean"), ("bam", "solve"), ("clean", "em"), ("solve", "em"), ("em", "windows"), ("windows", "fill"),
          ("em", "fill"), ("fill", "out")]


def _edge(a, b):
    """Start and end points of an arrow between two boxes: right edge to left edge, or bottom to top in a column."""
    xa, ya, wa, ha, *_ = BOXES[a]
    xb, yb, wb, hb, *_ = BOXES[b]
    if abs(xa - xb) < 0.01:
        top, bottom = (a, b) if ya > yb else (b, a)
        xt, yt, wt, ht, *_ = BOXES[top]
        xb2, yb2, wb2, hb2, *_ = BOXES[bottom]
        p, q = (xt + wt / 2, yt - ht / 2), (xb2 + wb2 / 2, yb2 + hb2 / 2)
        return (p, q) if top == a else (q, p)
    return (xa + wa, ya), (xb, yb)


def _overview_tikz():
    out = []
    for k, (x, y, w, h, kind, lines) in BOXES.items():
        style = {"in": "inbox", "out": "outbox", "step": "box"}[kind]
        text = "\\\\".join(lines)
        out.append(f"\\node[{style}, minimum width={w}cm, minimum height={h}cm, anchor=west] ({k}) at ({x},{y}) "
                   f"{{{text}}};")
    for a, b in ARROWS:
        (x1, y1), (x2, y2) = _edge(a, b)
        out.append(f"\\draw[arr] ({x1:.2f},{y1:.2f}) -- ({x2:.2f},{y2:.2f});")
    out.append("\\node[font=\\scriptsize\\itshape, anchor=west] at (7.0,5.0) {graph arm: known alleles};")
    out.append("\\node[font=\\scriptsize\\itshape, anchor=west] at (7.0,0.7) {union gap phasing: the sample's own hets};")
    return "\n".join(out)


def _overview_svg():
    s, H = 40, 5.6  # px per cm, height in cm
    W = max(x + w for x, _, w, *_ in BOXES.values()) + 0.3
    def X(x): return x * s + 8
    def Y(y): return (H - y) * s
    out = [f'<svg viewBox="0 0 {W * s + 16:.0f} {H * s + 10:.0f}" class="chart" role="img" aria-label="pgphase overview">',
           '<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto">'
           '<path d="M0,0 L7,3 L0,6 z" fill="currentColor"/></marker></defs>']
    for k, (x, y, w, h, kind, lines) in BOXES.items():
        fill = "var(--panel, #f2f2f2)" if kind in ("in", "out") else "none"
        sw = 2 if kind == "out" else 1
        out.append(f'<rect x="{X(x):.1f}" y="{Y(y + h / 2):.1f}" width="{w * s:.1f}" height="{h * s:.1f}" rx="4" '
                   f'fill="{fill}" stroke="currentColor" stroke-width="{sw}"/>')
        for i, line in enumerate(lines):
            ty = Y(y) + (i - (len(lines) - 1) / 2) * 13 + 4
            out.append(f'<text x="{X(x + w / 2):.1f}" y="{ty:.1f}" text-anchor="middle" font-size="11" fill="currentColor">{line}</text>')
    for a, b in ARROWS:
        (x1, y1), (x2, y2) = _edge(a, b)
        out.append(f'<line x1="{X(x1):.1f}" y1="{Y(y1):.1f}" x2="{X(x2) - 2:.1f}" y2="{Y(y2):.1f}" stroke="currentColor" '
                   f'stroke-width="1.3" marker-end="url(#ah)"/>')
    out.append(f'<text x="{X(7.0):.1f}" y="{Y(5.0):.1f}" font-size="10" font-style="italic" fill="currentColor">graph arm: known alleles</text>')
    out.append(f'<text x="{X(7.0):.1f}" y="{Y(0.7):.1f}" font-size="10" font-style="italic" fill="currentColor">union gap phasing: the sample\'s own hets</text>')
    out.append("</svg>")
    return "\n".join(out)

# ------------------------------------------------------------------------------- accuracy
LO_X, HI_X, LO_Y, HI_Y = 0.2, 20.0, 0.02, 20.0
TICKS = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]


def _ax(v, lo, hi, size):
    return (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * size


def _accuracy_tikz():
    w, h = 7.0, 5.6
    out = [f"\\draw[->] (0,0) -- ({w + 0.2},0);", f"\\draw[->] (0,0) -- (0,{h + 0.2});"]
    for t in TICKS:
        if LO_X <= t <= HI_X:
            x = _ax(t, LO_X, HI_X, w)
            out.append(f"\\draw[gray!25] ({x:.3f},0) -- ({x:.3f},{h});\\node[below,font=\\scriptsize] at ({x:.3f},0) {{{t:g}}};")
        if LO_Y <= t <= HI_Y:
            y = _ax(t, LO_Y, HI_Y, h)
            out.append(f"\\draw[gray!25] (0,{y:.3f}) -- ({w},{y:.3f});\\node[left,font=\\scriptsize] at (0,{y:.3f}) {{{t:g}}};")
    out.append(f"\\node[font=\\small] at ({w / 2},-0.75) {{discordant phaseable reads (\\%)}};")
    out.append(f"\\node[font=\\small,rotate=90] at (-1.0,{h / 2}) {{unphased phaseable reads (\\%)}};")
    marks = {"circle": "*", "square": "square*", "diamond": "diamond*", "triangle": "triangle*", "itriangle": "triangle*",
             "cross": "x"}
    for i, (t, lab) in enumerate(TOOLS):
        colour, _, shape = STYLE[t]
        rot = ",rotate=180" if shape == "itriangle" else ""
        for _, dx, uy in POINTS[t]:
            x, y = _ax(max(dx, LO_X), LO_X, HI_X, w), _ax(max(uy, LO_Y), LO_Y, HI_Y, h)
            out.append(f"\\node[{colour},inner sep=0pt,font=\\scriptsize{rot}] at ({x:.3f},{y:.3f}) {{\\pgfuseplotmark{{{marks[shape]}}}}};")
        ly = h - 0.35 - i * 0.4
        out.append(f"\\node[{colour},inner sep=0pt{rot}] at ({w + 0.6:.2f},{ly:.2f}) {{\\pgfuseplotmark{{{marks[shape]}}}}};"
                   f"\\node[right,font=\\scriptsize] at ({w + 0.75:.2f},{ly:.2f}) {{{lab}}};")
    return "\n".join(out)


def _accuracy_svg():
    w, h, s, ox, oy = 7.0, 5.6, 60, 70, 20
    def X(v): return ox + _ax(max(v, LO_X), LO_X, HI_X, w) * s
    def Y(v): return oy + (h - _ax(max(v, LO_Y), LO_Y, HI_Y, h)) * s
    W, H = ox + w * s + 170, oy + h * s + 55
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="Discordant against unphased reads">']
    for t in TICKS:
        if LO_X <= t <= HI_X:
            out.append(f'<line x1="{X(t):.1f}" y1="{oy}" x2="{X(t):.1f}" y2="{oy + h * s}" stroke="currentColor" stroke-opacity="0.15"/>'
                       f'<text x="{X(t):.1f}" y="{oy + h * s + 15}" font-size="11" text-anchor="middle" fill="currentColor">{t:g}</text>')
        if LO_Y <= t <= HI_Y:
            out.append(f'<line x1="{ox}" y1="{Y(t):.1f}" x2="{ox + w * s}" y2="{Y(t):.1f}" stroke="currentColor" stroke-opacity="0.15"/>'
                       f'<text x="{ox - 6}" y="{Y(t) + 4:.1f}" font-size="11" text-anchor="end" fill="currentColor">{t:g}</text>')
    out.append(f'<text x="{ox + w * s / 2:.0f}" y="{oy + h * s + 38}" font-size="12" text-anchor="middle" fill="currentColor">discordant phaseable reads (%)</text>')
    out.append(f'<text transform="translate(16,{oy + h * s / 2:.0f}) rotate(-90)" font-size="12" text-anchor="middle" fill="currentColor">unphased phaseable reads (%)</text>')

    def mark(shape, cx, cy, col):
        r = 4.5
        if shape == "circle":
            return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="#{col}"><title></title></circle>'
        if shape == "square":
            return f'<rect x="{cx - r:.1f}" y="{cy - r:.1f}" width="{2 * r}" height="{2 * r}" fill="#{col}"/>'
        if shape == "diamond":
            return f'<path d="M{cx},{cy - r - 1} L{cx + r + 1},{cy} L{cx},{cy + r + 1} L{cx - r - 1},{cy} z" fill="#{col}"/>'
        if shape in ("triangle", "itriangle"):
            d = -1 if shape == "triangle" else 1
            return f'<path d="M{cx},{cy + d * (r + 1)} L{cx + r + 1},{cy - d * r} L{cx - r - 1},{cy - d * r} z" fill="#{col}"/>'
        return (f'<path d="M{cx - r},{cy - r} L{cx + r},{cy + r} M{cx - r},{cy + r} L{cx + r},{cy - r}" stroke="#{col}" stroke-width="2"/>')

    for i, (t, lab) in enumerate(TOOLS):
        _, col, shape = STYLE[t]
        for p, dx, uy in POINTS[t]:
            out.append(f'<g><title>{lab}, {p}: {dx:.2f}% discordant, {uy:.2f}% unphased</title>{mark(shape, X(dx), Y(uy), col)}</g>')
        ly = oy + 12 + i * 20
        out.append(mark(shape, ox + w * s + 25, ly, col) +
                   f'<text x="{ox + w * s + 37}" y="{ly + 4}" font-size="12" fill="currentColor">{lab}</text>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "accuracy": _accuracy_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "accuracy": _accuracy_svg}[key]()
