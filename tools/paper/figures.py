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
# Two panels: (a) per tool, the phaseable arm reads misphased (dark) or unphased (light) over the eight
# chromosomes, as stacked bars; (b) per chromosome, each tool's misphased-or-unphased share as a dot (log scale).
# Okabe-Ito colours; pgphase in vermillion, drawn last and larger.
PALETTE = {"pgphase": "D55E00", "hiphase": "0072B2", "longphase": "E69F00", "whatshap": "009E73",
           "whatshap_opt": "56B4E9", "longcalld": "CC79A7"}
TIKZ_STYLES += "\n" + "\n".join(f"\\definecolor{{c{t.replace('_', '')}}}{{HTML}}{{{c}}}" for t, c in PALETTE.items())
WIDE["accuracy"] = True

from tables import DEV, HELD, TOTALS  # noqa: E402


def _ax(v, lo, hi, size):
    """Position of v on a log axis from lo to hi drawn over size."""
    return (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * size

A_W, ROW, BAR = 6.2, 0.62, 0.40         # panel a: width (cm) for 0-10%, row pitch, bar height
B_X0, B_W, B_LO, B_HI = 9.6, 6.8, 0.3, 20.0  # panel b: left edge, width, log range (%)
ORDER = sorted(TOTALS, key=lambda t: sum(TOTALS[t]))           # best first, top of panel a
ROWS = DEV + [None] + HELD                                     # panel b rows, a gap between the groups


def _b_x(v):
    return B_X0 + _ax(max(v, B_LO), B_LO, B_HI, B_W)


def _a_items():
    """Panel a: per tool its row's y and its two segments' widths."""
    n = len(ORDER)
    return [(t, (n - 1 - i) * ROW, TOTALS[t][0] / 10 * A_W, TOTALS[t][1] / 10 * A_W) for i, t in enumerate(ORDER)]


def _b_rows():
    """Panel b: per chromosome its row's y (top row first)."""
    n = len(ROWS)
    return [(p, (n - 1 - i) * ROW * 0.82) for i, p in enumerate(ROWS) if p]


def _accuracy_tikz():
    out = []
    n = len(ORDER)
    top = (n - 1) * ROW + 0.5
    # panel a
    out.append(f"\\node[font=\\bfseries] at (-2.3,{top + 0.35:.2f}) {{a}};")
    for x in range(0, 11, 2):
        X = x / 10 * A_W
        out.append(f"\\draw[gray!25] ({X:.2f},-0.35) -- ({X:.2f},{top - 0.2:.2f});"
                   f"\\node[below,font=\\scriptsize] at ({X:.2f},-0.35) {{{x}}};")
    out.append(f"\\node[font=\\small] at ({A_W / 2:.2f},-1.05) {{\\% of phaseable arm reads}};")
    for t, y, wd, wu in _a_items():
        c = "c" + t.replace("_", "")
        out.append(f"\\fill[{c}] (0,{y - BAR / 2:.3f}) rectangle ({wd:.3f},{y + BAR / 2:.3f});")
        out.append(f"\\fill[{c}!35] ({wd:.3f},{y - BAR / 2:.3f}) rectangle ({wd + wu:.3f},{y + BAR / 2:.3f});")
        weight = "\\bfseries" if t == "pgphase" else ""
        out.append(f"\\node[left,font=\\small{weight}] at (-0.05,{y:.3f}) {{{NAME[t]}}};")
        out.append(f"\\node[right,font=\\scriptsize{weight}] at ({wd + wu + 0.05:.3f},{y:.3f}) {{{sum(TOTALS[t]):.2f}\\%}};")
    ly = top + 0.05
    out.append(f"\\fill[black!75] (2.3,{ly - 0.1:.2f}) rectangle (2.6,{ly + 0.1:.2f});\\node[right,font=\\scriptsize] at (2.6,{ly:.2f}) {{misphased}};")
    out.append(f"\\fill[black!25] (4.2,{ly - 0.1:.2f}) rectangle (4.5,{ly + 0.1:.2f});\\node[right,font=\\scriptsize] at (4.5,{ly:.2f}) {{unphased}};")
    # panel b
    rows = _b_rows()
    btop = rows[0][1] + 0.5
    out.append(f"\\node[font=\\bfseries] at ({B_X0 - 1.3:.2f},{top + 0.35:.2f}) {{b}};")
    for v in (0.5, 1, 2, 5, 10, 20):
        X = _b_x(v)
        out.append(f"\\draw[gray!25] ({X:.3f},-0.35) -- ({X:.3f},{btop - 0.2:.2f});"
                   f"\\node[below,font=\\scriptsize] at ({X:.3f},-0.35) {{{v:g}}};")
    out.append(f"\\node[font=\\small] at ({B_X0 + B_W / 2:.2f},-1.05) {{misphased or unphased (\\%, log scale)}};")
    for p, y in rows:
        out.append(f"\\draw[gray!15, line width=3pt] ({B_X0:.2f},{y:.3f}) -- ({B_X0 + B_W:.2f},{y:.3f});"
                   f"\\node[left,font=\\small] at ({B_X0 - 0.05:.2f},{y:.3f}) {{{p}}};")
    gy = [y for p, y in rows]
    dev_mid = (rows[0][1] + rows[len(DEV) - 1][1]) / 2
    held_mid = (rows[len(DEV)][1] + rows[-1][1]) / 2
    for label, mid in (("development", dev_mid), ("not tuned", held_mid)):
        out.append(f"\\node[rotate=-90,font=\\scriptsize\\itshape,black!65] at ({B_X0 + B_W + 0.35:.2f},{mid:.3f}) {{{label}}};")
    for t in [t for t, _ in TOOLS if t != "pgphase"] + ["pgphase"]:
        c = "c" + t.replace("_", "")
        r = 0.11 if t == "pgphase" else 0.075
        for (p, d, u), (_, y) in zip(POINTS[t], [(q, yy) for q, yy in rows]):
            out.append(f"\\fill[{c}, draw=white, line width=0.3pt] ({_b_x(d + u):.3f},{y:.3f}) circle ({r});")
    return "\n".join(out)


def _accuracy_svg():
    s, pad = 52, 10
    n = len(ORDER)
    top = (n - 1) * ROW + 0.9
    W = (B_X0 + B_W + 0.9 + 2.6) * s
    H = (top + 1.4) * s
    def X(x): return (x + 2.6) * s
    def Y(y): return (top - y) * s
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="Misphased and unphased reads per tool">',
           '<style>.t{font-size:12px;fill:currentColor}.s{font-size:10.5px;fill:currentColor}.g{stroke:currentColor;stroke-opacity:.15}</style>']
    out.append(f'<text x="{X(-2.5):.1f}" y="{Y(top - 0.3):.1f}" class="t" font-weight="700">a</text>')
    for x in range(0, 11, 2):
        xx = X(x / 10 * A_W)
        out.append(f'<line x1="{xx:.1f}" y1="{Y(top - 0.55):.1f}" x2="{xx:.1f}" y2="{Y(-0.35):.1f}" class="g"/>'
                   f'<text x="{xx:.1f}" y="{Y(-0.35) + 14:.1f}" class="s" text-anchor="middle">{x}</text>')
    out.append(f'<text x="{X(A_W / 2):.1f}" y="{Y(-0.35) + 32:.1f}" class="t" text-anchor="middle">% of phaseable arm reads</text>')
    for t, y, wd, wu in _a_items():
        c = PALETTE[t]
        bold = ' font-weight="700"' if t == "pgphase" else ""
        out.append(f'<g><title>{NAME[t]}: {TOTALS[t][0]:.2f}% misphased, {TOTALS[t][1]:.2f}% unphased</title>'
                   f'<rect x="{X(0):.1f}" y="{Y(y + BAR / 2):.1f}" width="{wd * s:.1f}" height="{BAR * s:.1f}" fill="#{c}"/>'
                   f'<rect x="{X(wd):.1f}" y="{Y(y + BAR / 2):.1f}" width="{wu * s:.1f}" height="{BAR * s:.1f}" fill="#{c}" fill-opacity=".35"/></g>')
        out.append(f'<text x="{X(-0.08):.1f}" y="{Y(y) + 4:.1f}" class="t" text-anchor="end"{bold}>{NAME[t]}</text>')
        out.append(f'<text x="{X(wd + wu + 0.08):.1f}" y="{Y(y) + 4:.1f}" class="s"{bold}>{sum(TOTALS[t]):.2f}%</text>')
    ly = top - 0.45
    out.append(f'<rect x="{X(2.3):.1f}" y="{Y(ly) - 5:.1f}" width="14" height="10" fill="currentColor" fill-opacity=".75"/>'
               f'<text x="{X(2.3) + 18:.1f}" y="{Y(ly) + 4:.1f}" class="s">misphased</text>'
               f'<rect x="{X(4.2):.1f}" y="{Y(ly) - 5:.1f}" width="14" height="10" fill="currentColor" fill-opacity=".25"/>'
               f'<text x="{X(4.2) + 18:.1f}" y="{Y(ly) + 4:.1f}" class="s">unphased</text>')
    rows = _b_rows()
    out.append(f'<text x="{X(B_X0 - 1.3):.1f}" y="{Y(top - 0.3):.1f}" class="t" font-weight="700">b</text>')
    for v in (0.5, 1, 2, 5, 10, 20):
        xx = X(_b_x(v))
        out.append(f'<line x1="{xx:.1f}" y1="{Y(top - 0.55):.1f}" x2="{xx:.1f}" y2="{Y(-0.35):.1f}" class="g"/>'
                   f'<text x="{xx:.1f}" y="{Y(-0.35) + 14:.1f}" class="s" text-anchor="middle">{v:g}</text>')
    out.append(f'<text x="{X(B_X0 + B_W / 2):.1f}" y="{Y(-0.35) + 32:.1f}" class="t" text-anchor="middle">misphased or unphased (%, log scale)</text>')
    for p, y in rows:
        out.append(f'<line x1="{X(B_X0):.1f}" y1="{Y(y):.1f}" x2="{X(B_X0 + B_W):.1f}" y2="{Y(y):.1f}" stroke="currentColor" stroke-opacity=".08" stroke-width="7"/>'
                   f'<text x="{X(B_X0 - 0.08):.1f}" y="{Y(y) + 4:.1f}" class="t" text-anchor="end">{p}</text>')
    dev_mid = (rows[0][1] + rows[len(DEV) - 1][1]) / 2
    held_mid = (rows[len(DEV)][1] + rows[-1][1]) / 2
    for label, mid in (("development", dev_mid), ("not tuned", held_mid)):
        out.append(f'<text transform="translate({X(B_X0 + B_W + 0.35):.1f},{Y(mid):.1f}) rotate(90)" class="s" '
                   f'text-anchor="middle" font-style="italic" fill-opacity=".75">{label}</text>')
    for t in [t for t, _ in TOOLS if t != "pgphase"] + ["pgphase"]:
        r = 6 if t == "pgphase" else 4.2
        for (p, d, u), (_, y) in zip(POINTS[t], rows):
            out.append(f'<circle cx="{X(_b_x(d + u)):.1f}" cy="{Y(y):.1f}" r="{r}" fill="#{PALETTE[t]}" stroke="white" stroke-width=".6">'
                       f'<title>{NAME[t]}, {p}: {d + u:.2f}% misphased or unphased</title></circle>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "accuracy": _accuracy_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "accuracy": _accuracy_svg}[key]()
