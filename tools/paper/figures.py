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
# Inputs and read alignment (other tools) on the left; inside the dashed pgphase boundary, the graph arm (blue
# lane) and union gap phasing (orange lane, only with a BAM) meet in the global EM, then windows and a second EM
# pass, then stitching, rescue and islands, and the outputs. Coordinates in cm, y upwards; a box is
# (x, y-centre, width, height, kind, lines).
LANES = {  # x0, y0, x1, y1, fill (HTML hex), label
    "graph": (6.15, 2.45, 12.35, 4.95, "E3EEF9", "graph arm"),
    "union": (6.15, 0.15, 12.35, 1.95, "FCEBDD", "union gap phasing (with a BAM)"),
}
PGPHASE_BOX = (5.85, -0.2, 21.0, 5.35)
BOXES = {
    "graph": (0.0, 4.15, 2.3, 0.85, "data", ["Pangenome graph", "(GBZ)"]),
    "reads": (0.0, 2.05, 2.3, 0.85, "data", ["Long reads"]),
    "gaf": (2.95, 3.0, 2.5, 0.85, "ext", ["Graph alignments", "(vg Giraffe; GAF)"]),
    "bam": (2.95, 1.05, 2.5, 0.85, "ext", ["Linear alignments", "(BAM)"]),
    "catalog": (6.45, 4.15, 2.5, 0.85, "step", ["Site catalog", "(build-snarl-catalog)"]),
    "index": (6.45, 3.0, 2.5, 0.85, "step", ["Indexed GAF", "(index-gaf)"]),
    "clean": (9.45, 3.55, 2.6, 1.05, "step", ["Catalog alleles from", "read walks; k-means"]),
    "solve": (9.45, 1.05, 2.6, 1.05, "step", ["Sample's own hets:", "candidates, MSA"]),
    "em": (12.75, 2.3, 2.6, 1.3, "step", ["Global EM:", "site errors,", "verified joins"]),
    "windows": (15.75, 3.2, 2.45, 0.95, "step", ["Local windows,", "EM pass 2"]),
    "fill": (15.75, 1.4, 2.45, 0.95, "step", ["Stitching, rescue,", "islands"]),
    "vcf": (18.65, 3.2, 2.1, 0.85, "out", ["Phased VCF"]),
    "tags": (18.65, 1.4, 2.1, 0.85, "out", ["Read HP/PS tags"]),
}
ARROWS = [("graph", "catalog"), ("graph", "gaf"), ("reads", "gaf"), ("reads", "bam"), ("gaf", "index"),
          ("catalog", "clean"), ("index", "clean"), ("bam", "solve"), ("clean", "em"), ("solve", "em"),
          ("em", "windows"), ("windows", "fill"), ("fill", "vcf"), ("fill", "tags")]
FILL = {"data": "EDEDED", "ext": "FFFFFF", "step": "FFFFFF", "out": "FBE3D6"}


def _edge(a, b):
    """An arrow between two boxes: right edge to left edge, or bottom to top in one column."""
    xa, ya, wa, ha, *_ = BOXES[a]
    xb, yb, wb, hb, *_ = BOXES[b]
    if abs(xa - xb) < 0.01:
        return (xa + wa / 2, ya - ha / 2), (xb + wb / 2, yb + hb / 2)
    return (xa + wa, ya), (xb, yb)


def _overview_tikz():
    out = []
    x0, y0, x1, y1 = PGPHASE_BOX
    out.append(f"\\draw[gray!70, dashed, rounded corners=4pt] ({x0},{y0}) rectangle ({x1},{y1});"
               f"\\node[anchor=north east, font=\\footnotesize\\bfseries, black!75] at ({x1 - 0.1},{y1 - 0.05}) {{pgphase}};")
    for key, (lx0, ly0, lx1, ly1, fill, label) in LANES.items():
        out.append(f"\\fill[fill={{rgb,255:red,{int(fill[0:2], 16)};green,{int(fill[2:4], 16)};blue,{int(fill[4:6], 16)}}}, "
                   f"rounded corners=3pt] ({lx0},{ly0}) rectangle ({lx1},{ly1});"
                   f"\\node[anchor={'north west' if key == 'graph' else 'south west'}, font=\\scriptsize\\itshape] at "
                   f"({lx0 + 0.1},{ly1 - 0.05 if key == 'graph' else ly0 + 0.03}) {{{label}}};")
    out.append(f"\\node[font=\\scriptsize\\itshape, black!70] at (4.2,{y1 - 0.1}) {{read alignment}};")
    for k, (x, y, w, h, kind, lines) in BOXES.items():
        fill = FILL[kind]
        rgb = f"{{rgb,255:red,{int(fill[0:2], 16)};green,{int(fill[2:4], 16)};blue,{int(fill[4:6], 16)}}}"
        style = {"data": "draw=black!60", "ext": "draw=black!60, dashed", "step": "draw=black!70",
                 "out": "draw=pghl, thick"}[kind]
        text = "\\\\".join(lines)
        out.append(f"\\node[{style}, fill={rgb}, rounded corners=2pt, align=center, font=\\footnotesize, inner sep=2pt, "
                   f"minimum width={w}cm, minimum height={h}cm, anchor=west] ({k}) at ({x},{y}) {{{text}}};")
    for a, b in ARROWS:
        (xa, ya), (xb, yb) = _edge(a, b)
        out.append(f"\\draw[-latex, thick, black!70] ({xa:.2f},{ya:.2f}) -- ({xb:.2f},{yb:.2f});")
    return "\n".join(out)


def _overview_svg():
    s = 44
    x0, y0, x1, y1 = PGPHASE_BOX
    W, H = (x1 + 0.4) * s, (y1 - y0 + 0.5) * s
    def X(x): return (x + 0.15) * s
    def Y(y): return (y1 + 0.2 - y) * s
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="pgphase pipeline">',
           '<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto">'
           '<path d="M0,0 L7,3 L0,6 z" fill="currentColor"/></marker></defs>']
    out.append(f'<rect x="{X(x0):.1f}" y="{Y(y1):.1f}" width="{(x1 - x0) * s:.1f}" height="{(y1 - y0) * s:.1f}" rx="7" '
               f'fill="none" stroke="currentColor" stroke-opacity=".5" stroke-dasharray="5 4"/>'
               f'<text x="{X(x1 - 0.15):.1f}" y="{Y(y1 - 0.4):.1f}" text-anchor="end" font-size="12" font-weight="700" fill="currentColor" fill-opacity=".85">pgphase</text>')
    for key, (lx0, ly0, lx1, ly1, fill, label) in LANES.items():
        ty = Y(ly1 - 0.35) if key == "graph" else Y(ly0 + 0.15)
        out.append(f'<rect x="{X(lx0):.1f}" y="{Y(ly1):.1f}" width="{(lx1 - lx0) * s:.1f}" height="{(ly1 - ly0) * s:.1f}" rx="5" fill="#{fill}"/>'
                   f'<text x="{X(lx0 + 0.12):.1f}" y="{ty:.1f}" font-size="10" font-style="italic" fill="#333">{label}</text>')
    out.append(f'<text x="{X(4.2):.1f}" y="{Y(y1 - 0.35):.1f}" text-anchor="middle" font-size="10" font-style="italic" fill="currentColor" fill-opacity=".7">read alignment</text>')
    for k, (x, y, w, h, kind, lines) in BOXES.items():
        stroke = {"data": "#777", "ext": "#777", "step": "#555", "out": "#D55E00"}[kind]
        dash = ' stroke-dasharray="4 3"' if kind == "ext" else ""
        width = 2 if kind == "out" else 1
        out.append(f'<rect x="{X(x):.1f}" y="{Y(y + h / 2):.1f}" width="{w * s:.1f}" height="{h * s:.1f}" rx="4" '
                   f'fill="#{FILL[kind]}" stroke="{stroke}" stroke-width="{width}"{dash}/>')
        for i, line in enumerate(lines):
            ty = Y(y) + (i - (len(lines) - 1) / 2) * 13 + 4
            out.append(f'<text x="{X(x + w / 2):.1f}" y="{ty:.1f}" text-anchor="middle" font-size="11" fill="#111">{line}</text>')
    for a, b in ARROWS:
        (xa, ya), (xb, yb) = _edge(a, b)
        out.append(f'<line x1="{X(xa):.1f}" y1="{Y(ya):.1f}" x2="{X(xb) - (2 if abs(xa - xb) > 0.01 else 0):.1f}" '
                   f'y2="{Y(yb) - (2 if abs(xa - xb) <= 0.01 else 0):.1f}" stroke="currentColor" stroke-width="1.3" marker-end="url(#ah)"/>')
    out.append("</svg>")
    return "\n".join(out)

# ------------------------------------------------------------------------------- accuracy
# One row per tool, read across: (a) two bar charts, misphased and unphased over the eight chromosomes, each from
# zero with its value printed, only pgphase coloured; (b) a heatmap of the share misphased or unphased per
# chromosome, each cell printing its value and shaded by it (log scale), development and untuned chromosomes grouped.
WIDE["accuracy"] = True
HIGHLIGHT, MUTED = "D55E00", "9A9A9A"
TIKZ_STYLES += f"\n\\definecolor{{pghl}}{{HTML}}{{{HIGHLIGHT}}}\\definecolor{{pgmuted}}{{HTML}}{{{MUTED}}}"

from tables import DEV, HELD, TOTALS  # noqa: E402

ORDER = sorted(TOTALS, key=lambda t: sum(TOTALS[t]))  # best first: top row
COLS = DEV + HELD
ROW, BAR_H = 0.5, 0.32           # row pitch and bar height (cm)
LAB = 2.2                        # tool labels, left of x = 0
BARS = [("Misphased (%)", 0, 5.0, 3.4), ("Unphased (%)", 1, 10.0, 3.4)]  # title, TOTALS index, axis max, width
GAP = 1.0                        # after each bar chart (room for its value labels)
CELL_W, HM_LO, HM_HI = 0.95, 0.3, 15.0
HM_X0 = sum(w for _, _, _, w in BARS) + GAP * len(BARS) + 0.4
TITLE_Y, HEAD_Y = 0.55, 0.45     # bar titles; heatmap column headers


def _row_y(i):
    return -i * ROW


def _shade(v):
    """White to Okabe-Ito blue on a log scale from HM_LO to HM_HI: (r, g, b) in 0-255, and a text colour."""
    f = max(0.0, min(1.0, (math.log10(v) - math.log10(HM_LO)) / (math.log10(HM_HI) - math.log10(HM_LO))))
    rgb = tuple(round(l + (h - l) * f) for l, h in zip((255, 255, 255), (0, 114, 178)))
    return rgb, ("white" if f > 0.55 else "black")


def _hm_value(t, p):
    return next(d + u for q, d, u in POINTS[t] if q == p)


def _cell_x(j):
    return HM_X0 + (j + (0.35 if j >= len(DEV) else 0)) * CELL_W


def _groups():
    return [("development", 0, len(DEV) - 1), ("not used for tuning", len(DEV), len(COLS) - 1)]


def _accuracy_tikz():
    out = []
    n = len(ORDER)
    bottom = _row_y(n - 1) - ROW / 2
    out.append(f"\\node[font=\\bfseries] at ({-LAB:.2f},{TITLE_Y + 0.45:.2f}) {{a}};")
    out.append(f"\\node[font=\\bfseries] at ({HM_X0 - 0.35:.2f},{TITLE_Y + 0.45:.2f}) {{b}};")
    for i, t in enumerate(ORDER):
        bold = "\\bfseries" if t == "pgphase" else ""
        out.append(f"\\node[left,font=\\small{bold}] at (-0.1,{_row_y(i):.3f}) {{{NAME[t]}}};")
    x0 = 0.0
    for title, k, top, w in BARS:
        tex_title = title.replace("%", "\\%")
        out.append(f"\\node[font=\\small\\bfseries] at ({x0 + w / 2:.2f},{TITLE_Y + 0.1:.2f}) {{{tex_title}}};")
        for tick in range(0, int(top) + 1, 1 if top <= 5 else 2):
            X = x0 + tick / top * w
            out.append(f"\\draw[gray!25] ({X:.3f},{TITLE_Y - 0.2:.2f}) -- ({X:.3f},{bottom:.3f});"
                       f"\\node[below,font=\\scriptsize] at ({X:.3f},{bottom:.3f}) {{{tick}}};")
        for i, t in enumerate(ORDER):
            v, y = TOTALS[t][k], _row_y(i)
            colour = "pghl" if t == "pgphase" else "pgmuted"
            bold = "\\bfseries" if t == "pgphase" else ""
            out.append(f"\\fill[{colour}] ({x0:.3f},{y - BAR_H / 2:.3f}) rectangle ({x0 + v / top * w:.3f},{y + BAR_H / 2:.3f});"
                       f"\\node[right,font=\\scriptsize{bold}] at ({x0 + v / top * w:.3f},{y:.3f}) {{{v:.2f}}};")
        x0 += w + GAP
    for label, ja, jb in _groups():
        xa, xb = _cell_x(ja), _cell_x(jb) + CELL_W
        out.append(f"\\draw[gray!60] ({xa + 0.05:.3f},{HEAD_Y + 0.32:.2f}) -- ({xb - 0.05:.3f},{HEAD_Y + 0.32:.2f});"
                   f"\\node[above,font=\\scriptsize\\itshape] at ({(xa + xb) / 2:.3f},{HEAD_Y + 0.3:.2f}) {{{label}}};")
    for j, p in enumerate(COLS):
        out.append(f"\\node[font=\\scriptsize] at ({_cell_x(j) + CELL_W / 2:.3f},{HEAD_Y:.2f}) {{{p}}};")
    for i, t in enumerate(ORDER):
        y = _row_y(i)
        bold = "\\bfseries" if t == "pgphase" else ""
        for j, p in enumerate(COLS):
            v = _hm_value(t, p)
            (r, g, b2), text = _shade(v)
            x = _cell_x(j)
            out.append(f"\\fill[fill={{rgb,255:red,{r};green,{g};blue,{b2}}}, draw=white, line width=1pt] "
                       f"({x:.3f},{y - ROW / 2:.3f}) rectangle ({x + CELL_W:.3f},{y + ROW / 2:.3f});"
                       f"\\node[font=\\scriptsize{bold},text={text}] at ({x + CELL_W / 2:.3f},{y:.3f}) {{{v:.1f}}};")
    i_pg = ORDER.index("pgphase")
    out.append(f"\\fill[pghl] ({HM_X0 - 0.16:.3f},{_row_y(i_pg) - ROW / 2 + 0.04:.3f}) rectangle "
               f"({HM_X0 - 0.06:.3f},{_row_y(i_pg) + ROW / 2 - 0.04:.3f});")
    out.append(f"\\node[font=\\scriptsize] at ({(HM_X0 + _cell_x(len(COLS) - 1) + CELL_W) / 2:.2f},{bottom - 0.3:.2f}) "
               f"{{misphased or unphased (\\%); darker is worse}};")
    return "\n".join(out)


def _accuracy_svg():
    s = 52
    n = len(ORDER)
    bottom = _row_y(n - 1) - ROW / 2
    W = (LAB + 0.3 + _cell_x(len(COLS) - 1) + CELL_W + 0.2) * s
    H = (TITLE_Y + 1.2 - bottom + 0.6) * s
    def X(x): return (x + LAB + 0.3) * s
    def Y(y): return (TITLE_Y + 1.1 - y) * s
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="Misphased and unphased reads per tool and chromosome">',
           '<style>.t{font-size:12px;fill:currentColor}.s{font-size:10.5px;fill:currentColor}.g{stroke:currentColor;stroke-opacity:.15}</style>']
    out.append(f'<text x="{X(-LAB - 0.2):.1f}" y="{Y(TITLE_Y + 0.4):.1f}" class="t" font-weight="700">a</text>')
    out.append(f'<text x="{X(HM_X0 - 0.4):.1f}" y="{Y(TITLE_Y + 0.4):.1f}" class="t" font-weight="700">b</text>')
    for i, t in enumerate(ORDER):
        bold = ' font-weight="700"' if t == "pgphase" else ""
        out.append(f'<text x="{X(-0.1):.1f}" y="{Y(_row_y(i)) + 4:.1f}" class="t" text-anchor="end"{bold}>{NAME[t]}</text>')
    x0 = 0.0
    for title, k, top, w in BARS:
        out.append(f'<text x="{X(x0 + w / 2):.1f}" y="{Y(TITLE_Y + 0.1) + 4:.1f}" class="t" text-anchor="middle" font-weight="700">{title}</text>')
        for tick in range(0, int(top) + 1, 1 if top <= 5 else 2):
            xx = X(x0 + tick / top * w)
            out.append(f'<line x1="{xx:.1f}" y1="{Y(TITLE_Y - 0.2):.1f}" x2="{xx:.1f}" y2="{Y(bottom):.1f}" class="g"/>'
                       f'<text x="{xx:.1f}" y="{Y(bottom) + 13:.1f}" class="s" text-anchor="middle">{tick}</text>')
        for i, t in enumerate(ORDER):
            v, y = TOTALS[t][k], _row_y(i)
            colour = HIGHLIGHT if t == "pgphase" else MUTED
            bold = ' font-weight="700"' if t == "pgphase" else ""
            out.append(f'<rect x="{X(x0):.1f}" y="{Y(y + BAR_H / 2):.1f}" width="{v / top * w * s:.1f}" height="{BAR_H * s:.1f}" fill="#{colour}"/>'
                       f'<text x="{X(x0 + v / top * w) + 4:.1f}" y="{Y(y) + 4:.1f}" class="s"{bold}>{v:.2f}</text>')
        x0 += w + GAP
    for label, ja, jb in _groups():
        xa, xb = _cell_x(ja), _cell_x(jb) + CELL_W
        out.append(f'<line x1="{X(xa + 0.05):.1f}" y1="{Y(HEAD_Y + 0.32):.1f}" x2="{X(xb - 0.05):.1f}" y2="{Y(HEAD_Y + 0.32):.1f}" stroke="currentColor" stroke-opacity=".5"/>'
                   f'<text x="{X((xa + xb) / 2):.1f}" y="{Y(HEAD_Y + 0.4):.1f}" class="s" text-anchor="middle" font-style="italic">{label}</text>')
    for j, p in enumerate(COLS):
        out.append(f'<text x="{X(_cell_x(j) + CELL_W / 2):.1f}" y="{Y(HEAD_Y) + 4:.1f}" class="s" text-anchor="middle">{p}</text>')
    for i, t in enumerate(ORDER):
        y = _row_y(i)
        bold = ' font-weight="700"' if t == "pgphase" else ""
        for j, p in enumerate(COLS):
            v = _hm_value(t, p)
            (r, g, b2), text = _shade(v)
            out.append(f'<g><title>{NAME[t]}, {p}: {v:.2f}% misphased or unphased</title>'
                       f'<rect x="{X(_cell_x(j)):.1f}" y="{Y(y + ROW / 2):.1f}" width="{CELL_W * s:.1f}" height="{ROW * s:.1f}" '
                       f'fill="rgb({r},{g},{b2})" stroke="white" stroke-width="2"/>'
                       f'<text x="{X(_cell_x(j) + CELL_W / 2):.1f}" y="{Y(y) + 4:.1f}" class="s" text-anchor="middle" fill="{text}"{bold}>{v:.1f}</text></g>')
    out.append(f'<text x="{X((HM_X0 + _cell_x(len(COLS) - 1) + CELL_W) / 2):.1f}" y="{Y(bottom - 0.3) + 4:.1f}" class="s" '
               f'text-anchor="middle">misphased or unphased (%); darker is worse</text>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "accuracy": _accuracy_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "accuracy": _accuracy_svg}[key]()
