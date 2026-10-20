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
LANES = {  # x0, y0, x1, y1, fill (HTML hex), label at the top-left corner
    "align": (2.75, 0.15, 5.65, 3.92, "F1F1F1", "read alignment (other tools)"),
    "graph": (6.15, 2.45, 12.35, 4.95, "E3EEF9", "graph arm"),
    "union": (6.15, 0.15, 12.35, 2.05, "FCEBDD", "union gap phasing (with a BAM)"),
}
PGPHASE_BOX = (5.85, -0.2, 21.0, 5.35)
BOXES = {
    "graph": (0.0, 4.15, 2.3, 0.85, "data", ["Pangenome graph", "(GBZ)"]),
    "reads": (0.0, 2.05, 2.3, 0.85, "data", ["Long reads"]),
    "gaf": (2.95, 3.0, 2.5, 0.85, "ext", ["Graph alignments", "(vg Giraffe; GAF)"]),
    "bam": (2.95, 0.95, 2.5, 0.85, "ext", ["Linear alignments", "(BAM)"]),
    "catalog": (6.45, 4.15, 2.5, 0.85, "step", ["Site catalog", "(build-snarl-catalog)"]),
    "index": (6.45, 3.0, 2.5, 0.85, "step", ["Indexed GAF", "(index-gaf)"]),
    "clean": (9.45, 3.55, 2.6, 1.05, "step", ["Catalog alleles from", "read walks; k-means"]),
    "solve": (7.9, 0.95, 2.9, 1.05, "step", ["Sample's own hets:", "candidates, MSA"]),
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
                   f"\\node[anchor=north west, font=\\scriptsize\\itshape, inner sep=2pt] at "
                   f"({lx0 + 0.05},{ly1 - 0.02}) {{{label}}};")
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
        ty = Y(ly1 - 0.33)
        out.append(f'<rect x="{X(lx0):.1f}" y="{Y(ly1):.1f}" width="{(lx1 - lx0) * s:.1f}" height="{(ly1 - ly0) * s:.1f}" rx="5" fill="#{fill}"/>'
                   f'<text x="{X(lx0 + 0.12):.1f}" y="{ty:.1f}" font-size="10" font-style="italic" fill="#333">{label}</text>')
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
# (a) Two bar charts, misphased and unphased over the eight chromosomes, one row per tool, values printed, only
# pgphase coloured. (b) One row per chromosome: the share of phaseable reads misphased or unphased on a linear
# axis, pgphase's dot joined to the best other tool's, the other tools faint, and how many times fewer reads
# pgphase leaves misphased or unphased printed at the right.
WIDE["accuracy"] = True
HIGHLIGHT, MUTED, BEST, FAINT = "D55E00", "9A9A9A", "4D4D4D", "C8C8C8"
TIKZ_STYLES += (f"\n\\definecolor{{pghl}}{{HTML}}{{{HIGHLIGHT}}}\\definecolor{{pgmuted}}{{HTML}}{{{MUTED}}}"
                f"\\definecolor{{pgbest}}{{HTML}}{{{BEST}}}\\definecolor{{pgfaint}}{{HTML}}{{{FAINT}}}")

from tables import DEV, HELD, TOTALS  # noqa: E402

ORDER = sorted(TOTALS, key=lambda t: sum(TOTALS[t]))  # panel a rows, best first
ROW, BAR_H = 0.5, 0.32            # panel a row pitch and bar height (cm)
LAB = 2.2                         # tool labels, left of x = 0
BARS = [("Misphased (%)", 0, 5.0, 3.2), ("Unphased (%)", 1, 10.0, 3.2)]  # title, TOTALS index, axis max, width
GAP = 0.95
TITLE_Y = 0.55
B_X0 = sum(w for _, _, _, w in BARS) + GAP * len(BARS) + 1.5   # panel b axis start (chromosome labels left of it)
B_W, B_MAX = 6.4, 15.0            # panel b axis width (cm) and range (%)
B_ROW = 0.42


def _row_y(i):
    return -i * ROW


def _err(t, p):
    return next(d + u for q, d, u in POINTS[t] if q == p)


def _b_rows():
    """Panel b: (chromosome, y), development first, a gap, then the chromosomes not used for tuning."""
    rows, y = [], 0.0
    for group in (DEV, HELD):
        for p in group:
            rows.append((p, y))
            y -= B_ROW
        y -= 0.45
    return rows


def _b_x(v):
    return B_X0 + v / B_MAX * B_W


def _fold(p):
    rivals = [t for t, _ in TOOLS if t != "pgphase"]
    best = min(rivals, key=lambda t: _err(t, p))
    return best, _err(best, p) / _err("pgphase", p)


def _accuracy_tikz():
    out = []
    n = len(ORDER)
    bottom = _row_y(n - 1) - ROW / 2
    out.append(f"\\node[font=\\bfseries] at ({-LAB:.2f},{TITLE_Y + 0.45:.2f}) {{a}};")
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
    # panel b
    rows = _b_rows()
    b_bottom = rows[-1][1] - 0.3
    out.append(f"\\node[font=\\bfseries] at ({B_X0 - 1.4:.2f},{TITLE_Y + 0.45:.2f}) {{b}};")
    out.append(f"\\node[font=\\small\\bfseries] at ({B_X0 + B_W / 2:.2f},{TITLE_Y + 0.1:.2f}) "
               f"{{Misphased or unphased (\\%), per chromosome}};")
    for tick in range(0, int(B_MAX) + 1, 5):
        X = _b_x(tick)
        out.append(f"\\draw[gray!25] ({X:.3f},{TITLE_Y - 0.2:.2f}) -- ({X:.3f},{b_bottom:.3f});"
                   f"\\node[below,font=\\scriptsize] at ({X:.3f},{b_bottom:.3f}) {{{tick}}};")
    out.append(f"\\node[font=\\scriptsize\\bfseries, anchor=west, pghl, align=left] at ({B_X0 + B_W + 0.1:.2f},{TITLE_Y + 0.05:.2f}) {{times\\\\fewer}};")
    groups = [("development", rows[0][1]), ("not used for tuning", rows[len(DEV)][1])]
    for label, y in groups:
        out.append(f"\\node[font=\\scriptsize\\itshape, anchor=south west, black!65] at ({B_X0 - 1.35:.2f},{y + 0.12:.3f}) {{{label}}};")
    for p, y in rows:
        best, fold = _fold(p)
        out.append(f"\\node[left,font=\\small] at ({B_X0 - 0.1:.2f},{y:.3f}) {{{p}}};")
        for t, _ in TOOLS:
            if t not in ("pgphase", best):
                out.append(f"\\fill[pgfaint] ({_b_x(_err(t, p)):.3f},{y:.3f}) circle (0.06);")
        out.append(f"\\draw[pgbest!60, line width=1.6pt] ({_b_x(_err('pgphase', p)):.3f},{y:.3f}) -- ({_b_x(_err(best, p)):.3f},{y:.3f});")
        out.append(f"\\fill[pgbest] ({_b_x(_err(best, p)):.3f},{y:.3f}) circle (0.085);")
        out.append(f"\\fill[pghl, draw=white, line width=0.4pt] ({_b_x(_err('pgphase', p)):.3f},{y:.3f}) circle (0.115);")
        out.append(f"\\node[font=\\small\\bfseries, anchor=west, pghl] at ({B_X0 + B_W + 0.15:.2f},{y:.3f}) {{{fold:.1f}\\texttimes}};")
    ly = b_bottom - 0.85
    lx = B_X0
    out.append(f"\\fill[pghl] ({lx:.2f},{ly:.2f}) circle (0.1);\\node[right,font=\\scriptsize] at ({lx + 0.12:.2f},{ly:.2f}) {{pgphase}};")
    out.append(f"\\fill[pgbest] ({lx + 1.7:.2f},{ly:.2f}) circle (0.08);\\node[right,font=\\scriptsize] at ({lx + 1.82:.2f},{ly:.2f}) {{best other tool}};")
    out.append(f"\\fill[pgfaint] ({lx + 4.1:.2f},{ly:.2f}) circle (0.06);\\node[right,font=\\scriptsize] at ({lx + 4.22:.2f},{ly:.2f}) {{other tools}};")
    return "\n".join(out)


def _accuracy_svg():
    s = 52
    n = len(ORDER)
    rows = _b_rows()
    bottom = min(_row_y(n - 1) - ROW / 2, rows[-1][1] - 0.3) - 1.2
    W = (LAB + 0.3 + B_X0 + B_W + 1.3) * s
    H = (TITLE_Y + 1.2 - bottom) * s
    def X(x): return (x + LAB + 0.3) * s
    def Y(y): return (TITLE_Y + 1.0 - y) * s
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="Misphased and unphased reads per tool and chromosome">',
           '<style>.t{font-size:12px;fill:currentColor}.s{font-size:10.5px;fill:currentColor}.g{stroke:currentColor;stroke-opacity:.15}</style>']
    out.append(f'<text x="{X(-LAB - 0.2):.1f}" y="{Y(TITLE_Y + 0.4):.1f}" class="t" font-weight="700">a</text>')
    for i, t in enumerate(ORDER):
        bold = ' font-weight="700"' if t == "pgphase" else ""
        out.append(f'<text x="{X(-0.1):.1f}" y="{Y(_row_y(i)) + 4:.1f}" class="t" text-anchor="end"{bold}>{NAME[t]}</text>')
    x0, a_bottom = 0.0, _row_y(n - 1) - ROW / 2
    for title, k, top, w in BARS:
        out.append(f'<text x="{X(x0 + w / 2):.1f}" y="{Y(TITLE_Y + 0.1) + 4:.1f}" class="t" text-anchor="middle" font-weight="700">{title}</text>')
        for tick in range(0, int(top) + 1, 1 if top <= 5 else 2):
            xx = X(x0 + tick / top * w)
            out.append(f'<line x1="{xx:.1f}" y1="{Y(TITLE_Y - 0.2):.1f}" x2="{xx:.1f}" y2="{Y(a_bottom):.1f}" class="g"/>'
                       f'<text x="{xx:.1f}" y="{Y(a_bottom) + 13:.1f}" class="s" text-anchor="middle">{tick}</text>')
        for i, t in enumerate(ORDER):
            v, y = TOTALS[t][k], _row_y(i)
            colour = HIGHLIGHT if t == "pgphase" else MUTED
            bold = ' font-weight="700"' if t == "pgphase" else ""
            out.append(f'<rect x="{X(x0):.1f}" y="{Y(y + BAR_H / 2):.1f}" width="{v / top * w * s:.1f}" height="{BAR_H * s:.1f}" fill="#{colour}"/>'
                       f'<text x="{X(x0 + v / top * w) + 4:.1f}" y="{Y(y) + 4:.1f}" class="s"{bold}>{v:.2f}</text>')
        x0 += w + GAP
    b_bottom = rows[-1][1] - 0.3
    out.append(f'<text x="{X(B_X0 - 1.5):.1f}" y="{Y(TITLE_Y + 0.4):.1f}" class="t" font-weight="700">b</text>')
    out.append(f'<text x="{X(B_X0 + B_W / 2):.1f}" y="{Y(TITLE_Y + 0.1) + 4:.1f}" class="t" text-anchor="middle" font-weight="700">Misphased or unphased (%), per chromosome</text>')
    for tick in range(0, int(B_MAX) + 1, 5):
        xx = X(_b_x(tick))
        out.append(f'<line x1="{xx:.1f}" y1="{Y(TITLE_Y - 0.2):.1f}" x2="{xx:.1f}" y2="{Y(b_bottom):.1f}" class="g"/>'
                   f'<text x="{xx:.1f}" y="{Y(b_bottom) + 13:.1f}" class="s" text-anchor="middle">{tick}</text>')
    out.append(f'<text x="{X(B_X0 + B_W + 0.15):.1f}" y="{Y(TITLE_Y - 0.1) + 4:.1f}" class="s" font-weight="700" fill="#{HIGHLIGHT}">times fewer</text>')
    for label, y in (("development", rows[0][1]), ("not used for tuning", rows[len(DEV)][1])):
        out.append(f'<text x="{X(B_X0 - 1.35):.1f}" y="{Y(y + 0.18):.1f}" class="s" font-style="italic" fill-opacity=".7">{label}</text>')
    for p, y in rows:
        best, fold = _fold(p)
        out.append(f'<text x="{X(B_X0 - 0.1):.1f}" y="{Y(y) + 4:.1f}" class="t" text-anchor="end">{p}</text>')
        for t, lab in TOOLS:
            if t not in ("pgphase", best):
                out.append(f'<circle cx="{X(_b_x(_err(t, p))):.1f}" cy="{Y(y):.1f}" r="3.2" fill="#{FAINT}"><title>{lab}, {p}: {_err(t, p):.2f}%</title></circle>')
        out.append(f'<line x1="{X(_b_x(_err("pgphase", p))):.1f}" y1="{Y(y):.1f}" x2="{X(_b_x(_err(best, p))):.1f}" y2="{Y(y):.1f}" stroke="#{BEST}" stroke-opacity=".6" stroke-width="3"/>')
        out.append(f'<circle cx="{X(_b_x(_err(best, p))):.1f}" cy="{Y(y):.1f}" r="4.5" fill="#{BEST}"><title>{NAME[best]}, {p}: {_err(best, p):.2f}%</title></circle>')
        out.append(f'<circle cx="{X(_b_x(_err("pgphase", p))):.1f}" cy="{Y(y):.1f}" r="6" fill="#{HIGHLIGHT}" stroke="white" stroke-width=".8"><title>pgphase, {p}: {_err("pgphase", p):.2f}%</title></circle>')
        out.append(f'<text x="{X(B_X0 + B_W + 0.15):.1f}" y="{Y(y) + 4:.1f}" class="t" font-weight="700" fill="#{HIGHLIGHT}">{fold:.1f}×</text>')
    ly, lx = Y(b_bottom - 0.85), B_X0
    for dx, col, r, lab in ((0, HIGHLIGHT, 5, "pgphase"), (1.7, BEST, 4.2, "best other tool"), (4.1, FAINT, 3.2, "other tools")):
        out.append(f'<circle cx="{X(lx + dx):.1f}" cy="{ly:.1f}" r="{r}" fill="#{col}"/><text x="{X(lx + dx + 0.15):.1f}" y="{ly + 4:.1f}" class="s">{lab}</text>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "accuracy": _accuracy_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "accuracy": _accuracy_svg}[key]()
