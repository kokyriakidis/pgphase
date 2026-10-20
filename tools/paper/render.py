"""Renders the paper and its supplement: the LaTeX of paper/, an Overleaf project, and the HTML pages a draft is
read in.

    python3 tools/paper/render.py    # rewrites paper/main.tex, paper/sections/, paper/figures/, paper/tables/,
                                     # paper/references.bib and paper/supplement.tex; writes build/paper/paper.html
                                     # and supplement.html, and compiles build/paper/pgphase.pdf and supplement.pdf

paper/ holds only what Overleaf compiles: main.tex (the journal's front matter, then one \\input per section),
supplement.tex, the sections, figures and tables, references.bib, and OUP's General Template class
(oup-authoring-template.cls, January 2026) with its author-year bibliography style (oup-abbrvnat.bst), as in the
dinara-align paper. Everything in it but those two is generated: the text from main_text.py and supplement.py,
the tables and figure data from results.json through tables.py, the references from refs.py.
"""
import html as htmlmod
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PAPER = os.path.join(ROOT, "paper")
HTML_OUT = os.path.join(ROOT, "build", "paper")
sys.path.insert(0, HERE)
import figures as F
import main_text as M
import refs as R
import supplement as S

CITE = {k: R.cite(k) for k in R.REFS}

# ------------------------------------------------------------------------------------------- numbering
def number_items(body, tables):
    """Tables and figures numbered in order of first appearance."""
    tabs, figs = {}, {}
    for b in body:
        if b[0] == "table" and b[1] not in tabs:
            tabs[b[1]] = tables[b[1]].get("label", str(len(tabs) + 1))
        if b[0] == "figure":
            figs[b[1]] = str(len(figs) + 1)
    return tabs, figs

# ------------------------------------------------------------------------------------------- inline
UNIT = r"(µs|ms|ns|s|kbp|Mbp|kb|Mb|bp|MB|GB|KB|GHz)"

def split_math(text):
    """Pieces of text and maths, the maths kept as written."""
    parts = re.split(r"(\\\(.*?\\\)|\\\[.*?\\\])", text)
    return [(p, p.startswith("\\(") or p.startswith("\\[")) for p in parts if p]

def cite_html(keys, narrative):
    out = []
    for k in keys:
        name, year = CITE[k]
        out.append(f'<a class="cite" href="#ref-{k}">{name} ({year})</a>' if narrative else f'<a class="cite" href="#ref-{k}">{name}, {year}</a>')
    return ", ".join(out) if narrative else "(" + "; ".join(out) + ")"

def inline_html(text, tabs, figs):
    out = []
    for piece, is_math in split_math(text):
        if is_math:
            out.append(piece); continue
        t = htmlmod.escape(piece, quote=False)
        t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)
        t = re.sub(r"\[@!([\w]+)\]", lambda m: cite_html([m.group(1)], True), t)
        t = re.sub(r"\[(@[\w]+(?:;\s*@[\w]+)*)\]", lambda m: cite_html([k.strip()[1:] for k in m.group(1).split(";")], False), t)
        t = re.sub(r"\[#tab:(\w+)\]", lambda m: f'<a href="#tab-{m.group(1)}">Table {tabs[m.group(1)]}</a>', t)
        t = re.sub(r"\[#fig:(\w+)\]", lambda m: f'<a href="#fig-{m.group(1)}">Fig. {figs[m.group(1)]}</a>', t)
        t = re.sub(r"\^([\w/]+)", r"<sup>\1</sup>", t)
        t = re.sub(r"(\d) " + UNIT + r"\b", r"\1&nbsp;\2", t)
        t = re.sub(r"(https://[^\s<]+[^\s<.,;)])", r'<a href="\1">\1</a>', t)
        out.append(t)
    return "".join(out)

def tex_escape(t):
    t = t.replace("\\", "\\textbackslash{}")
    for a, b in (("&", "\\&"), ("%", "\\%"), ("#", "\\#"), ("_", "\\_"), ("$", "\\$")):
        t = t.replace(a, b)
    return t

def inline_tex(text):
    out = []
    for piece, is_math in split_math(text):
        if is_math:
            out.append(piece); continue
        codes = []
        def keep_code(m):
            codes.append(m.group(1)); return f"\x00{len(codes) - 1}\x00"
        t = re.sub(r"`([^`]+)`", keep_code, piece)
        refs = []
        def keep_ref(m):
            refs.append(m.group(0)); return f"\x01{len(refs) - 1}\x01"
        t = re.sub(r"\[@!?[\w]+(?:;\s*@[\w]+)*\]|\[#(?:tab|fig):\w+\]", keep_ref, t)
        urls = []
        def keep_url(m):
            urls.append(m.group(1)); return f"\x02{len(urls) - 1}\x02"
        t = re.sub(r"(https://[^\s]+[^\s.,;)])", keep_url, t)
        t = tex_escape(t)
        t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
        t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)
        t = re.sub(r"\^([\w/]+)", r"$^{\1}$", t)
        for a, b in (("µs", "\\textmu s"), ("×", "$\\times$"), ("≤", "$\\le$"), ("≥", "$\\ge$"), ("→", "$\\to$"),
                     ("≈", "$\\approx$"), ("−", "$-$"), ("±", "$\\pm$"), ("“", "``"), ("”", "''"), ("’", "'"), ("–", "--"),
                     ("—", "---"), ("·", "$\\cdot$")):
            t = t.replace(a, b)
        t = re.sub(r"(\d),(\d{3})", r"\1{,}\2", t)
        t = re.sub(r"(\d) " + UNIT.replace("µs", "\\\\textmu s") + r"(?![a-zA-Z])", r"\1~\2", t)
        def put_ref(m):
            r = refs[int(m.group(1))]
            if r.startswith("[#tab:"):
                return "Table~\\ref{tab:" + r[6:-1] + "}"
            if r.startswith("[#fig:"):
                return "Fig.~\\ref{fig:" + r[6:-1] + "}"
            if r.startswith("[@!"):
                return "\\citet{" + r[3:-1] + "}"
            keys = [k.strip()[1:] for k in r[1:-1].split(";")]
            return "\\citep{" + ",".join(keys) + "}"
        t = re.sub(r"\x01(\d+)\x01", put_ref, t)
        t = re.sub(r"\x00(\d+)\x00", lambda m: "\\texttt{" + tex_escape(codes[int(m.group(1))]).replace("--", "-{}-") + "}", t)
        t = re.sub(r"\x02(\d+)\x02", lambda m: "\\url{" + urls[int(m.group(1))] + "}", t)
        out.append(t)
    return "".join(out)

# ------------------------------------------------------------------------------------------- HTML
def section_numbers(body, prefix=""):
    sec = sub = 0
    nums = []
    for b in body:
        if b[0] == "sec":
            sec += 1; sub = 0; nums.append(f"{prefix}{sec}")
        elif b[0] == "sub":
            sub += 1; nums.append(f"{prefix}{sec}.{sub}")
        else:
            nums.append(None)
    return nums

def is_heading(row):
    """A row naming the group of rows below it: its first cell alone holds text, across the table."""
    return bool(row[0]) and all(c == "" for c in row[1:])

def html_table(key, t, tabs, figs):
    wrap = set(t.get("wrap", []))
    left = ' class="l"'
    head = "".join(f'<th{left if t["align"][i] == "l" else ""}>{inline_html(h, tabs, figs)}</th>' for i, h in enumerate(t["header"]))
    rows = []
    for r in t["rows"]:
        if is_heading(r):
            rows.append(f'<tr class="group"><td class="l" colspan="{len(r)}">{inline_html(r[0], tabs, figs)}</td></tr>')
            continue
        cells = []
        for i, c in enumerate(r):
            cls = []
            if t["align"][i] == "l": cls.append("l")
            if i in wrap: cls.append("wrap")
            if c.startswith("**"): cls.append("best")
            attr = f' class="{" ".join(cls)}"' if cls else ""
            cells.append(f"<td{attr}>{inline_html(c, tabs, figs)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<div class="table-wrap" id="tab-{key}"><table>\n<caption><b>Table {tabs[key]}.</b> {inline_html(t["caption"], tabs, figs)}</caption>\n'
            f'<thead><tr>{head}</tr></thead>\n<tbody>\n' + "\n".join(rows) + "\n</tbody></table></div>")

def html_figure(key, caption, num, tabs, figs):
    cap = f'<figcaption><b>Fig. {num}.</b> {inline_html(caption, tabs, figs)}</figcaption>'
    return f'<figure id="fig-{key}">\n  <div class="fig-box chart-wrap">{F.svg(key)}</div>\n  {cap}\n</figure>'

def html_body(body, tables, figures, prefix=""):
    tabs, figs = number_items(body, tables)
    nums = section_numbers(body, prefix)
    out, open_col = [], False
    def col():
        nonlocal open_col
        if not open_col:
            out.append('<div class="col">'); open_col = True
    def end_col():
        nonlocal open_col
        if open_col:
            out.append("</div>"); open_col = False
    for b, n in zip(body, nums):
        if b[0] == "sec":
            end_col(); col(); out.append(f'<h2 id="s{n}"><span class="num">{n}</span>{inline_html(b[1], tabs, figs)}</h2>')
        elif b[0] == "sub":
            col(); out.append(f'<h3 id="s{n}"><span class="num">{n}</span>{inline_html(b[1], tabs, figs)}</h3>')
        elif b[0] == "p":
            col()
            if b[1].startswith("\\["):
                out.append(f'<div class="math-block">{b[1]}</div>')
            else:
                out.append(f"<p>{inline_html(b[1], tabs, figs)}</p>")
        elif b[0] in ("list", "bullets"):
            col(); tag = "ol" if b[0] == "list" else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{inline_html(i, tabs, figs)}</li>" for i in b[1]) + f"</{tag}>")
        elif b[0] == "table":
            end_col(); out.append(html_table(b[1], tables[b[1]], tabs, figs))
        elif b[0] == "figure":
            end_col(); out.append(html_figure(b[1], figures[b[1]], figs[b[1]], tabs, figs))
    end_col()
    return "\n".join(out)

def cited_keys(*parts):
    text = json.dumps(parts)
    return sorted(set(re.findall(r"@!?(\w+)", text)) & set(R.REFS))

def reference_list_html(keys):
    items = []
    for k in sorted(keys, key=lambda k: (R.REFS[k]["authors"][0][0].lower(), R.REFS[k]["year"])):
        a, t, w, doi, url = R.html_entry(k)
        link = (f' <a class="doi" href="https://doi.org/{doi}">doi:{doi}</a>' if doi else
                f' <a class="doi" href="{url}">{url}</a>' if url else "")
        items.append(f'<li id="ref-{k}">{htmlmod.escape(a)} ({R.REFS[k]["year"]}) {htmlmod.escape(t)} <i>{htmlmod.escape(w)}</i>{link}</li>')
    return "<ul class=\"refs\">\n" + "\n".join(items) + "\n</ul>"

def style():
    return open(os.path.join(HERE, "assets", "style.css")).read()

MATHJAX = '''<script>
window.MathJax = { tex: { inlineMath: [["\\\\(", "\\\\)"]], displayMath: [["\\\\[", "\\\\]"]] }, svg: { fontCache: "global" } };
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>'''

def render_html_main():
    abstract = "\n".join(f'<p><b>{h}:</b> {inline_html(t, {}, {})}</p>' for h, t in M.ABSTRACT)
    body = html_body(M.BODY, M.TABLES, M.FIGURES)
    back = "\n".join(f'<h4>{h}</h4><p>{inline_html(t, {}, {})}</p>' for h, t in M.BACK)
    refs = reference_list_html(cited_keys(M.BODY, M.ABSTRACT, M.TABLES, M.FIGURES, M.BACK, S.BODY, S.TABLES))
    title = inline_html(M.TITLE, {}, {}).replace("pgphase:", '<span class="sys">pgphase</span>:', 1)
    authors = " and ".join(f"{name}<sup>1,*</sup>" for name, _ in M.AUTHORS)
    return f'''<title>pgphase paper</title>
<style>
{style()}
</style>
<main>
<header class="col masthead">
  <div class="journal">Bioinformatics · Original Paper · Genome analysis <span class="draft">Draft</span></div>
  <h1>{title}</h1>
  <div class="authors">{authors}</div>
  <div class="affil"><sup>1</sup>{M.AFFILIATION}</div>
  <div class="affil">*To whom correspondence should be addressed: {", ".join(email for _, email in M.AUTHORS)}</div>
</header>
<section class="col abstract" id="abstract">
<h4>Abstract</h4>
{abstract}
</section>
{body}
<section class="col back">
{back}
</section>
<section class="col" id="references">
<h2>References</h2>
{refs}
</section>
</main>
{MATHJAX}
'''

def render_html_supplement():
    body = html_body(S.BODY, S.TABLES, {}, prefix="S")
    return f'''<title>pgphase supplement</title>
<style>
{style()}
</style>
<main>
<header class="col masthead">
  <div class="journal">Bioinformatics · Supplementary material <span class="draft">Draft</span></div>
  <h1>{inline_html(S.TITLE, {}, {})}</h1>
  <div class="authors">{" and ".join(name for name, _ in M.AUTHORS)}</div>
</header>
{body}
</main>
{MATHJAX}
'''

# ------------------------------------------------------------------------------------------- LaTeX
def tex_table(key, t, wide):
    wrap = set(t.get("wrap", []))
    cols = []
    for i, a in enumerate(t["align"]):
        if i in wrap:
            cols.append(f">{{\\raggedright\\arraybackslash}}p{{{'5.2cm' if wide else '3.1cm'}}}")
        else:
            cols.append(a)
    env = "table*" if wide else "table"
    fit = ("\\resizebox{\\textwidth}{!}{", "}") if t.get("fit") and wide else (("\\resizebox{\\linewidth}{!}{", "}") if t.get("fit") else ("", ""))
    body = "\n".join(
        (f"\\multicolumn{{{len(r)}}}{{@{{}}l}}{{\\emph{{{inline_tex(r[0])}}}}} \\\\" if is_heading(r)
         else " & ".join(inline_tex(c) for c in r) + " \\\\")
        for r in t["rows"])
    head = " & ".join(inline_tex(h) for h in t["header"]) + " \\\\"
    return (f"\\begin{{{env}}}[t]\n\\caption{{{inline_tex(t['caption'])}}}\\label{{tab:{key}}}\n\\centering\\footnotesize\n"
            f"\\setlength{{\\tabcolsep}}{{4pt}}\n{fit[0]}\\begin{{tabular}}{{@{{}}{''.join(cols)}@{{}}}}\n\\toprule\n{head}\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}{fit[1]}\n\\end{{{env}}}\n")

def tex_figure(key, caption):
    env = "figure*" if F.WIDE.get(key) else "figure"
    width = "\\textwidth" if F.WIDE.get(key) else "\\linewidth"
    return (f"\\begin{{{env}}}[t]\n\\centering\n\\resizebox{{{width}}}{{!}}{{\\begin{{tikzpicture}}\n{F.tikz(key)}\n\\end{{tikzpicture}}}}\n"
            f"\\caption{{{inline_tex(caption)}}}\\label{{fig:{key}}}\n\\end{{{env}}}\n")

TEX_PREAMBLE = r"""\usepackage[a4paper,margin=1.8cm,columnsep=0.7cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage{graphicx}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{xcolor}
\usepackage{tikz}
\usepackage[round]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{url}
"""

def tex_body(body, tables, figures, wide_tables=(), parts=None):
    """The LaTeX of `body`; with `parts`, each table and figure goes there by its file name, and the body
    inputs it."""
    out = []

    def place(kind, key, text):
        if parts is None:
            out.append(text)
        else:
            parts[f"{kind}/{key}.tex"] = text
            out.append(f"\\input{{{kind}/{key}}}\n")

    for b in body:
        if b[0] == "sec":
            out.append(f"\\section{{{inline_tex(b[1])}}}\n")
        elif b[0] == "sub":
            out.append(f"\\subsection{{{inline_tex(b[1])}}}\n")
        elif b[0] == "p":
            if b[1].startswith("\\["):
                out.append("\\begin{equation}\n" + b[1][2:-2].strip() + "\n\\end{equation}\n")
            else:
                out.append(inline_tex(b[1]) + "\n")
        elif b[0] in ("list", "bullets"):
            env = "enumerate" if b[0] == "list" else "itemize"
            out.append(f"\\begin{{{env}}}\n" + "\n".join(f"  \\item {inline_tex(i)}" for i in b[1]) + f"\n\\end{{{env}}}\n")
        elif b[0] == "table":
            t = tables[b[1]]
            place("tables", b[1], tex_longtable(b[1], t) if t.get("long") else tex_table(b[1], t, b[1] in wide_tables))
        elif b[0] == "figure":
            place("figures", b[1], tex_figure(b[1], figures[b[1]]))
    return "\n".join(out)

def tex_longtable(key, t):
    wrap = set(t.get("wrap", []))
    cols = "".join(f">{{\\raggedright\\arraybackslash}}p{{5.0cm}}" if i in wrap else a for i, a in enumerate(t["align"]))
    head = " & ".join(inline_tex(h) for h in t["header"]) + " \\\\"
    rows = "\n".join((f"\\multicolumn{{{len(r)}}}{{@{{}}l}}{{\\emph{{{inline_tex(r[0])}}}}} \\\\" if is_heading(r)
                      else " & ".join(inline_tex(c) for c in r) + " \\\\") for r in t["rows"])
    return (f"{{\\footnotesize\n\\setlength{{\\LTcapwidth}}{{\\linewidth}}\n\\begin{{longtable}}{{@{{}}{cols}@{{}}}}\n"
            f"\\caption{{{inline_tex(t['caption'])}}}\\label{{tab:{key}}}\\\\\n\\toprule\n{head}\n\\midrule\n\\endfirsthead\n"
            f"\\caption[]{{(continued)}}\\\\\n\\toprule\n{head}\n\\midrule\n\\endhead\n{rows}\n\\bottomrule\n\\end{{longtable}}\n}}\n")

# The journal's class: OUP's General Template (oup-authoring-template, January 2026), which Bioinformatics
# recommends for submission, with the options and author-year style its sheet gives for the journal. The class
# loads natbib, hyperref, tikz, xcolor, url, amsmath and amssymb itself; it and the style sit beside the paper.
JOURNAL_CLASS = "oup-authoring-template"
JOURNAL_STYLE = "oup-abbrvnat"
GENERATED = "% Generated by tools/paper/render.py from main_text.py; edit the source there, not this file.\n"


def render_project():
    """The main paper as an Overleaf project, each file by its path in paper/: main.tex with the journal's front
    matter, one file a section, and the figures and tables the sections input."""
    parts = {}
    sections, current = [], None
    for block in M.BODY:
        if block[0] == "sec":
            current = (block[1].lower().replace(" ", "-"), [])
            sections.append(current)
        current[1].append(block)
    names = []
    for name, blocks in sections:
        parts[f"sections/{name}.tex"] = tex_body(blocks, M.TABLES, M.FIGURES, wide_tables=M.WIDE_TABLES, parts=parts)
        names.append(name)
    parts["sections/abstract.tex"] = (
        "\\abstract{" + "\\\\\n".join(f"\\textbf{{{h}:}} {inline_tex(t)}" for h, t in M.ABSTRACT) + "}\n")
    parts["sections/back-matter.tex"] = "\n".join(f"\\section*{{{h}}}\n{inline_tex(t)}\n" for h, t in M.BACK)
    authors = "\n".join(f"\\author[1,$\\ast$]{{{name}}}" for name, _ in M.AUTHORS)
    emails = ", ".join(f"\\href{{mailto:{email}}}{{{email}}}" for _, email in M.AUTHORS)
    organisation, *place = [part.strip() for part in M.AFFILIATION.split(",")]
    inputs = "\n".join(f"\\input{{sections/{name}}}" for name in names)
    corresp = "Corresponding author" if len(M.AUTHORS) == 1 else "Corresponding authors"
    parts["main.tex"] = f"""\\documentclass[numsec,webpdf,modern,large]{{{JOURNAL_CLASS}}}
\\usepackage{{booktabs}}
{F.TIKZ_STYLES}
\\begin{{document}}
\\journaltitle{{Bioinformatics}}
\\DOI{{DOI added during production}}
\\copyrightyear{{2026}}
\\pubyear{{2026}}
\\access{{Published: Date added during production}}
\\appnotes{{Original Paper}}
\\firstpage{{1}}
\\title[pgphase]{{{inline_tex(M.TITLE)}}}
{authors}
\\address[1]{{\\orgname{{{tex_escape(organisation)}}}, \\orgaddress{{\\state{{{tex_escape(", ".join(place[:-1]))}}}, \\country{{{tex_escape(place[-1])}}}}}}}
\\corresp[$\\ast$]{{{corresp}. {emails}}}
\\received{{Date}}{{0}}{{Year}}
\\revised{{Date}}{{0}}{{Year}}
\\accepted{{Date}}{{0}}{{Year}}
\\input{{sections/abstract}}
\\keywords{{{", ".join(M.KEYWORDS)}}}
\\maketitle
{inputs}
\\input{{sections/back-matter}}
\\bibliographystyle{{{JOURNAL_STYLE}}}
\\bibliography{{references}}
\\end{{document}}
"""
    out = {path: GENERATED + text for path, text in parts.items()}
    out["references.bib"] = R.bibtex()
    return out


def render_tex_supplement():
    # The supplement's tables stay where their section introduces them.
    body = tex_body(S.BODY, S.TABLES, {}).replace("\\begin{table}[t]", "\\begin{table}[h]")
    return f"""% Generated by tools/paper/render.py from supplement.py; edit the source there, not this file.
\\documentclass[10pt]{{article}}
{TEX_PREAMBLE}
\\renewcommand{{\\thesection}}{{S\\arabic{{section}}}}
\\renewcommand{{\\thetable}}{{S\\arabic{{table}}}}
\\begin{{document}}
{{\\Large\\bfseries {inline_tex(S.TITLE)}\\par}}\\medskip
{{\\large {" and ".join(name for name, _ in M.AUTHORS)}}}\\par\\bigskip
{body}
\\end{{document}}
"""

def build_pdfs():
    """Compiles paper/ as Overleaf would, in a copy under build/paper/latex so that paper/ keeps no build output,
    and leaves build/paper/pgphase.pdf and supplement.pdf; says so and returns if pdflatex is missing."""
    import shutil
    import subprocess

    if not (shutil.which("pdflatex") and shutil.which("bibtex")):
        print("pdflatex or bibtex not found: no PDFs built")
        return
    work = os.path.join(HTML_OUT, "latex")
    shutil.rmtree(work, ignore_errors=True)
    shutil.copytree(PAPER, work)
    for name, target in (("main", "pgphase"), ("supplement", "supplement")):
        def run(*command):
            return subprocess.run(command, cwd=work, capture_output=True, text=True)

        run("pdflatex", "-interaction=nonstopmode", f"{name}.tex")
        run("bibtex", name)
        run("pdflatex", "-interaction=nonstopmode", f"{name}.tex")
        last = run("pdflatex", "-interaction=nonstopmode", f"{name}.tex")
        log = open(os.path.join(work, f"{name}.log"), errors="replace").read()
        errors = [line for line in log.splitlines() if line.startswith("!")]
        unresolved = re.findall(r"(?:Citation|Reference) `[^']+' on page \d+ undefined", log)
        pdf = os.path.join(work, f"{name}.pdf")
        if last.returncode != 0 or errors or not os.path.exists(pdf):
            print(f"{name}.tex failed: {errors[:3] or 'no PDF'} (log: {work}/{name}.log)")
            continue
        shutil.copy(pdf, os.path.join(HTML_OUT, f"{target}.pdf"))
        pages = re.search(r"Output written on .*?\((\d+) pages", last.stdout)
        note = f", {len(unresolved)} unresolved references" if unresolved else ""
        print(f"build/paper/{target}.pdf: {pages.group(1) if pages else '?'} pages{note}")


if __name__ == "__main__":
    os.makedirs(HTML_OUT, exist_ok=True)
    open(os.path.join(HTML_OUT, "paper.html"), "w").write(render_html_main())
    open(os.path.join(HTML_OUT, "supplement.html"), "w").write(render_html_supplement())
    for kind in ("sections", "figures", "tables"):
        directory = os.path.join(PAPER, kind)
        os.makedirs(directory, exist_ok=True)
        # Only files this script wrote, so that a section renamed here leaves no stale copy behind.
        for name in os.listdir(directory):
            path = os.path.join(directory, name)
            if name.endswith(".tex") and open(path).readline() == GENERATED:
                os.remove(path)
    for path, text in render_project().items():
        open(os.path.join(PAPER, path), "w").write(text)
    open(os.path.join(PAPER, "supplement.tex"), "w").write(render_tex_supplement())
    print("rendered paper/ and build/paper/*.html")
    # The project as one zip, for Overleaf's New Project > Upload Project: its files at the zip's top level.
    import zipfile

    archive = os.path.join(HTML_OUT, "pgphase-overleaf.zip")
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as out:
        for folder, _, files in os.walk(PAPER):
            for name in files:
                if name != ".gitignore":
                    path = os.path.join(folder, name)
                    out.write(path, os.path.relpath(path, PAPER))
    print("build/paper/pgphase-overleaf.zip: the project, for Overleaf's Upload Project")
    build_pdfs()
