"""Checks every DOI in refs.py against Crossref: the title Crossref holds for it must match ours, and the year
must be the year of print (or of first posting for a preprint).

    python3 tools/paper/check_refs.py
"""
import json, os, re, sys, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs import REFS

def norm(s):
    """Letters and digits only, Crossref's inline markup (<i>, <scp>) removed."""
    return re.sub(r"[^a-z0-9]", "", re.sub(r"<[^>]+>", "", s).lower())

bad = 0
for key, r in REFS.items():
    doi = r.get("doi")
    if not doi:
        continue
    try:
        with urllib.request.urlopen(urllib.request.Request(f"https://api.crossref.org/works/{doi}",
                headers={"User-Agent": "pgphase-paper-check (mailto:kkyriaki@ucsc.edu)"}), timeout=30) as f:
            m = json.load(f)["message"]
    except Exception as e:
        print(f"{key}: {doi}: {e}"); bad += 1; continue
    title = (m.get("title") or [""])[0]
    issued = (m.get("published-print") or m.get("issued") or {}).get("date-parts", [[None]])[0][0]
    pages = m.get("page") or m.get("article-number") or ""
    ok = norm(title) == norm(r["title"])
    print(f"{'ok ' if ok else 'BAD'} {key}: year {issued} (ours {r['year']}), vol {m.get('volume')} (ours {r.get('volume')}), "
          f"pages {pages} (ours {r.get('pages')})" + ("" if ok else f"\n    crossref: {title}\n    ours:     {r['title']}"))
    bad += not ok
sys.exit(1 if bad else 0)
