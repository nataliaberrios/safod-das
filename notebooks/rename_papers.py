"""
Rename paper PDFs to author_year_journal_topic.pdf.

    python3 rename_papers.py somefile.pdf              # preview one
    python3 rename_papers.py ~/planpapers              # preview a directory
    python3 rename_papers.py ~/planpapers --apply      # do it
    python3 rename_papers.py file.pdf --into 02_tides --apply

Preview is the default; nothing moves without --apply.

HOW IT DECIDES, in order, because filenames lie and first pages are messy:

  1. DOI FROM THE PDF TEXT -> CrossRef. Authoritative: publisher-supplied author,
     year, container title and article title. This is the only route that cannot be
     fooled by a journal's cover page or by a first page that happens to start with
     the tail of the preceding article -- which is exactly what made
     science.285.5428.718.pdf look like a condensed-matter paper when it is
     Nadeau & McEvilly 1999.

  2. THE WILEY FILENAME PATTERN, "Journal - Year - Author - Title...", which most
     of this collection already follows.

  3. FIRST-PAGE HEURISTICS -- a four-digit year, a plausible surname. Weakest;
     always shown for review rather than applied silently.

Every applied rename is appended to MANIFEST.csv beside the files, so the
reorganisation stays reversible.

Needs pdftotext (ml system poppler). Network is optional -- without it, step 1 is
skipped and the other two still work.
"""
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import unicodedata
import urllib.error
import urllib.request

DOI_RE = re.compile(r'\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+\b')
JOURNALS = [
    # ORDER MATTERS: specific multi-word titles first. "Progress in Earth and
    # Planetary Science" and "Earth and Planetary Science Letters" both contain
    # the word "science", so a bare \bscience\b rule placed above them wins and
    # mislabels PEPS as Science.
    (r'progress in earth and planetary', 'peps'),
    (r'earth and planetary science letters', 'epsl'),
    (r'annual review', 'annurev'),
    (r'science advances', 'sciadv'),
    (r'nature geoscience', 'natgeo'),
    (r'nature communications', 'natcomm'),
    (r'geophysical research letters', 'grl'),
    (r'journal of geophysical research', 'jgr'),
    (r'bulletin of the seismological society', 'bssa'),
    (r'seismological research letters', 'srl'),
    (r'geophysical journal international', 'gji'),
    (r'the leading edge', 'tle'),
    (r'hydrology and earth system', 'hess'),
    (r'near surface geophysics', 'nsg'),
    (r'earth and space science', 'ess'),
    (r'earth, planets and space', 'eps'),
    (r'\bgeophysics\b', 'geophysics'),
    (r'tectonophysics', 'tectonophysics'),
    (r'seismica', 'seismica'),
    (r'solid earth', 'se'),
    (r'\bnature\b', 'nature'),
    (r'\bscience\b', 'science'),
]
STOP = {'a', 'an', 'the', 'of', 'in', 'on', 'at', 'to', 'for', 'from', 'and',
        'with', 'by', 'its', 'their', 'using', 'via', 'into', 'is', 'are', 'as',
        'new', 'towards', 'toward', 'case', 'study', 'evidence'}


def ascii_fold(s):
    """Bürgmann -> Burgmann, Lengliné -> Lengline. Also fixes the NFC/NFD split
    that makes accented filenames compare unequal to themselves."""
    s = unicodedata.normalize('NFKD', s)
    return ''.join(c for c in s if not unicodedata.combining(c))


def slug(text, maxwords=6):
    text = ascii_fold(text).lower()
    text = re.sub(r'[^a-z0-9\s-]', ' ', text)
    words = [w for w in text.split() if w and w not in STOP and len(w) > 2]
    return '-'.join(words[:maxwords]) or 'untitled'


def journal_code(name):
    n = ascii_fold(name or '').lower()
    for pat, code in JOURNALS:
        if re.search(pat, n):
            return code
    w = [x for x in re.split(r'[^a-z]+', n) if x and x not in STOP]
    return ''.join(x[0] for x in w[:4]) or 'journal'


def pdf_text(path, pages=2):
    try:
        return subprocess.run(['pdftotext', '-f', '1', '-l', str(pages),
                               path, '-'], capture_output=True, text=True,
                              timeout=60).stdout
    except Exception:
        return ''


def from_crossref(doi):
    url = f'https://api.crossref.org/works/{urllib.parse.quote(doi)}'
    req = urllib.request.Request(url, headers={
        'User-Agent': 'paper-renamer/1.0 (mailto:nberrios@stanford.edu)'})
    try:
        with urllib.request.urlopen(req, timeout=25) as r:
            m = json.load(r)['message']
    except Exception:
        return None
    auth = m.get('author') or []
    surname = auth[0].get('family') if auth and auth[0].get('family') else None
    if not surname:
        return None
    date = (m.get('issued', {}).get('date-parts') or [[None]])[0]
    year = date[0] if date and date[0] else None
    cont = (m.get('container-title') or [''])[0]
    title = (m.get('title') or [''])[0]
    if not (year and title):
        return None
    # "De Fazio" -> "defazio", "van der Elst" -> "vanderelst": surnames with
    # spaces or hyphens would otherwise break the underscore-delimited scheme
    surname = re.sub(r'[^a-z0-9]', '', ascii_fold(surname).lower())
    return (surname, str(year), journal_code(cont), slug(title), 'crossref')


def from_filename(stem):
    """Wiley export style: 'Journal - YYYY - Surname - Title words'."""
    parts = [p.strip() for p in stem.split(' - ')]
    if len(parts) < 4:
        return None
    journal, year, surname = parts[0], parts[1], parts[2]
    if not re.fullmatch(r'(19|20)\d{2}', year):
        return None
    title = ' - '.join(parts[3:])
    surname = re.sub(r'[^a-z0-9]', '', ascii_fold(surname).lower())
    return (surname, year, journal_code(journal), slug(title), 'filename')


def from_firstpage(text):
    yrs = re.findall(r'\b(19[5-9]\d|20[0-2]\d)\b', text[:3000])
    if not yrs:
        return None
    lines = [l.strip() for l in text.splitlines() if len(l.strip()) > 12][:25]
    title = lines[0] if lines else 'untitled'
    surname = 'unknown'
    for l in lines[:20]:
        m = re.search(r'\b([A-Z][a-z]{2,})\s*[,\d*†‡]', l)
        if m and m.group(1).lower() not in STOP:
            surname = m.group(1).lower(); break
    jr = None
    for pat, code in JOURNALS:
        if re.search(pat, ascii_fold(text[:3000]).lower()):
            jr = code; break
    return (surname, max(yrs, key=yrs.count), jr or 'journal',
            slug(title), 'firstpage')


def propose(path, offline=False):
    stem = os.path.splitext(os.path.basename(path))[0]
    text = pdf_text(path)
    if not offline:
        for d in DOI_RE.findall(text)[:4]:
            got = from_crossref(d.rstrip('.,;'))
            if got:
                return got
    return from_filename(stem) or from_firstpage(text) or \
        ('unknown', '0000', 'journal', slug(stem), 'fallback')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('target')
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--into', default=None, help='subfolder to move into')
    ap.add_argument('--offline', action='store_true', help='skip CrossRef')
    a = ap.parse_args()

    files = ([a.target] if a.target.lower().endswith('.pdf')
             else sorted(os.path.join(r, f)
                         for r, _, fs in os.walk(a.target) for f in fs
                         if f.lower().endswith('.pdf')))
    if not files:
        print('no PDFs found'); return

    rows = []
    print(f'{"src":>6}  {"proposed name":<74} from')
    for p in files:
        surname, year, jr, topic, src = propose(p, a.offline)
        new = f'{surname}_{year}_{jr}_{topic}.pdf'
        dest_dir = os.path.join(os.path.dirname(p), a.into) if a.into \
            else os.path.dirname(p)
        flag = ' <-- CHECK' if src in ('firstpage', 'fallback') else ''
        print(f'{src:>6}  {new:<74}{flag}')
        rows.append((p, os.path.join(dest_dir, new), src))

    if not a.apply:
        print(f'\n{len(rows)} file(s). Preview only -- add --apply to rename.')
        print('Anything marked CHECK was guessed from the first page; verify it.')
        return

    root = a.target if os.path.isdir(a.target) else os.path.dirname(a.target)
    man = os.path.join(root, 'MANIFEST.csv')
    exists = os.path.exists(man)
    done = 0
    with open(man, 'a', newline='') as fh:
        w = csv.writer(fh)
        if not exists:
            w.writerow(['folder', 'new_name', 'old_name'])
        for old, new, src in rows:
            if os.path.abspath(old) == os.path.abspath(new):
                continue
            os.makedirs(os.path.dirname(new), exist_ok=True)
            if os.path.exists(new):
                print(f'SKIP (exists): {os.path.basename(new)}'); continue
            shutil.move(old, new)
            w.writerow([os.path.basename(os.path.dirname(new)),
                        os.path.basename(new), os.path.basename(old)])
            done += 1
    print(f'\n{done} renamed; recorded in {man}')


if __name__ == '__main__':
    main()
