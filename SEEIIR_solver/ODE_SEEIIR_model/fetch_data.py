"""Fetch 2026 Bundibugyo Ebola (DRC / Uganda) surveillance data into the CSV
template used by seeiir_ebola.py.

Sources, in order of preference
  1. GitHub tracker CSV  (plain HTTP download, no browser needed)
       webbegole/ebola-bundibugyo-tracker-2026  data/country_breakdown.csv
     A third-party compilation of DRC INSP situation reports. Each row cites its
     primary source. Dates follow the tracker's convention (publication date,
     typically report date + 1 day for DRC).
  2. ECDC outbreak page   (cross-check / fallback; one "as of" snapshot per run)
       fetched with requests, or with Selenium when --selenium is given or the
       plain request is blocked.

Output columns (same as ebola_drc_uganda_template.csv):
  date,country,cum_confirmed_cases,cum_confirmed_deaths,cum_recovered,quality,source_note

Examples
  python fetch_data.py                       # tracker -> ebola_drc_uganda.csv
  python fetch_data.py --crosscheck          # also compare latest row with ECDC
  python fetch_data.py --crosscheck --selenium
  python fetch_data.py --source ecdc --append  # append an ECDC snapshot to the CSV

Dependencies: pandas, requests; optional: beautifulsoup4, selenium (Selenium 4.6+
downloads its own driver but needs Chrome/Chromium installed).
"""
import argparse
import datetime as dt
import io
import re
import sys
from pathlib import Path

import pandas as pd

TRACKER_URL = ("https://raw.githubusercontent.com/webbegole/"
               "ebola-bundibugyo-tracker-2026/main/data/country_breakdown.csv")
ECDC_URL = "https://www.ecdc.europa.eu/en/ebola-outbreak-democratic-republic-congo-and-uganda"
COUNTRIES = ("DRC", "Uganda")
OUT_COLS = ["date", "country", "cum_confirmed_cases", "cum_confirmed_deaths",
            "cum_recovered", "quality", "source_note"]
HEADERS = {"User-Agent": "Mozilla/5.0 (research script; ebola-seeiir)"}


# --------------------------------------------------------------------------
# 1. GitHub tracker (primary)
# --------------------------------------------------------------------------
def _short_note(text):
    """Condense the tracker's very long source field into a short note."""
    text = str(text)
    m = re.search(r"SitRep N[°o]\s*(\d+)", text)
    if m:
        return f"INSP SitRep {m.group(1)} via tracker"
    text = re.sub(r"\s+", " ", re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text))
    return text[:70].replace(",", ";")


def fetch_tracker(url=TRACKER_URL, timeout=60):
    import requests
    r = requests.get(url, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    raw = r.text
    df = pd.read_csv(io.StringIO(raw), parse_dates=["date"])
    df = df[df["country"].isin(COUNTRIES)].copy()
    out = pd.DataFrame({
        "date": df["date"],
        "country": df["country"],
        "cum_confirmed_cases": df["confirmed"],
        "cum_confirmed_deaths": df["confirmed_deaths"],
        "cum_recovered": df["recovered"],
        "quality": "compiled",
        "source_note": df["primary_source"].map(_short_note),
    })
    return out, raw


# --------------------------------------------------------------------------
# 2. ECDC page (cross-check / fallback)
# --------------------------------------------------------------------------
_WORDS = {w: i for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen "
    "fourteen fifteen sixteen seventeen eighteen nineteen twenty".split())}
_SP = r"[\s\u00a0\u202f]"                      # ECDC uses spaces as thousands separator
_NUMTOK = rf"(\d[\d\u00a0\u202f ]*|[a-z]+)"    # digits (with spaces) or a number word


def _to_int(tok):
    tok = tok.strip().lower()
    if tok in _WORDS:
        return _WORDS[tok]
    return int(re.sub(r"[\s\u00a0\u202f,]", "", tok))


def _parse_day(text, year):
    return dt.datetime.strptime(f"{text.strip()} {year}", "%d %B %Y").date()


def parse_ecdc_text(text):
    """Extract DRC and Uganda cumulative confirmed cases/deaths from ECDC prose.

    Returns a list of dict rows (possibly empty). Written against the sentence
    patterns ECDC has used so far, e.g.
      'the DRC Ministry of Health reported a total of 1 118 confirmed cases,
       including 291 confirmed related deaths ... (as of 23 June)'
      'As of 24 June, Uganda had reported a total of 20 confirmed cases,
       including two deaths.'
    If ECDC changes its wording this returns [] and the caller warns.
    """
    text = re.sub(r"\s+", " ", text.replace("\u00a0", " ").replace("\u202f", " "))
    ym = re.search(r"[Aa]s of \d{1,2} [A-Z][a-z]+ (\d{4})", text)
    year = int(ym.group(1)) if ym else dt.date.today().year
    rows = []

    m = re.search(r"Ministry of Health (?:published updated figures, )?reported a total of "
                  r"(\d[\d ]*) confirmed cases, including (\d[\d ]*) confirmed(?: related)? deaths"
                  r".{0,200}?\(as of (\d{1,2} [A-Z][a-z]+)\)", text)
    if m:
        rows.append(dict(date=_parse_day(m.group(3), year), country="DRC",
                         cum_confirmed_cases=_to_int(m.group(1)),
                         cum_confirmed_deaths=_to_int(m.group(2))))

    m = re.search(r"[Aa]s of (\d{1,2} [A-Z][a-z]+),? Uganda (?:had )?reported a total of "
                  rf"{_NUMTOK} confirmed cases,? including {_NUMTOK} deaths", text)
    if m:
        rows.append(dict(date=_parse_day(m.group(1), year), country="Uganda",
                         cum_confirmed_cases=_to_int(m.group(2)),
                         cum_confirmed_deaths=_to_int(m.group(3))))
    return rows


def _page_text_requests(url, timeout=60):
    import requests
    r = requests.get(url, headers=HEADERS, timeout=timeout)
    r.raise_for_status()
    try:
        from bs4 import BeautifulSoup
        return BeautifulSoup(r.text, "html.parser").get_text(" ")
    except ImportError:
        return re.sub(r"<[^>]+>", " ", r.text)


def _page_text_selenium(url, timeout=60):
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait

    opts = Options()
    opts.add_argument("--headless=new")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument(f"--user-agent={HEADERS['User-Agent']}")
    driver = webdriver.Chrome(options=opts)
    try:
        driver.set_page_load_timeout(timeout)
        driver.get(url)
        # Dismiss a cookie banner if one is present (best effort, selector may change)
        for sel in ("button#onetrust-accept-btn-handler", "button.cck-actions-button"):
            try:
                driver.find_element(By.CSS_SELECTOR, sel).click()
                break
            except Exception:
                pass
        WebDriverWait(driver, 30).until(EC.presence_of_element_located((By.TAG_NAME, "main")))
        return driver.find_element(By.TAG_NAME, "body").text
    finally:
        driver.quit()


def fetch_ecdc(url=ECDC_URL, use_selenium=False):
    text = None
    if not use_selenium:
        try:
            text = _page_text_requests(url)
        except Exception as e:
            print(f"[ecdc] plain request failed ({e}); trying Selenium", file=sys.stderr)
    if text is None:
        text = _page_text_selenium(url)
    rows = parse_ecdc_text(text)
    if not rows:
        print("[ecdc] WARNING: could not find case counts - page wording may have changed",
              file=sys.stderr)
    out = pd.DataFrame(rows)
    if len(out):
        out["date"] = pd.to_datetime(out["date"])
        out["cum_recovered"] = pd.NA
        out["quality"] = "snapshot"
        out["source_note"] = "ECDC outbreak page"
        out = out[OUT_COLS]
    return out


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
def validate(df, enforce_monotone=False):
    df = (df.dropna(subset=["cum_confirmed_cases"])
            .sort_values(["country", "date"])
            .drop_duplicates(["country", "date"], keep="last")
            .reset_index(drop=True))
    for c, g in df.groupby("country"):
        for col in ("cum_confirmed_cases", "cum_confirmed_deaths", "cum_recovered"):
            s = g[col].dropna()
            dips = int((s.diff() < 0).sum())
            if dips:
                print(f"[validate] {c}: {dips} decrease(s) in cumulative {col} "
                      f"(data revisions). Fit is on log scale so this is usually minor.",
                      file=sys.stderr)
    if enforce_monotone:
        for col in ("cum_confirmed_cases", "cum_confirmed_deaths", "cum_recovered"):
            df[col] = df.groupby("country")[col].cummax()
    return df[OUT_COLS]


def crosscheck(tracker_df, ecdc_df, tol=0.05):
    """Compare ECDC snapshot with the nearest-in-time tracker row."""
    for _, e in ecdc_df.iterrows():
        g = tracker_df[tracker_df["country"] == e["country"]]
        if g.empty:
            continue
        i = (g["date"] - e["date"]).abs().idxmin()
        t = g.loc[i]
        gap = abs((t["date"] - e["date"]).days)
        rel = abs(t["cum_confirmed_cases"] - e["cum_confirmed_cases"]) / max(e["cum_confirmed_cases"], 1)
        flag = "OK" if rel <= tol else "MISMATCH"
        print(f"[crosscheck] {e['country']}: ECDC {int(e['cum_confirmed_cases'])} on {e['date'].date()} "
              f"vs tracker {int(t['cum_confirmed_cases'])} on {t['date'].date()} "
              f"({gap} d apart, {rel:.1%}) -> {flag}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="ebola_drc_uganda.csv")
    ap.add_argument("--source", choices=["tracker", "ecdc"], default="tracker")
    ap.add_argument("--crosscheck", action="store_true", help="also fetch ECDC and compare")
    ap.add_argument("--selenium", action="store_true", help="use headless Chrome for the ECDC page")
    ap.add_argument("--append", action="store_true", help="append to / update an existing --out file")
    ap.add_argument("--enforce-monotone", action="store_true")
    ap.add_argument("--raw-dir", default="raw", help="where to archive the raw download")
    args = ap.parse_args(argv)

    if args.source == "tracker":
        df, raw = fetch_tracker()
        Path(args.raw_dir).mkdir(exist_ok=True)
        stamp = dt.datetime.now().strftime("%Y%m%d_%H%M")
        (Path(args.raw_dir) / f"country_breakdown_{stamp}.csv").write_text(raw, encoding="utf-8")
        if args.crosscheck:
            try:
                crosscheck(df, fetch_ecdc(use_selenium=args.selenium))
            except Exception as e:
                print(f"[crosscheck] skipped: {e}", file=sys.stderr)
    else:
        df = fetch_ecdc(use_selenium=args.selenium)

    if args.append and Path(args.out).exists():
        old = pd.read_csv(args.out, parse_dates=["date"])
        df = pd.concat([old, df], ignore_index=True)   # validate() keeps the newest duplicate

    df = validate(df, enforce_monotone=args.enforce_monotone)
    df.to_csv(args.out, index=False, date_format="%Y-%m-%d")

    print(f"Wrote {args.out}: {len(df)} rows")
    for c, g in df.groupby("country"):
        last = g.iloc[-1]
        print(f"  {c:7s} {g['date'].min().date()} -> {last['date'].date()}  "
              f"latest: {int(last['cum_confirmed_cases'])} confirmed, "
              f"{'' if pd.isna(last['cum_confirmed_deaths']) else int(last['cum_confirmed_deaths'])} deaths")
    return df


if __name__ == "__main__":
    main()
