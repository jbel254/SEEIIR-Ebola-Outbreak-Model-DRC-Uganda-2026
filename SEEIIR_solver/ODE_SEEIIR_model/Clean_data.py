"""Clean the raw Ebola surveillance CSV before it feeds the SEEIIR fit.

Fixes three issues in the raw data:
  1. Missing calendar dates -> insert a row and forward-fill cumulative columns.
  2. Missing values inside present rows -> forward-fill.
  3. Duplicate dates (if any) -> keep the last, assume it is the revised value.

Cumulative counts are monotonic, so forward-fill is correct.
Interpolation would invent growth that was never reported.

Also handles:
  * UTF-8 BOM in the header (common from Excel / Windows tools).
  * Column names with surrounding whitespace or inconsistent case.
  * Unparseable dates (reported explicitly rather than silently dropped).
"""
import sys

import pandas as pd


def _load_and_normalize(input_path):
    """Read the CSV, strip BOM, normalize column names, parse dates."""
    # utf-8-sig strips a BOM if present; harmless otherwise.
    df = pd.read_csv(input_path, encoding="utf-8-sig")

    # Normalize column names: strip whitespace, lowercase.
    df.columns = df.columns.str.strip().str.lower()

    # Fail loudly with a helpful message if 'date' still isn't there.
    if "date" not in df.columns:
        raise ValueError(
            f"'date' column not found in {input_path}. "
            f"Columns present: {df.columns.tolist()}"
        )

    # Parse dates explicitly now that we know the column exists.
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    n_bad = int(df["date"].isna().sum())
    if n_bad:
        bad_rows = df[df["date"].isna()].index.tolist()
        raise ValueError(
            f"{n_bad} row(s) had unparseable dates in {input_path}: "
            f"row indices {bad_rows[:10]}"
        )

    # Ensure the columns the model needs are numeric.
    for col in ["cum_confirmed_cases", "cum_confirmed_deaths", "cum_recovered"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def clean(input_path, output_path, country=None):
    df = _load_and_normalize(input_path)

    # Optional: restrict to a single country first so per-country
    # resampling doesn't smear one country's gaps with another's rows.
    if country is not None:
        df = df[df["country"] == country].copy()

    # Only keep the columns we need for the cleaning step.
    keep = ["date", "country", "cum_confirmed_cases",
            "cum_confirmed_deaths", "cum_recovered"]

    # Some raw files may be missing the optional columns entirely.
    for col in keep:
        if col not in df.columns:
            df[col] = pd.NA

    out = []
    for ctry, g in df.groupby("country"):
        g = g[keep].copy()
        g = g.sort_values("date")

        # 3. Duplicate dates: keep the latest revision.
        g = g.drop_duplicates(subset=["date"], keep="last")

        # 1. Insert missing calendar dates inside the observed range.
        g = g.set_index("date")
        full_range = pd.date_range(g.index.min(), g.index.max(), freq="D")
        g = g.reindex(full_range)

        # country column becomes NaN on inserted rows -> restore it.
        g["country"] = ctry

        # 2. Forward-fill cumulative columns.
        #    Leading NaNs (before the first report) stay NaN and are dropped.
        cum_cols = ["cum_confirmed_cases", "cum_confirmed_deaths", "cum_recovered"]
        g[cum_cols] = g[cum_cols].ffill()

        # Drop rows still missing the column the model absolutely needs.
        # (cum_recovered may legitimately stay NaN for long stretches.)
        g = g.dropna(subset=["cum_confirmed_cases"])

        # Enforce monotonic non-decreasing cumulative counts.
        # Guards against any source-to-source dip that would break the
        # log-residual in the fit.
        for c in cum_cols:
            g[c] = g[c].cummax()

        g = g.reset_index().rename(columns={"index": "date"})
        g = g[keep]
        out.append(g)

    if not out:
        raise ValueError(
            f"No rows left after cleaning. "
            f"Check the country filter (got {country!r}) and the input file."
        )

    cleaned = pd.concat(out, ignore_index=True)
    cleaned.to_csv(output_path, index=False)

    # ---- Report
    print(f"Wrote {output_path}")
    for ctry, g in cleaned.groupby("country"):
        print(f"  {ctry:8s}  {len(g):4d} daily rows  "
              f"({g['date'].min().date()} -> {g['date'].max().date()})")
    return cleaned


if __name__ == "__main__":
    # Defaults assume you run this from the repo folder. Override by
    # editing these three lines, or by passing arguments on the command line
    # if you run the script from a terminal.
    input_path = sys.argv[1] if len(sys.argv) > 1 else "ebola_drc_uganda.csv"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "ebola_drc_uganda_clean.csv"
    country = sys.argv[3] if len(sys.argv) > 3 else None

    clean(input_path, output_path, country=country)