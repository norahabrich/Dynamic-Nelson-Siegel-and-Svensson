"""Data loading.

`data/bonds.csv` holds the monthly US Treasury zero-coupon yields used by
Diebold & Li (2006): end-of-month unsmoothed Fama-Bliss yields, January 1972 to
December 2000, 17 maturities from 3 to 120 months, in percent.
"""

from pathlib import Path

import numpy as np
import pandas as pd

MATURITIES = np.array([3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120], dtype=float)

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "data" / "bonds.csv"


def load_yields(path=DEFAULT_PATH) -> pd.DataFrame:
    """Return a DataFrame of yields (rows: month-end dates, columns: maturities in months)."""
    df = pd.read_csv(path, sep=";")
    dates = pd.to_datetime(df.pop("date"), format="%Y:%m") + pd.offsets.MonthEnd(0)
    df.index = pd.DatetimeIndex(dates, name="date")
    df.columns = [int(c.lstrip("M")) for c in df.columns]
    if not np.array_equal(df.columns.to_numpy(dtype=float), MATURITIES):
        raise ValueError(f"Unexpected maturities in {path}: {list(df.columns)}")
    return df.astype(float)
