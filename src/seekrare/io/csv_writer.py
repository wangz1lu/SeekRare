"""CSV output utilities."""

from __future__ import annotations

import pandas as pd
from pathlib import Path


def write_results(df: pd.DataFrame, path: str | Path) -> None:
    """Write ranked results to CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
