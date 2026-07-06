"""Statistical factor models on S&P 500 daily returns — a Python walkthrough.

Mirrors the R version (`../hands-on_R_factors/sp500/pca_sp500.R`) using the
standard Python data-analysis stack: pandas for wrangling, scikit-learn for PCA
and factor analysis, matplotlib for the scree plot.

Pipeline:
    prices (long) -> daily simple returns -> wide (date x stock) matrix
    -> PCA (how many common factors?) -> 3-factor varimax model (what are they?)

Run:  python python/pca_sp500.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.preprocessing import StandardScaler

DATA = Path(__file__).resolve().parent.parent / "hands-on_R_factors" / "sp500" / "sp500_prices.csv"
FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def load_returns() -> tuple[pd.DataFrame, dict[str, str]]:
    """Long prices -> wide matrix of daily simple returns (date x stock).

    Also returns a {symbol: sector} map so we can *check* what each estimated
    factor represents instead of eyeballing tickers.
    """
    df = pd.read_csv(DATA, usecols=["symbol", "date", "adjusted", "sector"])
    df["date"] = pd.to_datetime(df["date"])
    sector = df.drop_duplicates("symbol").set_index("symbol")["sector"].to_dict()
    df = df.sort_values(["symbol", "date"])
    # simple return within each stock, then pivot to a wide panel
    df["ret"] = df.groupby("symbol")["adjusted"].pct_change()
    wide = df.pivot(index="date", columns="symbol", values="ret").dropna()
    return wide, sector


def how_many_factors(returns: pd.DataFrame) -> PCA:
    """PCA on standardized returns; report variance explained by the top PCs."""
    X = StandardScaler().fit_transform(returns.values)
    pca = PCA().fit(X)
    evr = pca.explained_variance_ratio_
    print(f"Panel: {returns.shape[0]:,} days x {returns.shape[1]} stocks")
    print(f"PC1 explains {evr[0]:.1%}; top 3 explain {evr[:3].sum():.1%}; "
          f"top 5 explain {evr[:5].sum():.1%}")

    # Scree plot
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(np.arange(1, 16), evr[:15] * 100, "o-")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("% of variance explained")
    ax.set_title("Scree plot — a few common factors drive most co-movement")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "scree.png", dpi=150)
    plt.close(fig)
    return pca


def name_the_factors(returns: pd.DataFrame, sector: dict[str, str],
                     k: int = 3) -> None:
    """3-factor varimax model; show the top-loading stocks + their sectors."""
    X = StandardScaler().fit_transform(returns.values)
    fa = FactorAnalysis(n_components=k, rotation="varimax", random_state=0).fit(X)
    loadings = pd.DataFrame(fa.components_.T, index=returns.columns,
                            columns=[f"Factor{i+1}" for i in range(k)])
    for f in loadings.columns:
        col = loadings[f]
        # varimax signs are arbitrary: flip so the dominant loadings are positive
        if col.loc[col.abs().idxmax()] < 0:
            col = -col
        top = col.reindex(col.abs().sort_values(ascending=False).index).head(8)
        tbl = pd.DataFrame({"loading": top.round(2),
                            "sector": [sector.get(s, "?") for s in top.index]})
        # the modal sector among the top loaders = what this factor "is"
        theme = tbl["sector"].mode().iat[0]
        print(f"\n{f}  (looks like: {theme})")
        print(tbl.to_string())


def main() -> None:
    returns, sector = load_returns()
    how_many_factors(returns)
    name_the_factors(returns, sector, k=3)
    print(f"\nScree plot saved to {FIG_DIR / 'scree.png'}")


if __name__ == "__main__":
    main()
