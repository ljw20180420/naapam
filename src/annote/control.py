import os

import matplotlib.pyplot as plt
import pandas as pd


def b_synerr(df_control: pd.DataFrame) -> pd.Series:
    b_synerr = 1 - df_control["count_wt"] / df_control["count_tot"]
    os.makedirs("figures/hists/control/barcode", exist_ok=True)
    df_control.assign(b_synerr=b_synerr)[
        ["chip", "rep", "barcode", "b_synerr"]
    ].drop_duplicates()["b_synerr"].plot.hist(bins=100).get_figure().savefig(
        "figures/hists/control/barcode/b_synerr.pdf"
    )
    plt.close("all")
    return b_synerr


def b_temNerr(df_control: pd.DataFrame, tem: int) -> pd.Series:
    b_temNerr = df_control[f"count_tem{tem}"] / df_control["count_tot"]
    os.makedirs("figures/hists/control/barcode", exist_ok=True)
    df_control.assign(b_temNerr=b_temNerr)[
        ["chip", "rep", "barcode", "b_temNerr"]
    ].drop_duplicates()["b_temNerr"].plot.hist(bins=100).get_figure().savefig(
        f"figures/hists/control/barcode/b_tem{tem}err.pdf"
    )
    plt.close("all")
    return b_temNerr


def barcode_legal(
    df_control: pd.DataFrame,
    max_b_synerr: float,
    max_b_temNerr: dict[int, float],
) -> pd.DataFrame:
    """
    A proper barcode legal condition satisfies that if a barcode is missing in control, then the barcode is legal. A non-proper example is the lower bound of total read counts of a barcode. A missing barcode has 0 count, which is always illegal.
    """
    df_control["barcode_legal"] = b_synerr(df_control) <= max_b_synerr
    for tem in range(1, 5):
        df_control["barcode_legal"] = df_control["barcode_legal"] & (
            b_temNerr(df_control, tem) <= max_b_temNerr[tem]
        )
    return df_control
