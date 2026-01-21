import os

import matplotlib.pyplot as plt
import pandas as pd


def m_synerr(df_control: pd.DataFrame) -> pd.Series:
    m_synerr = df_control["count"] / df_control["count_tot"]
    m_synerr.loc[df_control["is_wt"]] = float("nan")
    os.makedirs("figures/hists/control/mutant", exist_ok=True)
    m_synerr.plot.hist(bins=100).get_figure().savefig(
        "figures/hists/control/mutant/m_synerr.pdf"
    )
    plt.close("all")
    return m_synerr


def mutant_legal(df_control: pd.DataFrame, max_m_synerr: float) -> pd.DataFrame:
    """
    A proper mutant legal condition satisfies that if a mutant is missing in control, then the mutant is legal. A non-proper example is the lower bound of total read counts of a mutant. A missing mutant has 0 count, which is always illegal.

    If a mutant legal condition is only determined by the mutant itself but not its count (like mutant indel size), then it is redundant to apply the condition to both control and treat.
    """
    df_control["mutant_legal"] = (m_synerr(df_control) <= max_m_synerr) | df_control[
        "is_wt"
    ]

    return df_control


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
