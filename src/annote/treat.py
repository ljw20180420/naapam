import os

import matplotlib.pyplot as plt
import pandas as pd

from . import common


def count_tot(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat = df_treat.assign(
        count_tot=df_treat.groupby(["cas", "polq", "chip", "sample", "barcode"])[
            "count"
        ].transform("sum")
    )
    os.makedirs("figures/hists/treat/barcode", exist_ok=True)
    df_treat[
        ["cas", "polq", "chip", "sample", "barcode", "count_tot"]
    ].drop_duplicates().query("count_tot <= 500")["count_tot"].plot.hist(
        bins=100
    ).get_figure().savefig(
        "figures/hists/treat/barcode/count_tot.pdf"
    )
    plt.close("all")

    return df_treat


def count_wt(df_treat: pd.DataFrame) -> pd.DataFrame:
    return df_treat.merge(
        right=df_treat.query("is_wt")[
            ["cas", "polq", "chip", "sample", "barcode", "count"]
        ].rename(columns={"count": "count_wt"}),
        how="left",
        on=["cas", "polq", "chip", "sample", "barcode"],
        validate="many_to_one",
    ).assign(count_wt=lambda df: df["count_wt"].fillna(0).astype(int))


def mutant_legal(
    df_treat: pd.DataFrame,
    max_up_del_size: int,
    max_down_del_size: int,
    max_rand_ins_size: int,
    cut1: int = 50,
    cut2: int = 60,
) -> pd.DataFrame:
    df_treat["mutant_legal"] = (
        (common.up_del_size(df_treat, "treat", cut1) <= max_up_del_size)
        & (common.down_del_size(df_treat, "treat", cut2) <= max_down_del_size)
        & (common.rand_ins_size(df_treat, "treat") <= max_rand_ins_size)
    )

    return df_treat


def barcode_legal(
    df_treat: pd.DataFrame,
    min_count_tot: float,
) -> pd.DataFrame:
    df_treat["barcode_legal"] = df_treat["count_tot"] >= min_count_tot

    return df_treat
