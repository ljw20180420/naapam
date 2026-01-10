import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def up_del_size(df: pd.DataFrame, dir: str, cut1: int = 50) -> pd.Series:
    up_del_size: pd.Series = np.maximum(cut1 - df["ref_end1"], 0)
    os.makedirs(f"figures/hists/{dir}/mutant", exist_ok=True)
    up_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        f"figures/hists/{dir}/mutant/up_del_size.pdf"
    )
    plt.close("all")
    return up_del_size


def down_del_size(df: pd.DataFrame, dir: str, cut2: int = 60) -> pd.Series:
    down_del_size: pd.Series = np.maximum(df["ref_start2"] - cut2, 0)
    os.makedirs(f"figures/hists/{dir}/mutant", exist_ok=True)
    down_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        f"figures/hists/{dir}/mutant/down_del_size.pdf"
    )
    plt.close("all")
    return down_del_size


def rand_ins_size(df: pd.DataFrame, dir: str) -> pd.Series:
    rand_ins_size: pd.Series = df["random_insertion"].str.len()
    os.makedirs(f"figures/hists/{dir}/mutant", exist_ok=True)
    rand_ins_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        f"figures/hists/{dir}/mutant/rand_ins_size.pdf"
    )
    plt.close("all")
    return rand_ins_size


def is_wt(df: pd.DataFrame, cut1: int = 50, cut2: int = 60) -> pd.DataFrame:
    return df.assign(
        is_wt=lambda df: (
            (df["ref_end1"] == cut1)
            & (df["ref_start2"] == cut2)
            & (df["random_insertion"] == "")
        )
    )
