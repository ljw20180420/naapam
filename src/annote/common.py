import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def up_del_size(df: pd.DataFrame) -> pd.Series:
    up_del_size: pd.Series = np.maximum(df["cut1"] - df["ref_end1"], 0)
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    up_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "up_del_size.pdf"
    )
    plt.close("all")
    return up_del_size


def down_del_size(df: pd.DataFrame) -> pd.Series:
    down_del_size: pd.Series = np.maximum(df["ref_start2"] - df["cut2"], 0)
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    down_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "down_del_size.pdf"
    )
    plt.close("all")
    return down_del_size


def rand_ins_size(df: pd.DataFrame, dir: str) -> pd.Series:
    rand_ins_size: pd.Series = df["random_insertion"].str.len()
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    rand_ins_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "rand_ins_size.pdf"
    )
    plt.close("all")
    return rand_ins_size


def is_wt(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        is_wt=lambda df: (
            (df["ref_end1"] == df["cut1"])
            & (df["ref_start2"] == df["cut2"])
            & (df["random_insertion"] == "")
        )
    )
