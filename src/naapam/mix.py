import os
import pathlib

import pandas as pd

from . import utils


def duplicate_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "dup", exist_ok=True)
    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "filter" / "ref" / "treat.feather"
    )
    df_treat = df_treat.assign(
        chip=lambda df: df["stem"].map(utils.infer_chip),
        time=lambda df: df["stem"].map(utils.infer_time),
        wt="wt1 wt2",
    )
    df_treat.loc[(df_treat["chip"] == "g2n") & (df_treat["time"] == 2), "wt"] = (
        "wt1 wt2 wt11 wt21"
    )
    (
        df_treat.assign(wt=lambda df: df["wt"].str.split())
        .explode("wt")
        .reset_index(drop=True)
        .assign(stem=lambda df: df["stem"] + "_" + df["wt"])
        .to_feather(root_dir / "anaylze" / "treat" / "dup" / "treat.feather")
    )


def duplicate_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "control" / "dup", exist_ok=True)
    df_control = pd.read_feather(
        root_dir / "analyze" / "control" / "full" / "control.feather"
    )

    df_control = df_control.assign(
        chip=lambda df: df["stem"].map(utils.infer_chip),
        time=lambda df: df["stem"].map(utils.infer_time).astype(str),
        wt=lambda df: df["stem"].map(utils.infer_wt),
    )

    df_control.loc[
        (df_control["wt"] == "wt1")
        & (df_control["chip"] == "a2")
        & (df_control["time"] == "1"),
        "time",
    ] = "1 3"
    df_control.loc[
        (df_control["wt"] == "wt1")
        & (df_control["chip"] == "a2")
        & (df_control["time"] == "2"),
        "time",
    ] = "2 4"
    df_control.loc[
        (df_control["wt"] == "wt1")
        & (df_control["chip"] == "g1n")
        & (df_control["time"] == "2"),
        "time",
    ] = "2 1"
    df_control.loc[
        (df_control["wt"] == "wt1")
        & (df_control["chip"] == "g2n")
        & (df_control["time"] == "4"),
        "time",
    ] = "4 3"
    df_control.loc[
        (df_control["wt"] == "wt1")
        & (df_control["chip"] == "g3n")
        & (df_control["time"] == "1"),
        "time",
    ] = "1 2 3 4"
    df_control.loc[
        (df_control["wt"] == "wt2")
        & (df_control["chip"] == "a2")
        & (df_control["time"] == "3"),
        "time",
    ] = "3 4"
    df_control.loc[
        (df_control["wt"] == "wt2")
        & (df_control["chip"] == "a3")
        & (df_control["time"] == "3"),
        "time",
    ] = "3 4"
    df_control.loc[
        (df_control["wt"] == "wt2")
        & (df_control["chip"] == "g3n")
        & (df_control["time"] == "1"),
        "time",
    ] = "1 2 3 4"

    (
        df_control.assign(time=lambda df: df["time"].str.split())
        .explode("time")
        .reset_index(drop=True)
        .assign(stem=lambda df: df["stem"] + "_" + df["time"])
        .astype({"time": int})
        .to_feather(root_dir / "analyze" / "control" / "dup" / "control.feather")
    )


def merge(root_dir: os.PathLike):
    """
    count_wt_ctl and count_tot_ctl must be calculated before left merge because information about control will loss after left merge.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "merge", exist_ok=True)
    df_treat = pd.read_feather(root_dir / "analyze" / "treat" / "dup" / "treat.feather")
    df_control = (
        pd.read_feather(root_dir / "analyze" / "control" / "dup" / "treat.feather")
        .assign(
            count_wt_ctl=lambda df: utils.count_wt(df),
            count_tot_ctl=lambda df: df.groupby(["stem", "ref_id"])["count"].transform(
                "sum"
            ),
        )
        .rename(columns={"count": "count_ctl"})
    )
    on = [
        "chip",
        "time",
        "wt",
        "ref_id",
        "ref_end1",
        "ref_start2",
        "random_insertion",
    ]
    df_treat = df_treat.merge(
        right=df_control[on + ["count_ctl", "count_wt_ctl", "count_tot_ctl"]],
        how="left",
        on=on,
        validate="many_to_one",
    ).assign(
        count_ctl=lambda df: df["count_ctl"].fillna(0),
        count_wt_ctl=lambda df: df["count_wt_ctl"].fillna(0),
        count_tot_ctl=lambda df: df["count_tot_ctl"].fillna(0),
    )

    df_treat.to_feather(root_dir / "analyze" / "treat" / "merge" / "treat.feather")
