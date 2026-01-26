import os
import pathlib

import pandas as pd

from . import utils


def duplicate_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "main" / "treat" / "dup", exist_ok=True)
    df_treat = pd.read_feather(
        root_dir / "main" / "treat" / "filter" / "ref" / "treat.feather"
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
        .to_feather(root_dir / "main" / "treat" / "dup" / "treat.feather")
    )


def duplicate_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "main" / "control" / "dup", exist_ok=True)
    df_control = pd.read_feather(
        root_dir / "main" / "control" / "full" / "control.feather"
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
        .astype({"time": int})
        .to_feather(root_dir / "main" / "control" / "dup" / "control.feather")
    )


def merge_mutant(df_treat: pd.DataFrame, df_control: pd.DataFrame) -> pd.DataFrame:
    return df_treat.merge(
        df_control[
            [
                "chip",
                "rep",
                "barcode",
                "ref_end1",
                "ref_start2",
                "random_insertion",
                "count",
                "mutant_legal",
            ]
        ].rename(
            columns={
                "count": "count_ctl",
                "mutant_legal": "mutant_legal_ctl",
            }
        ),
        how="left",
        on=[
            "chip",
            "rep",
            "barcode",
            "ref_end1",
            "ref_start2",
            "random_insertion",
        ],
        validate="many_to_one",
    ).assign(
        count_ctl=lambda df: df["count_ctl"].fillna(0).astype(int),
        mutant_legal_ctl=lambda df: df["mutant_legal_ctl"].fillna(1).astype(bool),
    )


def merge_barcode(df_treat: pd.DataFrame, df_control: pd.DataFrame) -> pd.DataFrame:
    return df_treat.merge(
        df_control[
            [
                "chip",
                "rep",
                "barcode",
                "count_tot",
                "barcode_legal",
            ]
        ]
        .drop_duplicates()
        .rename(
            columns={
                "count_tot": "count_tot_ctl",
                "barcode_legal": "barcode_legal_ctl",
            }
        ),
        how="left",
        on=[
            "chip",
            "rep",
            "barcode",
        ],
        validate="many_to_one",
    ).assign(
        count_tot_ctl=lambda df: df["count_tot_ctl"].fillna(0).astype(int),
        barcode_legal_ctl=lambda df: df["barcode_legal_ctl"].fillna(1).astype(bool),
    )
