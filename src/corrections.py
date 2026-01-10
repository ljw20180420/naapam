import numpy as np
import pandas as pd


def no_correction(df: pd.DataFrame, cut1: int, cut2: int) -> pd.DataFrame:
    mask_wt = (
        (df["ref_end1"] == cut1)
        & (df["ref_start2"] == cut2)
        & (df["random_insertion"] == "")
    )
    df = df.assign(
        frequency_edit=lambda df: df["count"] / (df["count_total"] - df["count_wt"]),
        frequency_total=lambda df: df["count"] / df["count_total"],
    )
    df.loc[mask_wt, "frequency_edit"] = np.nan

    return df


def kim_correction(df: pd.DataFrame, cut1: int, cut2: int) -> pd.DataFrame:
    mask_wt = (
        (df["ref_end1"] == cut1)
        & (df["ref_start2"] == cut2)
        & (df["random_insertion"] == "")
    )
    df = df.assign(
        frequency_total_kim=lambda df: (
            df["count"]
            - df["count_total"]
            * (df["count_control"] / (df["count_total_control"] + 1e-6))
        )
        / (
            df["count_total"]
            - df["count_total"]
            * (df["count_control"] / (df["count_total_control"] + 1e-6))
        ),
        frequency_edit_kim=lambda df: df["frequency_total_kim"]
        / (1 - df["count_wt"] / (df["count_total"] + 1e-6)),
    )
    df["frequency_total_kim"] = np.maximum(0, df["frequency_total_kim"])
    df.loc[mask_wt, "frequency_total_kim"] = np.nan
    df["frequency_edit_kim"] = np.maximum(0, df["frequency_edit_kim"])
    df.loc[mask_wt, "frequency_edit_kim"] = np.nan

    return df


def syn_err_correction(df: pd.DataFrame, cut1: int, cut2: int) -> pd.DataFrame:
    mask_wt = (
        (df["ref_end1"] == cut1)
        & (df["ref_start2"] == cut2)
        & (df["random_insertion"] == "")
    )
    mask = (df["count_total_control"] != 0) & ~mask_wt
    df["count_correct"] = df["count"].astype(float)
    df.loc[mask, "count_correct"] = np.maximum(
        0,
        df.loc[mask, "count"]
        - df.loc[mask, "count_total"]
        * (df.loc[mask, "count_control"] / df.loc[mask, "count_total_control"]),
    )
    df = df.assign(
        count_total_correct=lambda df: df.groupby(
            ["cas", "polq", "rep", "chip", "sample", "barcode"]
        )["count_correct"].transform("sum"),
        frequency_edit_correct=lambda df: df["count_correct"]
        / (df["count_total_correct"] - df["count_wt"]),
        frequency_total_correct=lambda df: df["count_correct"]
        / df["count_total_correct"],
    )
    df.loc[mask_wt, "frequency_edit_correct"] = np.nan
    return df
