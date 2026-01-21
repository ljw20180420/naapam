import numpy as np
import pandas as pd


def frequency_tot_kim(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat = df_treat.assign(
        frequency_tot_kim=lambda df: np.maximum(
            0,
            (
                df["count"]
                - df["count_tot"] * (df["count_ctl"] / (df["count_tot_ctl"] + 1e-6))
            )
            / (
                df["count_tot"]
                - df["count_tot"] * (df["count_ctl"] / (df["count_tot_ctl"] + 1e-6))
            ),
        )
    )
    df_treat.loc[df_treat["is_wt"], "frequency_tot_kim"] = float("nan")

    return df_treat


def frequency_nowt_kim(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat = df_treat.assign(
        frequency_nowt_kim=lambda df: np.maximum(
            0, df["frequency_tot_kim"] / (1 - df["count_wt"] / (df["count_tot"] + 1e-6))
        )
    )
    df_treat.loc[df_treat["is_wt"], "frequency_nowt_kim"] = float("nan")

    return df_treat


def frequency_tot_adj(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat = df_treat.assign(
        count_adj=lambda df: np.maximum(
            0, (df["count"] - df["count_tot"] * (df["count_ctl"] / df["count_tot_ctl"]))
        ),
    )
    df_treat.loc[df_treat["is_wt"], "count_adj"] = df_treat.loc[
        df_treat["is_wt"], "count_wt"
    ]
    df_treat = df_treat.assign(
        count_tot_adj=lambda df: df.groupby(
            ["cas", "polq", "chip", "sample", "rep", "barcode"]
        )["count_adj"].transform("sum"),
        frequency_tot_adj=lambda df: df["count_adj"] / (df["count_tot_adj"] + 1e-6),
    )

    return df_treat


def frequency_nowt_adj(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat = df_treat.assign(
        frequency_nowt_adj=lambda df: df["count_adj"]
        / (df["count_tot_adj"] - df["count_wt"] + 1e-6),
    )
    df_treat.loc[df_treat["is_wt"], "frequency_nowt_adj"] = float("nan")
    return df_treat
