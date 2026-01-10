import pandas as pd


def duplicate_treat(df_treat: pd.DataFrame) -> pd.DataFrame:
    df_treat["rep"] = "wt1 wt2"
    df_treat.loc[df_treat["chip"] == "g2n-2", "rep"] = "wt1 wt2 wt11 wt21"
    return (
        df_treat.assign(rep=lambda df: df["rep"].str.split())
        .explode("rep")
        .reset_index(drop=True)
    )


def duplicate_control(df_control: pd.DataFrame) -> pd.DataFrame:
    df_control.loc[
        (df_control["rep"] == "wt1") & (df_control["chip"] == "a2-1"), "chip"
    ] = "a2-1 a2-3"
    df_control.loc[
        (df_control["rep"] == "wt1") & (df_control["chip"] == "a2-2"), "chip"
    ] = "a2-2 a2-4"
    df_control.loc[
        (df_control["rep"] == "wt1") & (df_control["chip"] == "g1n-2"), "chip"
    ] = "g1n-2 g1n-1"
    df_control.loc[
        (df_control["rep"] == "wt1") & (df_control["chip"] == "g2n-4"), "chip"
    ] = "g2n-4 g2n-3"
    df_control.loc[
        (df_control["rep"] == "wt1") & (df_control["chip"] == "g3n-1"), "chip"
    ] = "g3n-1 g3n-2 g3n-3 g3n-4"
    df_control.loc[
        (df_control["rep"] == "wt2") & (df_control["chip"] == "a2-3"), "chip"
    ] = "a2-3 a2-4"
    df_control.loc[
        (df_control["rep"] == "wt2") & (df_control["chip"] == "a3-3"), "chip"
    ] = "a3-3 a3-4"
    df_control.loc[
        (df_control["rep"] == "wt2") & (df_control["chip"] == "g3n-1"), "chip"
    ] = "g3n-1 g3n-2 g3n-3 g3n-4"

    return (
        df_control.assign(chip=lambda df: df["chip"].str.split())
        .explode("chip")
        .reset_index(drop=True)
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
