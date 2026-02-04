import os
import pathlib

import pandas as pd
from plotnine import aes, geom_line, ggplot, scale_color_manual

from . import utils

grid_colors = [
    "#FF0000",
    "#00FF00",
    "#0000FF",
    "#FFFF00",
    "#FF00FF",
    "#00FFFF",
    "#000000",
    "#888888",
    "#888800",
    "#008888",
    "#880088",
    "#8888FF",
    "#FF8888",
    "#88FF88",
    "#880000",
    "#88FF00",
    "#88FFFF",
    "#8800FF",
    "#008800",
    "#FF8800",
    "#FF88FF",
    "#0088FF",
    "#000088",
    "#FF0088",
    "#FFFF88",
    "#00FF88",
]


def mean_over_up_del_size_on_tem(root_dir: os.PathLike, target: str):
    root_dir = pathlib.Path(os.fspath(root_dir))
    df_treat = pd.read_csv(
        root_dir / "analyze" / "treat" / "correct" / "treat.csv",
        header=0,
        na_values=["NA"],
        keep_default_na=False,
    )
    for tem in range(1, 5):
        save_dir = (
            root_dir
            / "figures"
            / "analyze"
            / "mean_freq_over_up_del_size_on_tem"
            / str(tem)
        )
        os.makedirs(save_dir, exist_ok=True)
        mean_over_up_del_size_on_tem_inner(
            df_treat,
            tem,
            target,
            save_dir / "all.pdf",
        )
        for stem in df_treat["stem"].unique().tolist():
            df_treat_stem = df_treat.query("stem == @stem").reset_index(drop=True)
            mean_over_up_del_size_on_tem_inner(
                df_treat_stem,
                tem,
                target,
                save_dir / f"{stem}.pdf",
            )


def mean_over_up_del_size_on_tem_inner(
    df_treat: pd.DataFrame, tem: int, target: str, save_file: os.PathLike
):
    df_tem = (
        df_treat.query(
            """
                cut2 - ref_start2 == @tem and \
                ref_end1 <= cut1
            """
        )
        .reset_index(drop=True)
        .assign(
            up_del_size=lambda df: df["cut1"] - df["ref_end1"],
        )
    )
    if not df_tem["legal"].any():
        return

    df_value = (
        utils.pivot_value(df=df_tem, value=target, column="up_del_size")
        .assign(cas=lambda df: df["stem"].map(utils.infer_cas))
        .drop(columns=["stem", "ref_id"])
        .groupby("cas")
        .sum()
    )
    df_legal = (
        utils.pivot_legal(df=df_tem, column="up_del_size")
        .assign(cas=lambda df: df["stem"].map(utils.infer_cas))
        .drop(columns=["stem", "ref_id"])
        .groupby("cas")
        .sum()
    )
    value_vars = df_legal.columns

    df_mean = (
        (df_value / df_legal)
        .reset_index()
        .melt(
            id_vars="cas",
            value_vars=value_vars,
            var_name="up_del_size",
            value_name=target,
        )
    )

    (
        ggplot(df_mean, aes(x="up_del_size", y=target, group="cas", color="cas"))
        + geom_line()
        + scale_color_manual(values=grid_colors)
    ).save(save_file)
