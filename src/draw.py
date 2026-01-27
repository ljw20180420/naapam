import os
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
from plotnine import aes, geom_line, ggplot, scale_color_manual

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


def over_up_del_size_on_tem(
    df_treat: pd.DataFrame,
    tem: int,
    targets: list[str],
    aggfunc: str,
) -> None:
    df_treat = (
        df_treat.query(
            """
                @cut2 - ref_start2 == @tem and \
                ref_end1 <= @cut1
            """
        )
        .assign(up_del_size=lambda df: cut1 - df["ref_end1"])
        .groupby(["cas", "up_del_size"])
        .agg(
            **{
                target: pd.NamedAgg(column=target, aggfunc=aggfunc)
                for target in targets
            }
        )
        .reset_index()
        .melt(
            id_vars=["cas", "up_del_size"],
            value_vars=targets,
            var_name="target",
            value_name="value",
        )
        .assign(group=lambda df: df["cas"] + "_" + df["target"])
    )
    if len(df_treat) == 0:
        return

    save_dir = pathlib.Path("figures/over_up_del_size_on_tem")
    os.makedirs(save_dir, exist_ok=True)
    filename = "_".join(targets) + f"_{aggfunc}_{tem}.pdf"
    (
        ggplot(df_treat, aes(x="up_del_size", y="value", color="group"))
        + geom_line()
        + scale_color_manual(values=grid_colors)
    ).save(save_dir / filename)


def draw_deletion_size(
    df: pd.DataFrame,
    targets: list[str],
    aggfunc: str,
    min_edit_total: int,
    min_edit_freq: float,
    cut1: int,
    cut2: int,
):
    df = df.query(
        """
            ref_end1 <= @cut1 and \
            count_total_correct - count_wt >= @min_edit_total and \
            (count_total_correct - count_wt) / count_total_correct - (count_total_control - count_wt_control) / count_total_control >= @min_edit_freq
        """
    )
    df = (
        df.assign(deletion_size=lambda df: cut1 - df["ref_end1"])
        .groupby(["cas", "deletion_size"])
        .agg(
            **{
                target: pd.NamedAgg(column=target, aggfunc=aggfunc)
                for target in targets
            }
        )
        .reset_index()
        .melt(
            id_vars=["cas", "deletion_size"],
            value_vars=targets,
            var_name="target",
            value_name="value",
        )
        .assign(group=lambda df: df["cas"] + "_" + df["target"])
    )
    if len(df) == 0:
        return

    save_dir = pathlib.Path("figures/draw_deletion_size")
    os.makedirs(save_dir, exist_ok=True)
    filename = "_".join(targets) + f"_{aggfunc}.pdf"
    (
        ggplot(df, aes(x="deletion_size", y="value", color="group"))
        + geom_line()
        + scale_color_manual(
            values=[
                "#FF0000",
                "#00FF00",
                "#0000FF",
                "#FFFF00",
                "#FF00FF",
                "#00FFFF",
                "#000000",
                "#888888",
                "#880000",
            ]
        )
    ).save(save_dir / filename)
