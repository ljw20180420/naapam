import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import special

from . import utils


def correct_alg(root_dir: os.PathLike, temperature: float):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "align" / "correct", exist_ok=True)
    for alg_file in os.listdir(root_dir / "align" / "raw"):
        df_alg = utils.read_alg(root_dir / "align" / "raw" / alg_file)
        df_query = pd.read_csv(
            root_dir / "query" / "found" / f"{pathlib.Path(alg_file).stem}.query",
            sep="\t",
            names=["query", "count", "ref_id"],
            keep_default_na=False,
        )
        df_alg["index"] = df_query.groupby("query").transform("ngroup")

        chip = utils.infer_chip(alg_file)
        df_ref = pd.read_feather(root_dir / "control" / "hq_mut" / f"{chip}.feather")[
            ["count"]
        ].reset_index(names="ref_id")
        df_alg = df_alg.merge(
            right=df_ref[["ref_id", "count"]].rename(columns={"count": "count_ref"}),
            how="left",
            on=["ref_id"],
            validate="many_to_one",
        ).assign(
            count_distri=lambda df: df["count"]
            * df.groupby("index")["score"].transform(
                lambda score, temperature=temperature: pd.Series(
                    special.softmax(score / temperature), name=score.name
                )
            )
        )

        with open(root_dir / "align" / "correct" / alg_file, "w") as fd:
            for _, row in df_alg.iterrows():
                fd.write(
                    "\t".join(
                        row[
                            [
                                "index",
                                "count",
                                "score",
                                "ref_id",
                                "updangle",
                                "ref_start1",
                                "query_start1",
                                "ref_end1",
                                "query_end1",
                                "random_insertion",
                                "ref_start2",
                                "query_start2",
                                "ref_end2",
                                "query_end2",
                                "downdangle",
                                "cut1",
                                "cut2",
                                "count_ref",
                                "count_distri",
                            ]
                        ].astype(str)
                    )
                    + "\n"
                )
                fd.write(row["ref"] + "\n")
                fd.write(row["query"] + "\n")


def stat_read(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path("figures/main/stat_read")
    os.makedirs(save_dir, exist_ok=True)
    df_algs = []
    for alg_file in os.listdir(root_dir / "align" / "correct"):
        df_algs.append(
            utils.read_alg(root_dir / "align" / "correct" / alg_file)
            .groupby("score")["count_distri"]
            .sum()
            .rename("count")
            .reset_index()
        )

    df_alg = pd.concat(df_algs)
    df_alg["score"].plot.hist(bins=300, weights=df_alg["count"]).get_figure().savefig(
        save_dir / "score.pdf"
    )


def collect_data(
    root_dir: os.PathLike,
    min_score: int,
):
    """
    Only collect data. Do not apply any annotation. Only apply read-wise filter such as score.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "main" / "treat" / "full", exist_ok=True)
    os.makedirs(root_dir / "main" / "control" / "full", exist_ok=True)
    df_algs = []
    for alg_file in root_dir / "align" / "correct":
        df_alg = utils.read_alg(root_dir / "align" / "correct" / alg_file)
        if df_alg.shape[0] == 0:
            continue

        df_alg = (
            df_alg.query("score >= @min_score")
            .groupby(
                [
                    "ref_id",
                    "cut1",
                    "cut2",
                    "ref_end1",
                    "random_insertion",
                    "ref_start2",
                ]
            )["count_distri"]
            .sum()
            .rename("count")
            .reset_index()
            .astype(
                {
                    "ref_end1": "int8",
                    "ref_start2": "int8",
                    "cut1": "int8",
                    "cut2": "int8",
                }
            )
            .assign(stem=pathlib.Path(alg_file).stem)
        )
        df_algs.append(df_alg)

    df_alg = pd.concat(df_algs).assign(cas=lambda df: df["stem"].map(utils.infer_cas))
    df_alg.query("cas != 'control'").drop(columns="cas").reset_index(
        drop=True
    ).to_feather(root_dir / "main" / "treat" / "full" / "treat.feather")
    df_alg.query("cas == 'control'").drop(columns="cas").reset_index(
        drop=True
    ).to_feather(root_dir / "main" / "control" / "full" / "treat.feather")


def stat_mutant(root_dir: os.PathLike):
    save_dir = pathlib.Path("figures/main/stat_mutant")
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(root_dir / "main" / "treat" / "full" / "treat.feather")

    utils.up_del_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "up_del_size.pdf")
    plt.close("all")

    utils.down_del_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "down_del_size.pdf")
    plt.close("all")

    utils.rand_ins_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "rand_ins_size.pdf")
    plt.close("all")

    utils.freq_mutant(df_treat).plot.hist(bins=100).get_figure().savefig(
        save_dir / "freq_mutant.pdf"
    )
    plt.close("all")


def filter_mutant(
    root_dir: os.PathLike,
    max_up_del_size: int,
    max_down_del_size: int,
    max_rand_ins_size: int,
    max_freq_mutant: float,
):
    """
    Do not filter mutant because missing mutant are treated as count 0. Set mutant count to nan to exclude it from all statistics involving count.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "main" / "treat" / "filter" / "mutant", exist_ok=True)

    df_treat = pd.read_feather(root_dir / "main" / "treat" / "full" / "treat.feather")

    mask = (
        (utils.up_del_size(df_treat) <= max_up_del_size)
        & (utils.down_del_size(df_treat) <= max_down_del_size)
        & (utils.rand_ins_size(df_treat) <= max_rand_ins_size)
        & ((utils.freq_mutant(df_treat) <= max_freq_mutant) | utils.is_wt(df_treat))
    )
    df_treat.loc[~mask, "count"] = float("nan")

    df_treat.to_feather(
        root_dir / "main" / "treat" / "filter" / "mutant" / "treat.feather"
    )


def stat_ref(root_dir: os.PathLike, min_count_tot: int):
    save_dir = pathlib.Path("figures/hists/ref")
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "main" / "treat" / "filter" / "mutant" / "treat.feather"
    )

    df_treat.groupby(["stem", "ref_id"])["count"].sum().plot.hist(
        bins=100
    ).get_figure().savefig(save_dir / "count_tot.pdf")
    plt.close("all")

    df_treat.assign(freq_nowt=utils.freq_nowt(df_treat))[
        ["stem", "ref_id", "freq_nowt"]
    ].drop_duplicates()["freq_nowt"].plot.hist(bins=100).get_figure().savefig(
        save_dir / "freq_nowt.pdf"
    )
    plt.close("all")

    for tem in range(1, 5):
        df_treat_freq = (
            df_treat.assign(
                **{
                    "count_tot": lambda df: df.groupby(["stem", "ref_id"])[
                        "count"
                    ].transform("sum"),
                    f"freq_tem{tem}": lambda df: utils.freq_temN(df, tem),
                    f"freq_tem{tem}_blunt": lambda df: utils.freq_temN_blunt(df, tem),
                    f"freq_tem{tem}_dummy_rel_blunt": lambda df: utils.freq_temN_dummy_rel_blunt(
                        df, tem
                    ),
                }
            )
            .query("count_tot >= @min_count_tot")[
                [
                    "stem",
                    "ref_id",
                    f"freq_tem{tem}",
                    f"freq_tem{tem}_blunt",
                    f"freq_tem{tem}_dummy_rel_blunt",
                ]
            ]
            .drop_duplicates()
        )

        for column in [
            f"freq_tem{tem}",
            f"freq_tem{tem}_blunt",
            f"freq_tem{tem}_dummy_rel_blunt",
        ]:
            df_treat_freq[column].plot.hist(bins=100).get_figure().savefig(
                save_dir / f"{column}.pdf"
            )
            plt.close("all")


def filter_ref(
    root_dir: os.PathLike,
    min_count_tot: int,
    max_freq_nowt: float,
    max_freq_temN: dict[int, float],
    max_freq_temN_blunt: dict[int, float],
    max_freq_temN_dummy_rel_blunt: dict[int, float],
):
    """
    Use positive mask because nan compare always return False.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "main" / "treat" / "filter" / "ref", exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "main" / "treat" / "filter" / "mutant" / "treat.feather"
    )

    mask = (
        df_treat.groupby(["stem", "ref_id"])["count"].transform("sum") >= min_count_tot
    )
    mask = mask & (utils.freq_nowt(df_treat) <= max_freq_nowt)
    for tem in range(1, 5):
        mask = mask & (utils.freq_temN(df_treat, tem) <= max_freq_temN[tem])
        mask = mask & (utils.freq_temN_blunt(df_treat, tem) <= max_freq_temN_blunt[tem])
        mask = mask & (
            utils.freq_temN_dummy_rel_blunt(df_treat, tem)
            <= max_freq_temN_dummy_rel_blunt[tem]
        )

    df_treat.loc[mask].reset_index(drop=True).to_feather(
        root_dir / "main" / "treat" / "filter" / "ref" / "treat.feather"
    )
