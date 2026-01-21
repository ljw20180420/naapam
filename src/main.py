import os
import pathlib
import subprocess

import matplotlib.pyplot as plt
import pandas as pd

from . import utils


def read_alg(alg_file: os.PathLike):
    with subprocess.Popen(
        args=["sed", "-e", r"N;N;s/\n/\t/g", os.fspath(alg_file)],
        stdout=subprocess.PIPE,
    ) as process:
        df_alg = pd.read_csv(
            process.stdout,
            sep="\t",
            names=[
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
                "ref",
                "query",
            ],
            keep_default_na=False,
        )

    return df_alg


def correct_index(
    root_dir: os.PathLike,
):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "align_correct", exist_ok=True)
    for alg_file in os.listdir(root_dir / "align"):
        df_alg = read_alg(root_dir / "align" / alg_file)
        df_query = pd.read_csv(
            root_dir / "align" / f"{pathlib.Path(alg_file).stem}.query",
            sep="\t",
            names=["query", "count", "ref_id"],
            keep_default_na=False,
        )
        df_alg["index"] = df_query.groupby("query").transform("ngroup")
        with open(root_dir / "align_correct" / alg_file, "w") as fd:
            for _, row in df_alg.iterrows():
                fd.write("\t".join(row.drop(["ref", "query"]).astype(str)) + "\n")
                fd.write(row["ref"] + "\n")
                fd.write(row["query"] + "\n")


def stat_read(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path("figures/hists")
    df_algs = []
    for alg_file in os.listdir(root_dir / "align_correct"):
        df_algs.append(
            read_alg(root_dir / "align_correct" / alg_file)
            .groupby("score")["count"]
            .sum()
            .reset_index()
        )

    df_alg = pd.concat(df_algs)
    df_alg["score"].plot.hist(bins=300, weights=df_alg["count"]).get_figure().savefig(
        save_dir / "score.pdf"
    )


def collect_data(
    root_dir: os.PathLike,
    min_score: int,
) -> pd.DataFrame:
    """
    Only collect data. Do not apply any annotation. Only apply read-wise filter such as score.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    df_algs = []
    for alg_file in root_dir / "align_correct":
        df_alg = read_alg(root_dir / "align_correct" / alg_file)
        if df_alg.shape[0] == 0:
            continue
        df_alg = (
            df_alg.query("score >= @min_score")
            .assign(
                size=lambda df: df.groupby("index").transform("size"),
                count=lambda df: df["count"] / df["size"],
            )
            .groupby(
                [
                    "ref_id",
                    "cut1",
                    "cut2",
                    "ref_end1",
                    "random_insertion",
                    "ref_start2",
                ]
            )["count"]
            .sum()
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

    return pd.concat(df_algs).reset_index(drop=True)


def stat_mutant(df_alg: pd.DataFrame):
    df_alg = (
        utils.up_del_size(df_alg)
        .plot.hist(bins=30, weights=df_alg["count"])
        .get_figure()
        .savefig("figures/hists/up_del_size.pdf")
    )
    plt.close("all")
    df_alg = (
        utils.down_del_size(df_alg)
        .plot.hist(bins=30, weights=df_alg["count"])
        .get_figure()
        .savefig("figures/hists/down_del_size.pdf")
    )
    plt.close("all")
    df_alg = (
        utils.rand_ins_size(df_alg)
        .plot.hist(bins=30, weights=df_alg["count"])
        .get_figure()
        .savefig("figures/hists/rand_ins_size.pdf")
    )
    plt.close("all")


def filter_mutant(
    df_alg: pd.DataFrame,
    max_up_del_size: int,
    max_down_del_size: int,
    max_rand_ins_size: int,
) -> pd.DataFrame:
    mask = (
        (utils.up_del_size(df_alg) <= max_up_del_size)
        & (utils.down_del_size(df_alg) <= max_down_del_size)
        & (utils.rand_ins_size(df_alg) <= max_rand_ins_size)
    )
    return df_alg.loc[mask].reset_index(drop=True)
