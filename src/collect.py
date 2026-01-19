import os
import pathlib
import subprocess

import pandas as pd

from . import utils


def score(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path("figures/hists")
    df_treats = []
    for file in os.listdir(root_dir / "align"):
        with subprocess.Popen(
            args=["sed", "-e", r"N;N;s/\n/\t/g", root_dir / "align" / file],
            stdout=subprocess.PIPE,
        ) as process:
            df_treats.append(
                pd.read_csv(
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
                    ],
                )
                .groupby("score")["count"]
                .sum()
                .reset_index()
            )

    df_treat: pd.DataFrame = pd.concat(df_treats)
    df_treat["score"].plot.hist(
        bins=300, weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "score.pdf")


def collect_data(
    tables: list[os.PathLike],
    cut1: int,
    cut2: int,
    min_score: int,
) -> pd.DataFrame:
    """
    Only collect data. Do not apply any annotation. Only apply read-wise filter.
    """

    # Collect data
    dfs = []
    for table in tables:
        df = pd.read_csv(table, sep="\t", header=0, keep_default_na=False)
        if len(df) == 0:
            continue
        df = (
            df.query("score >= @min_score")
            .groupby(
                [
                    "barcode",
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
                }
            )
        )

        df = df.assign(
            cas=utils.infer_cas(table),
            polq=utils.infer_polq(table),
            chip=utils.infer_chip(table) + "-" + utils.infer_time(table),
            sample=utils.infer_sample(table),
        )

        dfs.append(df)

    df = pd.concat(dfs).reset_index(drop=True)

    return df


def read_feather(file: os.PathLike):
    return pd.read_feather(file).assign(
        random_insertion=lambda df: df["random_insertion"].fillna("")
    )
