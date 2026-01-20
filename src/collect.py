import os
import pathlib
import subprocess

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
            ],
            keep_default_na=False,
        )

    return df_alg


def score(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path("figures/hists")
    df_algs = []
    for alg_file in os.listdir(root_dir / "align"):
        df_algs.append(
            read_alg(root_dir / "align" / alg_file)
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
    for alg_file in root_dir / "align":
        df_alg = read_alg(root_dir / "align" / alg_file)
        if df_alg.shape[0] == 0:
            continue
        df_alg = (
            df_alg.query("score >= @min_score")
            .groupby(
                [
                    "ref_id",
                    "ref_end1",
                    "random_insertion",
                    "ref_start2",
                    "cut1",
                    "cut2",
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
        )

        df_alg = df_alg.assign(stem=pathlib.Path(alg_file).stem)

        df_algs.append(df_alg)

    return pd.concat(df_algs).reset_index(drop=True)
