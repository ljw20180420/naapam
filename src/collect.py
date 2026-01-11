import os
import pathlib
import re

import pandas as pd

from . import utils


def score(tables: list[os.PathLike]):
    pass


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
