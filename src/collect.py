import os
import pathlib
import re

import pandas as pd


def score(tables: list[os.PathLike]):
    pass


def infer_file_property(table: os.PathLike) -> tuple[str, str, int, str]:
    file_stem = pathlib.Path(os.fspath(table)).stem
    chip = re.search(r"-(.+)\.fq$", file_stem).group(1)
    sample = re.search(r"^(.*?)-", file_stem).group(1)
    if (
        file_stem.startswith("A2")
        or file_stem.startswith("A7")
        or file_stem.startswith("D2")
    ):
        cas = "spycas9"
        polq = "wt"
    elif (
        file_stem.startswith("36t")
        or file_stem.startswith("B2")
        or file_stem.startswith("x")
    ):
        cas = "spymac"
        polq = "wt"
    elif file_stem.startswith("i10t") or file_stem.startswith("i83"):
        cas = "ispymac"
        polq = "wt"
    elif file_stem.startswith("h73") or file_stem.startswith("h80"):
        cas = "spycas9"
        polq = "hel"
    elif file_stem.startswith("p79") or file_stem.startswith("p125"):
        cas = "spycas9"
        polq = "pol"
    elif (
        file_stem.startswith("Q101")
        or file_stem.startswith("Q102")
        or file_stem.startswith("Q192")
    ):
        cas = "spycas9"
        polq = "ko"
    elif file_stem.startswith("wt1-"):
        cas = "control"
        polq = "wt"
    elif file_stem.startswith("wt2-"):
        cas = "control"
        polq = "wt"
    elif file_stem.startswith("wt11"):
        cas = "control"
        polq = "wt"
    elif file_stem.startswith("wt21"):
        cas = "control"
        polq = "wt"
    else:
        raise Exception("unknown stem")

    return cas, polq, chip, sample


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

        cas, polq, chip, sample = infer_file_property(table)
        df = df.assign(
            cas=cas,
            polq=polq,
            chip=chip,
            sample=sample,
        )

        dfs.append(df)

    df = pd.concat(dfs).reset_index(drop=True)

    return df


def read_feather(file: os.PathLike):
    return pd.read_feather(file).assign(
        random_insertion=lambda df: df["random_insertion"].fillna("")
    )
