import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import aes, after_stat, geom_histogram, ggplot

from . import utils


def stat_indel(
    df_new: pd.DataFrame, df_old: pd.DataFrame, save_dir: os.PathLike
) -> pd.DataFrame:
    save_dir = pathlib.Path(os.fspath(save_dir))

    df_new = df_new.assign(
        del_size=lambda df: np.maximum(0, df["cut1"] - df["ref_end1"])
        + np.maximum(0, df["ref_start2"] - df["cut2"]),
        tem_ins_size=lambda df: np.maximum(0, df["cut2"] - df["ref_start2"]),
        rand_ins_size=lambda df: df["random_insertion"].str.len(),
        is_wt=lambda df: (df["ref_end1"] == df["cut1"])
        & (df["ref_start2"] == df["cut2"])
        & (df["random_insertion"] == ""),
    )
    df_old = df_old.assign(
        del_size=lambda df: np.maximum(0, df["cut1"] - df["ref_end1"])
        + np.maximum(0, df["ref_start2"] - df["cut2"]),
        tem_ins_size=lambda df: np.maximum(0, df["cut2"] - df["ref_start2"]),
        rand_ins_size=lambda df: df["random_insertion"].str.len(),
        is_wt=lambda df: (df["ref_end1"] == df["cut1"])
        & (df["ref_start2"] == df["cut2"])
        & (df["random_insertion"] == ""),
    )

    df_agg = {}
    for m_size in [
        "del_size",
        "tem_ins_size",
        "rand_ins_size",
    ]:
        df_new[m_size] = df_new[m_size].clip(upper=30)
        df_old[m_size] = df_old[m_size].clip(upper=30)
        df_agg[m_size] = pd.concat(
            [
                df_new.groupby(["is_wt", m_size])["count"]
                .sum()
                .reset_index()
                .assign(epoch="new"),
                df_old.groupby(["is_wt", m_size])["count"]
                .sum()
                .reset_index()
                .assign(epoch="old"),
            ]
        ).reset_index(drop=True)

        (
            ggplot(
                df_agg[m_size].query("not is_wt"),
                aes(x=m_size, fill="epoch", weight="count"),
            )
            + geom_histogram(
                aes(y=after_stat("density")), position="dodge", binwidth=1, bins=30
            )
        ).save(save_dir / f"{m_size}.pdf")

    return df_agg


def read_new(alg_file: os.PathLike) -> pd.DataFrame:
    return utils.read_alg(
        alg_file,
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


def read_old(alg_file: os.PathLike) -> pd.DataFrame:
    return pd.read_csv(
        alg_file,
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


def collect_indel(root_dir: os.PathLike) -> pd.DataFrame:
    root_dir = pathlib.Path(os.fspath(root_dir))

    epochs = []
    stems = []
    df_cat = {
        "del_size": [],
        "tem_ins_size": [],
        "rand_ins_size": [],
    }
    new_dir = root_dir / "align" / "raw"
    old_dir = pathlib.Path("/home/ljw/sdc1/sx_data/SX/algs")
    for new_file in os.listdir(new_dir):
        new_stem = pathlib.Path(new_file).stem.lower()
        for old_file in os.listdir(old_dir):
            old_stem = old_file.replace(".R2.fq.alg.gz", "").lower()
            if old_stem == new_stem:
                break
        if old_stem != new_stem:
            sys.stderr.write(f"cannot find old stem for {new_stem}\n")
            continue
        stem = new_stem

        save_dir = root_dir / "ts" / "collect_indel" / stem
        os.makedirs(save_dir, exist_ok=True)

        df_new = read_new(new_dir / new_file)
        df_old = read_old(old_dir / old_file)

        df_agg = stat_indel(df_new, df_old, save_dir)
        for m_size, df in df_agg.items():
            df_cat[m_size].append(df.assign(stem=stem))

    for m_size in df_cat.keys():
        df_cat[m_size] = pd.concat(df_cat[m_size]).reset_index(drop=True)
        df_cat[m_size].to_csv(
            root_dir / "ts" / "collect_indel" / f"{m_size}.csv", index=False
        )

        df_summary = (
            df_cat[m_size]
            .query("not is_wt")
            .groupby(["epoch", m_size])["count"]
            .sum()
            .reset_index()
        )
        (
            ggplot(df_summary, aes(x=m_size, fill="epoch", weight="count"))
            + geom_histogram(
                aes(y=after_stat("density")), position="dodge", binwidth=1, bins=30
            )
        ).save(root_dir / "ts" / "collect_indel" / f"{m_size}.pdf")

    return df


def original_count(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = root_dir / "ts" / "original_count"
    os.makedirs(save_dir, exist_ok=True)

    stems = []
    counts = []
    for unique_file in os.listdir(root_dir / "unique"):
        stems.append(pathlib.Path(unique_file).stem)
        counts.append(
            pd.read_csv(unique_file, names=["R1", "R2", "count"])["count"].sum()
        )

    pd.DataFrame(
        {
            "stem": stems,
            "count": counts,
        }
    ).to_csv(save_dir / "original_count.csv", index=False)


def stat_treat_full(root_dir: os.PathLike, clip: int):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = root_dir / "ts" / "stat_treat_full"
    os.makedirs(save_dir, exist_ok=True)

    df = pd.read_feather(root_dir / "analyze" / "treat" / "full" / "treat.feather")
    df = df.assign(
        del_size=lambda df: np.maximum(0, df["cut1"] - df["ref_end1"])
        + np.maximum(0, df["ref_start2"] - df["cut2"]),
        tem_ins_size=lambda df: np.maximum(0, df["cut2"] - df["ref_start2"]),
        rand_ins_size=lambda df: df["random_insertion"].str.len(),
    )
    for m_size in [
        "del_size",
        "tem_ins_size",
        "rand_ins_size",
    ]:
        df[m_size].clip(upper=clip).plot.hist(
            bins=np.linspace(0, clip + 1, clip + 2), weights=df["count"]
        ).get_figure().savefig(save_dir / f"{m_size}.pdf")
        plt.close("all")
