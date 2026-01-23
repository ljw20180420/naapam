import os
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Seq


def infer_cas(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    if (
        name.startswith("a2")
        or name.startswith("a7")
        or name.startswith("d2")
        or name.startswith("h73")
        or name.startswith("h80")
        or name.startswith("p79")
        or name.startswith("p125")
        or name.startswith("q101")
        or name.startswith("q102")
        or name.startswith("q192")
    ):
        return "spycas9"
    if name.startswith("36t") or name.startswith("b2") or name.startswith("x"):
        return "spymac"
    if name.startswith("i10t") or name.startswith("i83"):
        return "ispymac"
    if name.startswith("wt"):
        return "control"
    raise Exception("unknown cas")


def infer_polq(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    if (
        name.startswith("a2")
        or name.startswith("a7")
        or name.startswith("d2")
        or name.startswith("36t")
        or name.startswith("b2")
        or name.startswith("x")
        or name.startswith("i10t")
        or name.startswith("i83")
        or name.startswith("wt")
    ):
        return "wt"
    if name.startswith("h73") or name.startswith("h80"):
        return "hel"
    if name.startswith("p79") or name.startswith("p125"):
        return "pol"
    if name.startswith("q101") or name.startswith("q102") or name.startswith("q192"):
        return "ko"
    raise Exception("unknown polq")


def infer_chip(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"-(a1|a2|a3|g1n|g2n|g3n)-", name).group(1)


def infer_time(file: os.PathLike) -> int:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return int(re.search(r"-(?:a1|a2|a3|g1n|g2n|g3n)-(\d)", name).group(1))


def infer_sample(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"^(.+?)-", name).group(1)


def infer_rep(file: os.PathLike) -> int:
    name = pathlib.Path(os.fspath(file)).name.lower()
    if re.search(r"re21-", name):
        return 3
    if (
        name.startswith("d21-")
        or name.startswith("a71-")
        or name.startswith("wt11-")
        or name.startswith("wt21-")
        or re.search(r"re2-", name)
    ):
        return 2
    return 1


def p5primer() -> str:
    return "GTGGAAAGGACGAAACACC"


def p7primer() -> str:
    return str(Seq.Seq("AATTCGCTAGCTAGGTCTTGA").reverse_complement())


def normal_scaffold() -> str:
    return "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"


def nbt_scaffold() -> str:
    return "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"


def scaffold(chip: str) -> str:
    if chip in ["a1", "a2", "a3"]:
        return nbt_scaffold()
    else:
        return normal_scaffold()


def scaffold_alt(chip: str) -> str:
    if chip in ["a1", "a2", "a3"]:
        return normal_scaffold()
    else:
        return nbt_scaffold()


def RCscaffold(chip: str) -> str:
    return str(Seq.Seq(scaffold(chip)).reverse_complement())


def RCscaffold_alt(chip: str) -> str:
    return str(Seq.Seq(scaffold_alt(chip)).reverse_complement())


def rev_comp(seq: str) -> str:
    return str(Seq.Seq(seq).reverse_complement())


#################################################
# main
#################################################


def up_del_size(df: pd.DataFrame) -> pd.Series:
    up_del_size: pd.Series = np.maximum(df["cut1"] - df["ref_end1"], 0)
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    up_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "up_del_size.pdf"
    )
    plt.close("all")
    return up_del_size


def down_del_size(df: pd.DataFrame) -> pd.Series:
    down_del_size: pd.Series = np.maximum(df["ref_start2"] - df["cut2"], 0)
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    down_del_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "down_del_size.pdf"
    )
    plt.close("all")
    return down_del_size


def rand_ins_size(df: pd.DataFrame) -> pd.Series:
    rand_ins_size: pd.Series = df["random_insertion"].str.len()
    save_dir = pathlib.Path(f"figures/hists")
    os.makedirs(save_dir, exist_ok=True)
    rand_ins_size.plot.hist(bins=30, weights=df["count"]).get_figure().savefig(
        save_dir / "rand_ins_size.pdf"
    )
    plt.close("all")
    return rand_ins_size


def is_wt(df: pd.DataFrame) -> pd.Series:
    return (
        (df["ref_end1"] == df["cut1"])
        & (df["ref_start2"] == df["cut2"])
        & (df["random_insertion"] == "")
    )


def freq_mutant(df: pd.DataFrame) -> pd.Series:
    freq_mutant = df["count"] / df.groupby(["stem", "ref_id"])["count"].transform("sum")
    freq_mutant.loc[is_wt(df)] = float("nan")

    return freq_mutant


def count_wt(df: pd.DataFrame) -> pd.Series:
    return (
        df[["stem", "ref_id"]]
        .merge(
            right=df.loc[is_wt(df)][["stem", "ref_id", "count"]],
            how="left",
            on=["stem", "ref_id"],
            validate="many_to_one",
        )["count"]
        .fillna(0)
    )


def freq_nowt(df: pd.DataFrame) -> pd.Series:
    freq_nowt = 1 - count_wt(df) / df.groupby(["stem", "ref_id"])["count"].transform(
        "sum"
    )

    return freq_nowt


def count_temN(df: pd.DataFrame, tem: int) -> pd.Series:
    return (
        df[["stem", "ref_id"]]
        .merge(
            right=df.query(
                'ref_end1 + @tem + 1 == cut1 and ref_start2 + @tem == cut2 and random_insertion == ""'
            )[["stem", "ref_id", "count"]],
            how="left",
            on=["stem", "ref_id"],
            validate="many_to_one",
        )["count"]
        .fillna(0)
    )


def freq_temN(df: pd.DataFrame, tem: int) -> pd.Series:
    freq_temN = count_temN(df, tem) / df.groupby(["stem", "ref_id"])["count"].transform(
        "sum"
    )

    return freq_temN
