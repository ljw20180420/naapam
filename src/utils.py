import os
import pathlib
import re
import subprocess
import tempfile
from typing import Iterable

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
    return up_del_size


def down_del_size(df: pd.DataFrame) -> pd.Series:
    down_del_size: pd.Series = np.maximum(df["ref_start2"] - df["cut2"], 0)
    return down_del_size


def rand_ins_size(df: pd.DataFrame) -> pd.Series:
    rand_ins_size: pd.Series = df["random_insertion"].str.len()
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
                """
                    ref_end1 + @tem + 1 == cut1 and \
                    ref_start2 + @tem == cut2 and \
                    random_insertion == ""
                """
            )[["stem", "ref_id", "count"]],
            how="left",
            on=["stem", "ref_id"],
            validate="many_to_one",
        )["count"]
        .fillna(0)
    )


def count_temN_blunt(df: pd.DataFrame, tem: int) -> pd.Series:
    return (
        df[["stem", "ref_id"]]
        .merge(
            right=df.query(
                'ref_end1 == cut1 and ref_start2 + @tem == cut2 and random_insertion == ""'
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


def freq_temN_blunt(df: pd.DataFrame, tem: int) -> pd.Series:
    freq_temN_blunt = count_temN_blunt(df, tem) / df.groupby(["stem", "ref_id"])[
        "count"
    ].transform("sum")

    return freq_temN_blunt


def freq_temN_dummy_rel_blunt(df: pd.DataFrame, tem: int) -> pd.Series:
    freq_temN_dummy_rel_blunt = count_temN(df, tem) / (count_temN_blunt(df, tem) + 1e-6)

    return freq_temN_dummy_rel_blunt


###########################################
# rearr
###########################################


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


def call_rearr(
    ref1s: Iterable[str],
    ref2s: Iterable[str],
    cuts: Iterable[int],
    exts: Iterable[int],
    queries: Iterable[str],
    counts: Iterable[int],
) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(os.fspath(tmpdir))
        with open(tmpdir / "ref", "w") as fd_ref, open(
            tmpdir / "query", "w"
        ) as fd_query, open(tmpdir / "correct", "w") as fd_correct:
            for ref_id, (ref1, ref2, cut, ext, query, count) in enumerate(
                zip(ref1s, ref2s, cuts, exts, queries, counts)
            ):
                fd_ref.write(f"{0}\t{ref1}\t{cut}\t{ext}\t{ref2}\t{len(ref2)}\n")
                fd_query.write(f"{query}\t{count}\t{ref_id}\n")
                fd_correct.write(f"up\n")

        subprocess.run(
            args=" ".join(
                [
                    "rearrangement",
                    "<",
                    (tmpdir / "query").as_posix(),
                    "3<",
                    (tmpdir / "ref").as_posix(),
                    "|",
                    "gawk",
                    "-f",
                    "correct_micro_homology.awk",
                    "--",
                    (tmpdir / "ref").as_posix(),
                    (tmpdir / "correct").as_posix(),
                    ">",
                    (tmpdir / "alg").as_posix(),
                ]
            ),
            shell=True,
        )

        df_alg = read_alg(tmpdir / "alg")

    return df_alg


def infer_cut(row: pd.Series, ext: int) -> pd.Series:
    ref_line = row["ref"]
    query_line = row["query"]
    ref_end1 = row["ref_end1"]
    ref_start2 = row["ref_start2"]
    cut = row["cut_query"]
    first_lower_end = re.search(r"[acgtn]", ref_line).span()[1]
    before_rand_ins = (
        re.search(r"[acgtn]", ref_line[first_lower_end:]).span()[1] + first_lower_end
    )
    after_rand_ins = (
        re.search(r"[acgtn]", ref_line[before_rand_ins:]).span()[0] + before_rand_ins
    )
    ref_line_array = np.array(list(ref_line))
    temp = (ref_line_array[:before_rand_ins] != "-").cumsum()
    ref_end1_alg = np.where(temp == ref_end1)[0].item() + 1
    ref1_len = temp[-1]
    ref_start2_alg = (
        np.where(
            (ref_line_array[after_rand_ins:] != "-").cumsum() == ref_start2 - ref1_len
        )[0].item()
        + 1
        + after_rand_ins
    )

    query_line_array = np.array(list(query_line))
    cut_alg = np.where((query_line_array != "-").cumsum() == cut)[0].item() + 1
    if cut_alg <= ref_end1_alg:
        ref1_correct = ref_line[:cut_alg]
    elif cut_alg <= after_rand_ins:
        ref1_correct = ref_line[:ref_end1_alg] + query_line[before_rand_ins:cut_alg]
    else:
        ref1_correct = (
            ref_line[:ref_end1_alg]
            + query_line[before_rand_ins:after_rand_ins]
            + ref_line[ref_start2_alg:cut_alg]
        )

    cut_correct = len(ref1_correct.replace("-", ""))
    ref_correct = (
        (
            ref_line[:ref_end1_alg]
            + query_line[before_rand_ins:after_rand_ins]
            + ref_line[ref_start2_alg:]
        )
        .replace("-", "")
        .upper()
    )
    assert cut_correct >= ext, "ref1 too short"
    assert len(ref_correct) - cut_correct >= ext, "ref2 too short"
    ref1_correct = ref_correct[: cut_correct + ext]
    ref2_correct = ref_correct[cut_correct - ext :]

    breakpoint()

    return pd.Series(
        {
            "ref1": ref1_correct,
            "ref2": ref2_correct,
            "cut": cut_correct,
        }
    )
