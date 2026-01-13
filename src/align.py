import os
import pathlib
import subprocess

import matplotlib.pyplot as plt
import pandas as pd

from . import utils


def get_file_pairs(data_dir: os.PathLike) -> list[dict[str, os.PathLike]]:
    data_dir = pathlib.Path(data_dir)
    file_pairs = []
    for file in os.listdir(data_dir):
        if file.startswith("p-") or not file.endswith(".R2.fq.gz"):
            continue
        stem = file.replace(".R2.fq.gz", "")
        file_pairs.append(
            {"R1": data_dir / f"{stem}.fq.gz", "R2": data_dir / f"{stem}.R2.fq.gz"}
        )
    return file_pairs


def remove_duplicates(
    file_pairs: list[dict[str, os.PathLike]], unique_dir: os.PathLike
):
    unique_dir = pathlib.Path(os.fspath(unique_dir))
    os.makedirs(unique_dir, exist_ok=True)
    for file_pair in file_pairs:
        R1 = file_pair["R1"]
        R2 = file_pair["R2"]
        stem = R1.name.replace(".fq.gz", "")
        with open(unique_dir / f"{stem}.unique", "w") as fd:
            subprocess.run(
                args=[
                    "removeDuplicates.sh",
                    R1.as_posix(),
                    R2.as_posix(),
                ],
                stdout=fd,
            )


def collect_control(
    parse_dir: os.PathLike,
    control_dir: os.PathLike,
):
    """
    R1_barcode,R1_primer,G,R1_sgRNA,R1_scaffold_prefix,R1_tail,R1_primer_score,R1_scaffold_prefix_score,R2_barcode,R2_primer,barcode_CTG_target_prefix,R2_sgRNA,pam,target_suffix,C,R2_scaffold_prefix,R2_tail,R2_primer_score,R2_scaffold_prefix_score,R2_sgRNA_score,count
    """
    parse_dir = pathlib.Path(os.fspath(parse_dir))
    control_dir = pathlib.Path(os.fspath(control_dir))
    os.makedirs(control_dir, exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_controls = []
        for parse_file in os.listdir(parse_dir):
            if (
                utils.infer_cas(parse_file) != "control"
                or utils.infer_chip(parse_file) != chip
            ):
                continue
            df_controls.append(
                pd.read_csv(parse_dir / parse_file, sep="\t", header=0).drop(
                    columns=["R1_barcode", "R1_barcode"]
                )
            )

        df_control = (
            pd.concat(df_controls)
            .groupby(
                [
                    "R1_primer",
                    "G",
                    "R1_sgRNA",
                    "R1_scaffold_prefix",
                    "R1_tail",
                    "R2_primer",
                    "barcode_CTG_target_prefix",
                    "R2_sgRNA",
                    "pam",
                    "target_suffix",
                    "C",
                    "R2_scaffold_prefix",
                    "R2_tail",
                ]
            )
            .agg(
                R1_primer_score=pd.NamedAgg(column="R1_primer_score", aggfunc="first"),
                R1_scaffold_prefix_score=pd.NamedAgg(
                    column="R1_scaffold_prefix_score", aggfunc="first"
                ),
                R2_primer_score=pd.NamedAgg(column="R2_primer_score", aggfunc="first"),
                R2_scaffold_prefix_score=pd.NamedAgg(
                    column="R2_scaffold_prefix_score", aggfunc="first"
                ),
                R2_sgRNA_score=pd.NamedAgg(column="R2_sgRNA_score", aggfunc="first"),
                count=pd.NamedAgg(column="count", aggfunc="sum"),
            )
            .sort_values("count", ascending=False)
            .reset_index()
        )

        df_control.to_feather(control_dir / f"{chip}.feather")


def stat_control(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(control_dir / f"{chip}.feather")

        for column in [
            "R1_primer",
            "R1_sgRNA",
            "R1_scaffold_prefix",
            "R1_tail",
            "R2_primer",
            "barcode_CTG_target_prefix",
            "R2_sgRNA",
            "target_suffix",
            "R2_scaffold_prefix",
            "R2_tail",
        ]:
            df_control[column].str.len().plot.hist(
                bins=150, weights=df_control["count"]
            ).get_figure().savefig(save_dir / f"{column}_length.pdf")
            plt.close("all")

        for column in [
            "R1_primer_score",
            "R1_scaffold_prefix_score",
            "R2_primer_score",
            "R2_sgRNA_score",
            "R2_scaffold_prefix_score",
        ]:
            df_control[column].plot.hist(
                bins=300, weights=df_control["count"]
            ).get_figure().savefig(save_dir / f"{column}.pdf")
            plt.close("all")

        for column in ["G", "pam", "C"]:
            df_control.groupby(column)["count"].sum().plot.bar().get_figure().savefig(
                save_dir / f"{column}.pdf"
            )
            plt.close("all")

        df_control.query("count <= 100")["count"].plot.hist(
            bins=100
        ).get_figure().savefig(save_dir / "count.pdf")
        plt.close("all")


def filter_control(
    control_dir: os.PathLike,
    # R1_primer,
    # R1_sgRNA,
    # R1_scaffold_prefix,
    # R1_tail,
    # R2_primer,
    # barcode_CTG_target_prefix,
    # R2_sgRNA,
    # target_suffix,
    # R2_scaffold_prefix,
    # R2_tail,
    min_R1_primer_score: int,
    # min_R1_scaffold_prefix_score: int,
    min_R2_primer_score: int,
    # min_R2_sgRNA_score: int,
    # min_R2_scaffold_prefix_score: int,
    min_count: int,
):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = (
            pd.read_feather(control_dir / f"{chip}.feather")
            .query(
                """
                    R1_primer_score >= @min_R1_primer_score and \
                    R2_primer_score >= @min_R2_primer_score and \
                    G == "G" and \
                    C == "C" and \
                    count >= @min_count
                """
            )
            .reset_index(drop=True)
        )
        df_control.to_feather(control_dir / f"{chip}_filter.feather")


def find_reference(
    treat_file: os.PathLike, control_file: os.PathLike, chip: str, prime_len: int
):
    """
    |<barcodeP5|9/18><primerP5|19>G<sgRNA|20><scaffold|83/93><16|CCN|sgRNARC|5><CAG><barcodeRC|18><primerP7RC|21><barcodeP7RC|9/18>|
    """
    #
    scaffold_len = 93 if chip.startswith("a") else 83
    p5_len = 18
    p7_len = 16
    prime_len = 21
    barcode_len = 18

    df_treat = pd.read_csv(treat_file, sep="\t", names=["R1", "R2", "count"])
    df_control = pd.read_csv(control_file, sep="\t", names=["R1", "R2", "count"])
