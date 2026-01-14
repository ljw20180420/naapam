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


def remove_duplicates(data_dir: os.PathLike, root_dir: os.PathLike):
    file_pairs = get_file_pairs(data_dir)
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "unique", exist_ok=True)
    for file_pair in file_pairs:
        R1 = file_pair["R1"]
        R2 = file_pair["R2"]
        stem = R1.name.replace(".fq.gz", "")
        with open(root_dir / "unique" / f"{stem}.unique", "w") as fd:
            subprocess.run(
                args=[
                    "removeDuplicates.sh",
                    R1.as_posix(),
                    R2.as_posix(),
                ],
                stdout=fd,
            )


def agg(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(
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


def collect_control(
    root_dir: os.PathLike,
):
    """
    R1_barcode,R1_primer,G,R1_sgRNA,R1_scaffold_prefix,R1_tail,R1_primer_score,R1_scaffold_prefix_score,R2_barcode,R2_primer,barcode_CTG_target_prefix,R2_sgRNA,pam,target_suffix,C,R2_scaffold_prefix,R2_tail,R2_primer_score,R2_scaffold_prefix_score,R2_sgRNA_score,count
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "full", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_controls = []
        for parse_file in os.listdir(root_dir / "parse"):
            if (
                utils.infer_cas(parse_file) != "control"
                or utils.infer_chip(parse_file) != chip
            ):
                continue
            df_controls.append(
                pd.read_csv(
                    root_dir / "parse" / parse_file,
                    sep="\t",
                    header=0,
                    keep_default_na=False,
                ).drop(columns=["R1_barcode", "R1_barcode"])
            )

        df_control = agg(pd.concat(df_controls))
        df_control.to_feather(root_dir / "control" / "full" / f"{chip}.feather")


def stat(df: pd.DataFrame, save_dir: os.PathLike):
    save_dir = pathlib.Path(os.fspath(save_dir))
    os.makedirs(save_dir, exist_ok=True)

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
        df[column].str.len().plot.hist(
            bins=150, weights=df["count"]
        ).get_figure().savefig(save_dir / f"{column}_length.pdf")
        plt.close("all")

    for column in [
        "R1_primer_score",
        "R1_scaffold_prefix_score",
        "R2_primer_score",
        "R2_sgRNA_score",
        "R2_scaffold_prefix_score",
    ]:
        df[column].plot.hist(bins=300, weights=df["count"]).get_figure().savefig(
            save_dir / f"{column}.pdf"
        )
        plt.close("all")

    for column in ["G", "pam", "C"]:
        df.groupby(column)["count"].sum().plot.bar().get_figure().savefig(
            save_dir / f"{column}.pdf"
        )
        plt.close("all")

    df.query("count <= 100")["count"].plot.hist(bins=100).get_figure().savefig(
        save_dir / "count.pdf"
    )
    plt.close("all")


def stat_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        stat(df=df_control, save_dir=f"figures/align/control/{chip}")


def filter_control(
    root_dir: os.PathLike,
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
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "filter", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        pam_tail = "AA" if chip in ["a1", "a2", "a3"] else "GG"
        df_control = (
            pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
            .query(
                """
                    R1_primer_score >= @min_R1_primer_score and \
                    R2_primer_score >= @min_R2_primer_score and \
                    pam.str.len() == 3 and \
                    pam.str.endswith(@pam_tail) and \
                    G == "G" and \
                    C == "C" and \
                    count >= @min_count
                """
            )
            .reset_index(drop=True)
        )
        df_control.to_feather(root_dir / "control" / "filter" / f"{chip}.feather")


def collect_treat(root_dir: os.PathLike):
    """
    R1_barcode,R1_primer,G,R1_sgRNA,R1_scaffold_prefix,R1_tail,R1_primer_score,R1_scaffold_prefix_score,R2_barcode,R2_primer,barcode_CTG_target_prefix,R2_sgRNA,pam,target_suffix,C,R2_scaffold_prefix,R2_tail,R2_primer_score,R2_scaffold_prefix_score,R2_sgRNA_score,count
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "treat" / "full", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse"):
        if utils.infer_cas(parse_file) == "control":
            continue
        df_treat = pd.read_csv(
            root_dir / "parse" / parse_file, sep="\t", header=0, keep_default_na=False
        ).drop(columns=["R1_barcode", "R1_barcode"])

        df_treat = agg(df_treat)
        df_treat.to_feather(root_dir / "treat" / "full" / f"{parse_file.stem}.feather")


def stat_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path(f"figures/align/treat")
    os.makedirs(save_dir, exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        df_treat = pd.read_feather(treat_file)
        stat(
            df=df_treat, save_dir=f"figures/align/treat/{pathlib.Path(treat_file).stem}"
        )


def filter_treat(
    root_dir: os.PathLike, min_R1_primer_score: int, min_R2_primer_score: int
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "treat" / "filter", exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        df_treat = (
            pd.read_feather(treat_file)
            .query(
                """
                    R1_primer_score >= @min_R1_primer_score and \
                    R2_primer_score >= @min_R2_primer_score and \
                    G == "G" and \
                    C == "C"
                """
            )
            .reset_index(drop=True)
        )
        df_treat.to_feather(root_dir / "treat" / "filter" / treat_file)


def reference(
    root_dir: os.PathLike,
    ext: int = 10,
):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "ref" / "ref", exist_ok=True)
    os.makedirs(root_dir / "ref" / "control", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_ref = (
            pd.read_feather(root_dir / "control" / "filter" / f"{chip}.feather")
            .assign(
                ref=lambda df: df["R2_primer"]
                + df["barcode_CTG_target_prefix"]
                + df["R2_sgRNA"]
                + df["pam"]
                + df["target_suffix"]
                + df["C"]
                + df["R2_scaffold_prefix"]
                + df["R2_tail"],
                cut=lambda df: df["R2_primer"].str.len()
                + df["barcode_CTG_target_prefix"].str.len()
                + df["R2_sgRNA"].str.len()
                - 3,
            )
            .query(
                """
                    cut + @ext <= ref.str.len() and \
                    cut >= @ext
                """
            )
            .reset_index(drop=True)
        )
        df_ref.assign(
            ref1=lambda df: df.apply(
                lambda row: row["ref"][: row["cut"] + ext], axis=1
            ),
            ref2=lambda df: df.apply(
                lambda row: row["ref"][row["cut"] - ext :], axis=1
            ),
            zero=0,
            ext=ext,
            ref2len=lambda df: df["ref2"].str.len(),
        )[["zero", "ref1", "cut", "ext", "ref2", "ref2len"]].to_csv(
            root_dir / "ref" / "ref" / f"{chip}.ref", sep="\t", header=False
        )
        df_ref.drop(columns=["ref", "cut"]).to_feather(
            root_dir / "ref" / "control" / f"{chip}.feather"
        )


def demultiplex(
    root_dir: os.PathLike,
):
    """
    |<p5barcode|9-18><p5primer|19>G<sgRNA|20><scaffold|82-92>G<16|CCN|sgRNARC|5>CAG<barcodeRC|18><p7primerRC|21><p7barcodeRC|9-18>|
    """
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "query", exist_ok=True)
    for treat_filter_file in os.listdir(root_dir / "treat" / "filter"):
        chip = utils.infer_chip(treat_filter_file)
        save_dir = pathlib.Path(f"figures/align/demultiplex/{chip}")
        os.makedirs(save_dir, exist_ok=True)

        df_ref = pd.read_feather(
            root_dir / "ref" / "control" / f"{chip}.feather"
        ).reset_index(names="ref_id")

        on = [
            "R1_primer",
            "G",
            "R1_sgRNA",
            "R1_scaffold_prefix",
            "R2_primer",
            "C",
            "R2_scaffold_prefix",
        ]
        df_query = (
            pd.read_feather(treat_filter_file)
            .assign(
                query=lambda df: df["R2_primer"]
                + df["barcode_CTG_target_prefix"]
                + df["R2_sgRNA"]
                + df["pam"]
                + df["target_suffix"]
                + df["C"]
                + df["R2_scaffold_prefix"]
                + df["R2_tail"]
            )[on + ["query", "count"]]
            .merge(
                right=df_ref[on + ["ref_id"]],
                how="left",
                on=on,
            )
        )

        df_demultiplex = df_query.groupby("query").agg(
            ref_num=pd.NamedAgg(column="query", aggfunc="size"),
            ref_id=pd.NamedAgg(column="ref_id", aggfunc="first"),
            count=pd.NamedAgg(column="count", aggfunc="first"),
        )
        df_demultiplex.loc[df_demultiplex["ref_id"].isna(), "ref_num"] = 0
        df_demultiplex["ref_num"].plot.hist(
            bins=100, weights=df_demultiplex["count"]
        ).get_figure().savefig(save_dir / "ref_num.pdf")

        df_query.query("not ref_id.isna()").reset_index(drop=True).astype(
            {"ref_id": int}
        )[["query", "count", "ref_id"]].to_csv(
            root_dir / "query" / f"{treat_filter_file.stem}.query",
            sep="\t",
            header=False,
        )
