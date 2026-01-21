import os
import pathlib
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import utils


def get_file_pairs(data_dir: os.PathLike) -> list[dict[str, os.PathLike]]:
    data_dir = pathlib.Path(data_dir)
    file_pairs = []
    for file in os.listdir(data_dir):
        if not file.endswith(".R2.fq.gz"):
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
                ).drop(columns=["R1_barcode", "R2_barcode"])
            )

        df_control = agg(pd.concat(df_controls))
        df_control.to_feather(root_dir / "control" / "full" / f"{chip}.feather")


def stat(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    stats = {}
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
        "pam",
    ]:
        stats[f"{column}_length"] = (
            df.assign(**{f"{column}_length": lambda df: df[column].str.len()})
            .groupby(f"{column}_length")["count"]
            .sum()
            .reset_index()
        )

    for column in [
        "R1_primer_score",
        "R1_scaffold_prefix_score",
        "R2_primer_score",
        "R2_sgRNA_score",
        "R2_scaffold_prefix_score",
        "G",
        "C",
        "pam_tail",
    ]:
        if column == "pam_tail":
            df = df.assign(pam_tail=lambda df: df["pam"].str.slice(start=-2))
        stats[column] = df.groupby(column)["count"].sum().reset_index()

    stats["count_full"] = (
        df.groupby("count")
        .size()
        .reset_index()
        .rename(columns={"count": "count_full", 0: "count"})
    )
    stats["count_small"] = (
        df.query("count <= 100").groupby("count").size().reset_index()
    ).rename(columns={"count": "count_small", 0: "count"})

    return stats


def draw(stats: dict[str, pd.DataFrame], save_dir: os.PathLike):
    save_dir = pathlib.Path(os.fspath(save_dir))
    os.makedirs(save_dir, exist_ok=True)
    for column, df in stats.items():
        if column in ["G", "C", "pam_tail"]:
            df.set_index(column)["count"].plot.bar().get_figure().savefig(
                save_dir / f"{column}.pdf"
            )
        else:
            bins = 300 if column.endswith("score") else 150
            df[column].plot.hist(bins=bins, weights=df["count"]).get_figure().savefig(
                save_dir / f"{column}.pdf"
            )

        plt.close("all")


def stat_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        stats = stat(df=df_control)
        draw(stats, save_dir=f"figures/align/stat/control/{chip}")


def filter_control(
    root_dir: os.PathLike,
    range_R1_primer_length: list[int],
    range_R1_sgRNA_length: list[int],
    min_R1_scaffold_prefix_length: int,
    max_R1_tail_length: int,
    range_R2_primer_length: list[int],
    range_barcode_CTG_target_prefix_length: int,
    range_R2_sgRNA_length: list[int],
    range_target_suffix_length: list[int],
    min_R2_scaffold_prefix_length: int,
    max_R2_tail_length: int,
    range_pam_length: int,
    min_R1_primer_score: int,
    min_R1_scaffold_prefix_score: int,
    min_R2_primer_score: int,
    min_R2_sgRNA_score: int,
    min_R2_scaffold_prefix_score: int,
    G: list[str],
    C: list[str],
    a_pam_tail: list[str],
    g_pam_tail: list[str],
    min_count: int,
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "filter", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/filter/control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_stat = pd.DataFrame(columns=["row_num", "count"], index=["full", "filter"])

        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        df_stat.loc["full", "row_num"] = df_control.shape[0]
        df_stat.loc["full", "count"] = df_control["count"].sum()

        pam_tail = a_pam_tail if chip in ["a1", "a2", "a3"] else g_pam_tail
        df_control = df_control.query(
            """
                @range_R1_primer_length[0] <= R1_primer.str.len() <= @range_R1_primer_length[1] and \
                @range_R1_sgRNA_length[0] <= R1_sgRNA.str.len() <= @range_R1_sgRNA_length[1] and \
                R1_scaffold_prefix.str.len() >= @min_R1_scaffold_prefix_length and \
                R1_tail.str.len() <= @max_R1_tail_length and \
                @range_R2_primer_length[0] <= R2_primer.str.len() <= @range_R2_primer_length[1] and \
                @range_barcode_CTG_target_prefix_length[0] <= barcode_CTG_target_prefix.str.len() <= @range_barcode_CTG_target_prefix_length[1] and \
                @range_R2_sgRNA_length[0] <= R2_sgRNA.str.len() <= @range_R2_sgRNA_length[1] and \
                @range_target_suffix_length[0] <= target_suffix.str.len() <= @range_target_suffix_length[1] and \
                R2_scaffold_prefix.str.len() >= @min_R2_scaffold_prefix_length and \
                R2_tail.str.len() <= @max_R2_tail_length and \
                @range_pam_length[0] <= pam.str.len() <= @range_pam_length[1] and \
                R1_primer_score >= @min_R1_primer_score and \
                R1_scaffold_prefix_score >= @min_R1_scaffold_prefix_score and \
                R2_primer_score >= @min_R2_primer_score and \
                R2_sgRNA_score >= @min_R2_sgRNA_score and \
                R2_scaffold_prefix_score >= @min_R2_scaffold_prefix_score and \
                G.isin(@G) and \
                C.isin(@C) and \
                pam.str.slice(start=-2).isin(@pam_tail) and \
                count >= @min_count
            """
        ).reset_index(drop=True)
        df_stat.loc["filter", "row_num"] = df_control.shape[0]
        df_stat.loc["filter", "count"] = df_control["count"].sum()

        df_stat["row_num"].plot.bar().get_figure().savefig(save_dir / "row_num.pdf")
        df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")

        df_control.to_feather(root_dir / "control" / "filter" / f"{chip}.feather")
        del df_control


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
        ).drop(columns=["R1_barcode", "R2_barcode"])

        df_treat = agg(df_treat)
        df_treat.to_feather(
            root_dir / "treat" / "full" / f"{pathlib.Path(parse_file).stem}.feather"
        )


def stat_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    stats_collect = {}
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        df_treat = pd.read_feather(root_dir / "treat" / "full" / treat_file)
        for column, df in stat(df=df_treat).items():
            if column not in stats_collect:
                stats_collect[column] = []
            stats_collect[column].append(df)

    for column, dfs in stats_collect.items():
        stats_collect[column] = (
            pd.concat(dfs).groupby(column)["count"].sum().reset_index()
        )
    draw(stats_collect, save_dir=f"figures/align/stat/treat")


def filter_treat(
    root_dir: os.PathLike,
    range_R1_primer_length: list[int],
    range_R1_sgRNA_length: list[int],
    # min_R1_scaffold_prefix_length: int,
    # max_R1_tail_length: int,
    range_R2_primer_length: list[int],
    # range_barcode_CTG_target_prefix_length: int,
    # range_R2_sgRNA_length: list[int],
    # range_target_suffix_length: list[int],
    # min_R2_scaffold_prefix_length: int,
    # max_R2_tail_length: int,
    # range_pam_length: int,
    min_R1_primer_score: int,
    # min_R1_scaffold_prefix_score: int,
    min_R2_primer_score: int,
    # min_R2_sgRNA_score: int,
    # min_R2_scaffold_prefix_score: int,
    G: list[str],
    # C: list[str],
    # a_pam_tail: list[str],
    # g_pam_tail: list[str],
    min_count: int,
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "treat" / "filter", exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        save_dir = pathlib.Path(
            f"figures/align/filter/treat/{pathlib.Path(treat_file).stem}"
        )
        os.makedirs(save_dir, exist_ok=True)
        df_stat = pd.DataFrame(columns=["row_num", "count"], index=["full", "filter"])

        df_treat = pd.read_feather(root_dir / "treat" / "full" / treat_file)
        df_stat.loc["full", "row_num"] = df_treat.shape[0]
        df_stat.loc["full", "count"] = df_treat["count"].sum()

        df_treat = df_treat.query(
            """
                @range_R1_primer_length[0] <= R1_primer.str.len() <= @range_R1_primer_length[1] and \
                @range_R1_sgRNA_length[0] <= R1_sgRNA.str.len() <= @range_R1_sgRNA_length[1] and \
                @range_R2_primer_length[0] <= R2_primer.str.len() <= @range_R2_primer_length[1] and \
                R1_primer_score >= @min_R1_primer_score and \
                R2_primer_score >= @min_R2_primer_score and \
                G.isin(@G) and \
                count >= @min_count
            """
        ).reset_index(drop=True)
        df_stat.loc["filter", "row_num"] = df_treat.shape[0]
        df_stat.loc["filter", "count"] = df_treat["count"].sum()

        df_stat["row_num"].plot.bar().get_figure().savefig(save_dir / "row_num.pdf")
        df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")

        df_treat.to_feather(root_dir / "treat" / "filter" / treat_file)


def reference(
    root_dir: os.PathLike,
    ext: int = 10,
):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "ref", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_ref = pd.read_feather(root_dir / "control" / "filter" / f"{chip}.feather")
        assert (
            df_ref["R2_sgRNA"].str.len() >= 3
        ).all(), "R2_sgRNA too short (less than 3bp)"
        df_ref = df_ref.assign(
            ref1=lambda df: df["R2_primer"]
            + df["barcode_CTG_target_prefix"]
            + df["R2_sgRNA"].str.slice(stop=-3),
            ref2=lambda df: df["R2_sgRNA"].str.slice(start=-3)
            + df["pam"]
            + df["target_suffix"]
            + df["C"]
            + df["R2_scaffold_prefix"]
            + df["R2_tail"],
            cut=lambda df: df["ref1"].str.len(),
            zero=0,
            ext=ext,
            ref2len=lambda df: df["ref2"].str.len() + ext,
        )
        assert (df_ref["cut"] >= ext).all(), f"ref1 too short (less than {ext})"
        assert (df_ref["ref2len"] >= 2 * ext).all(), f"ref2 too short (less than {ext})"
        df_ref = df_ref.assign(
            ref1=lambda df: df["ref1"] + df["ref2"].str.slice(stop=ext),
            ref2=lambda df: df["ref1"].str.slice(start=-2 * ext, stop=-ext)
            + df["ref2"],
        )
        df_ref[["zero", "ref1", "cut", "ext", "ref2", "ref2len"]].to_csv(
            root_dir / "ref" / f"{chip}.ref", sep="\t", header=False, index=False
        )


def demultiplex(
    root_dir: os.PathLike,
):
    """
    |<p5barcode|9-18><p5primer|19>G<sgRNA|20><scaffold|82-92>G<16|CCN|sgRNARC|5>CAG<barcodeRC|18><p7primerRC|21><p7barcodeRC|9-18>|
    """
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "query", exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "filter"):
        chip = utils.infer_chip(treat_file)
        save_dir = pathlib.Path(
            f"figures/align/demultiplex/{chip}/{pathlib.Path(treat_file).stem}"
        )
        os.makedirs(save_dir, exist_ok=True)

        df_ref = pd.read_feather(
            root_dir / "control" / "filter" / f"{chip}.feather"
        ).reset_index(names="ref_id")

        on = [
            "R1_primer",
            "G",
            "R1_sgRNA",
            "R1_scaffold_prefix",
            "R2_primer",
            "C",
        ]
        df_query = (
            pd.read_feather(root_dir / "treat" / "filter" / treat_file)
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

        # Get ref number distribution for each query.
        df_ref_num = (
            df_query.groupby("query")
            .agg(
                ref_num=pd.NamedAgg(column="query", aggfunc="size"),
                ref_id=pd.NamedAgg(column="ref_id", aggfunc="first"),
                count=pd.NamedAgg(column="count", aggfunc="first"),
            )
            .reset_index(drop=True)
        )
        df_ref_num.loc[df_ref_num["ref_id"].isna(), "ref_num"] = 0

        df_ref_num.query("ref_num == 0")["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302), logy=True
        ).get_figure().savefig(save_dir / "vanish_count.pdf")
        plt.close("all")

        df_ref_num = df_ref_num.groupby("ref_num")["count"].sum().reset_index()
        df_ref_num.to_csv(save_dir / "ref_num.csv", index=False)
        df_ref_num["ref_num"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302), weights=df_ref_num["count"], logy=True
        ).get_figure().savefig(save_dir / "ref_num.pdf")
        plt.close("all")

        # Get query count distribution for each ref. Queries with multiple refs are distributed to each ref evenly.
        df_count = (
            df_query.assign(
                ref_num=lambda df: df.groupby("query").transform("size"),
                count=lambda df: df["count"] / df["ref_num"],
                ref_id=lambda df: df["ref_id"].fillna(-1),
            )
            .groupby("ref_id")["count"]
            .sum()
            .reset_index()
        )
        df_count.to_csv(save_dir / "count.csv", index=False)
        # df_count.set_index("ref_id")["count"].plot.bar().get_figure().savefig(
        #     save_dir / "count.pdf"
        # )
        # plt.close("all")

        df_query.query("not ref_id.isna()").reset_index(drop=True).astype(
            {"ref_id": int}
        )[["query", "count", "ref_id"]].to_csv(
            root_dir / "query" / f"{pathlib.Path(treat_file).stem}.query",
            sep="\t",
            header=False,
            index=False,
        )
