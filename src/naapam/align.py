import os
import pathlib
import subprocess
from importlib import resources

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
                "barcode_head",
                "barcode",
                "CTG_target_prefix",
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
            barcode_score=pd.NamedAgg(column="barcode_score", aggfunc="first"),
            R1_sgRNA_bowtie2_score=pd.NamedAgg(
                column="R1_sgRNA_bowtie2_score", aggfunc="first"
            ),
            R2_sgRNA_bowtie2_score=pd.NamedAgg(
                column="R2_sgRNA_bowtie2_score", aggfunc="first"
            ),
            barcode_id=pd.NamedAgg(column="barcode_id", aggfunc="first"),
            R1_sgRNA_id=pd.NamedAgg(column="R1_sgRNA_id", aggfunc="first"),
            R2_sgRNA_id=pd.NamedAgg(column="R2_sgRNA_id", aggfunc="first"),
            count=pd.NamedAgg(column="count", aggfunc="sum"),
        )
        .sort_values("count", ascending=False)
        .reset_index()
    )


def stat(df: pd.DataFrame, save_dir: os.PathLike):
    save_dir = pathlib.Path(os.fspath(save_dir))
    os.makedirs(save_dir, exist_ok=True)
    for column in [
        "R1_primer",
        "R1_sgRNA",
        "R1_scaffold_prefix",
        "R1_tail",
        "R2_primer",
        "barcode_head",
        "barcode",
        "CTG_target_prefix",
        "R2_sgRNA",
        "pam",
        "target_suffix",
        "R2_scaffold_prefix",
        "R2_tail",
    ]:
        df_stat = (
            df[[column, "count"]]
            .assign(**{f"{column}_length": lambda df: df[column].str.len()})
            .groupby(f"{column}_length")["count"]
            .sum()
            .reset_index()
        )
        df_stat.to_csv(save_dir / f"{column}_length.csv", index=False)

    for column in [
        "G",
        "C",
        "pam_tail",
        "R1_primer_score",
        "R1_scaffold_prefix_score",
        "R2_primer_score",
        "R2_sgRNA_score",
        "R2_scaffold_prefix_score",
        "barcode_score",
        "R1_sgRNA_bowtie2_score",
        "R2_sgRNA_bowtie2_score",
        "barcode_id",
        "R1_sgRNA_id",
        "R2_sgRNA_id",
    ]:
        if column == "pam_tail":
            df = df.assign(pam_tail=lambda df: df["pam"].str.slice(start=-2))
        df_stat = df[[column, "count"]].groupby(column)["count"].sum().reset_index()
        df_stat.to_csv(save_dir / f"{column}.csv", index=False)

    df_stat = (
        df[["count"]]
        .groupby("count")
        .size()
        .reset_index()
        .rename(columns={"count": "count_full", 0: "count"})
    )
    df_stat.to_csv(save_dir / f"count_full.csv", index=False)

    df_stat = (
        df[["count"]].clip(upper=300).groupby("count").size().reset_index()
    ).rename(columns={"count": "count_small", 0: "count"})
    df_stat.to_csv(save_dir / f"count_small.csv", index=False)


def draw(save_dirs: list[os.PathLike], summary_dir: os.PathLike):
    save_dirs = [pathlib.Path(os.fspath(save_dir)) for save_dir in save_dirs]
    summary_dir = pathlib.Path(os.fspath(summary_dir))
    os.makedirs(summary_dir, exist_ok=True)
    for csv_file in os.listdir(save_dirs[0]):
        if not csv_file.endswith(".csv"):
            continue
        df_stat = pd.concat(
            [pd.read_csv(save_dir / csv_file, header=0) for save_dir in save_dirs]
        ).reset_index(drop=True)
        name = (
            df_stat.columns[0] if df_stat.columns[0] != "count" else df_stat.columns[1]
        )
        if name in [
            "R1_primer_length",
            "R2_primer_length",
            "barcode_id",
            "R1_sgRNA_id",
            "R2_sgRNA_id",
            "count_full",
            "count_small",
        ]:
            logy = True
        else:
            logy = False

        if name in ["G", "C", "pam_tail"]:
            df_stat.groupby(name)["count"].sum().plot.bar(
                logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        elif name in ["barcode_id", "R1_sgRNA_id", "R2_sgRNA_id"]:
            df_stat.groupby(name)["count"].sum().reset_index().plot.scatter(
                x=name, y="count", logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        elif name.endswith("_length") or name in ["count_full", "count_small"]:
            bins = 150 if name.endswith("_length") else 300
            df_stat[name].clip(upper=bins).plot.hist(
                bins=np.linspace(0, bins + 1, bins + 2),
                weights=df_stat["count"],
                logy=logy,
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        elif name.endswith("_score"):
            bins = 300
            df_stat[name].plot.hist(
                bins=bins, weights=df_stat["count"], logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        else:
            raise ValueError("Unknown column")
        plt.close("all")


def collect_control(
    root_dir: os.PathLike,
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "full", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_controls = []
        for parse_file in os.listdir(root_dir / "parse" / "sgRNA"):
            if (
                utils.infer_cas(parse_file) != "control"
                or utils.infer_chip(parse_file) != chip
            ):
                continue
            df_controls.append(
                pd.read_csv(
                    root_dir / "parse" / "sgRNA" / parse_file,
                    sep="\t",
                    header=0,
                    keep_default_na=False,
                ).drop(columns=["R1_barcode", "R2_barcode"])
            )

        df_controls = pd.concat(df_controls)
        df_controls = agg(df_controls)
        df_controls.to_feather(root_dir / "control" / "full" / f"{chip}.feather")


def stat_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/stat_control/{chip}")
        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        stat(df=df_control, save_dir=save_dir)
        draw(save_dirs=[save_dir], summary_dir=save_dir)
        del df_control


def filter_nofunc_control(
    root_dir: os.PathLike,
    range_R1_primer_length: list[int],
    range_R1_sgRNA_length: list[int],
    min_R1_scaffold_prefix_length: int,
    max_R1_tail_length: int,
    range_R2_primer_length: list[int],
    range_barcode_head_length: int,
    range_barcode_length: int,
    range_CTG_target_prefix_length: int,
    range_R2_sgRNA_length: list[int],
    range_pam_length: int,
    range_target_suffix_length: list[int],
    min_R2_scaffold_prefix_length: int,
    max_R2_tail_length: int,
    min_R1_primer_score: int,
    min_R1_scaffold_prefix_score: int,
    min_R2_primer_score: int,
    min_R2_sgRNA_score: int,
    min_R2_scaffold_prefix_score: int,
    min_barcode_score: int,
    min_R1_sgRNA_bowtie2_score: int,
    min_R2_sgRNA_bowtie2_score: int,
    G: list[str],
    C: list[str],
    a_pam_tail: list[str],
    g_pam_tail: list[str],
    min_count: int,
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "func", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/filter_nofunc_control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_stat = pd.DataFrame(columns=["row_num", "count"], index=["full", "filter"])

        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        df_stat.loc["full", "row_num"] = df_control.shape[0]
        df_stat.loc["full", "count"] = df_control["count"].sum()

        pam_tail = a_pam_tail if chip in ["a1", "a2", "a3"] else g_pam_tail
        df_control.query(
            """
                @range_R1_primer_length[0] <= R1_primer.str.len() <= @range_R1_primer_length[1] and \
                @range_R1_sgRNA_length[0] <= R1_sgRNA.str.len() <= @range_R1_sgRNA_length[1] and \
                R1_scaffold_prefix.str.len() >= @min_R1_scaffold_prefix_length and \
                R1_tail.str.len() <= @max_R1_tail_length and \
                @range_R2_primer_length[0] <= R2_primer.str.len() <= @range_R2_primer_length[1] and \
                @range_barcode_head_length[0] <= barcode_head.str.len() <= @range_barcode_head_length[1] and \
                @range_barcode_length[0] <= barcode.str.len() <= @range_barcode_length[1] and \
                @range_CTG_target_prefix_length[0] <= CTG_target_prefix.str.len() <= @range_CTG_target_prefix_length[1] and \
                @range_R2_sgRNA_length[0] <= R2_sgRNA.str.len() <= @range_R2_sgRNA_length[1] and \
                @range_pam_length[0] <= pam.str.len() <= @range_pam_length[1] and \
                @range_target_suffix_length[0] <= target_suffix.str.len() <= @range_target_suffix_length[1] and \
                R2_scaffold_prefix.str.len() >= @min_R2_scaffold_prefix_length and \
                R2_tail.str.len() <= @max_R2_tail_length and \
                R1_primer_score >= @min_R1_primer_score and \
                R1_scaffold_prefix_score >= @min_R1_scaffold_prefix_score and \
                R2_primer_score >= @min_R2_primer_score and \
                R2_sgRNA_score >= @min_R2_sgRNA_score and \
                R2_scaffold_prefix_score >= @min_R2_scaffold_prefix_score and \
                barcode_score >= @min_barcode_score and \
                R1_sgRNA_bowtie2_score >= @min_R1_sgRNA_bowtie2_score and \
                R2_sgRNA_bowtie2_score >= @min_R2_sgRNA_bowtie2_score and \
                G.isin(@G) and \
                C.isin(@C) and \
                pam.str.slice(start=-2).isin(@pam_tail) and \
                count >= @min_count and \
                barcode_id >= 0 and \
                barcode_id == R1_sgRNA_id == R2_sgRNA_id
            """,
            inplace=True,
        )
        df_control.reset_index(drop=True, inplace=True)
        df_stat.loc["filter", "row_num"] = df_control.shape[0]
        df_stat.loc["filter", "count"] = df_control["count"].sum()

        df_stat.to_csv(save_dir / "stat.csv")
        df_stat["row_num"].plot.bar(logy=True).get_figure().savefig(
            save_dir / "row_num.pdf"
        )
        df_stat["count"].plot.bar(logy=True).get_figure().savefig(
            save_dir / "count.pdf"
        )

        df_control.to_feather(root_dir / "control" / "func" / f"{chip}.feather")


def cluster_func_control_by_mutant(
    root_dir: os.PathLike, ext: int, plasmid_file: os.PathLike | None
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "cluster", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        if plasmid_file is None:
            if chip in ["a1", "a2", "a3"]:
                plasmid_file = "final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv"
            else:
                plasmid_file = "plasmids/final_hgsgrna_libb_all_0811-NGG.csv"
        with resources.as_file(
            resources.files().parent / "plasmids" / plasmid_file
        ) as pf:
            df_plasmid = pd.read_csv(pf, header=0)
        df_plasmid = df_plasmid.assign(
            ref2=lambda df: (
                "AAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTG"
                + df["Target sequence"].str.slice(stop=22 + ext)
            ).map(utils.rev_comp),
            ref1=lambda df: (
                df["Target sequence"].str.slice(start=22 - ext)
                + "CAG"
                + df["Barcode2"]
                + df["Primer binding sites 20nt"]
            ).map(utils.rev_comp),
            cut=lambda df: df["ref1"].str.len() - ext,
            ext=ext,
        ).reset_index(names="barcode_id")

        df_control = pd.read_feather(root_dir / "control" / "func" / f"{chip}.feather")
        df_control = df_control.assign(
            query=lambda df: df["R2_primer"]
            + df["barcode_head"]
            + df["barcode"]
            + df["CTG_target_prefix"]
            + df["R2_sgRNA"]
            + df["pam"]
            + df["target_suffix"]
            + df["C"]
            + df["R2_scaffold_prefix"]
            + df["R2_tail"],
        )[["barcode_id", "query", "count"]].merge(
            right=df_plasmid[["barcode_id", "ref1", "ref2", "cut", "ext"]],
            how="left",
            on=["barcode_id"],
            validate="many_to_one",
        )

        df_alg = utils.call_rearr(
            df_control["ref1"],
            df_control["ref2"],
            df_control["cut"],
            df_control["ext"],
            df_control["query"],
            df_control["count"],
        )

        df_control = (
            df_control.assign(
                ref_end1=df_alg["ref_end1"],
                ref_start2=df_alg["ref_start2"],
                random_insertion=df_alg["random_insertion"],
                cut1=lambda df: df["cut"],
                cut2=lambda df: df["cut1"] + 2 * ext,
            )
            .groupby(
                [
                    "barcode_id",
                    "ref_end1",
                    "ref_start2",
                    "random_insertion",
                    "cut1",
                    "cut2",
                    "ref1",
                    "ref2",
                ]
            )["count"]
            .sum()
            .sort_values()
            .reset_index()
            .assign(
                percentage=lambda df: df.groupby("barcode_id")["count"].transform(
                    "cumsum"
                )
                / df.groupby("barcode_id")["count"].transform("sum"),
                rank=lambda df: df.groupby("barcode_id")["count"].rank(ascending=False),
            )
        )

        # count_wt
        df_control = (
            df_control.merge(
                right=df_control.query(
                    """
                    ref_end1 == cut1 and \
                    ref_start2 == cut2 and \
                    random_insertion == ""
                """
                )[["barcode_id", "count"]].rename(columns={"count": "count_wt"}),
                how="left",
                on=["barcode_id"],
                validate="many_to_one",
            )
            .assign(**{"count_wt": lambda df: df["count_wt"].fillna(0)})
            .astype({"count_wt": int})
        )

        # count_temN
        for tem in range(1, 5):
            df_control = (
                df_control.merge(
                    right=df_control.query(
                        """
                        ref_end1 + @tem + 1 == cut1 and \
                        ref_start2 + @tem == cut2 and \
                        random_insertion == ""
                    """
                    )[["barcode_id", "count"]].rename(
                        columns={"count": f"count_tem{tem}"}
                    ),
                    how="left",
                    on=["barcode_id"],
                    validate="many_to_one",
                )
                .assign(
                    **{f"count_tem{tem}": lambda df: df[f"count_tem{tem}"].fillna(0)}
                )
                .astype({f"count_tem{tem}": int})
            )

        # count_tot, first, second
        df_control = df_control.assign(
            count_tot=lambda df: df.groupby("barcode_id")["count"].transform("sum"),
            first=lambda df: df[
                ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
            ].max(axis=1),
            second=lambda df: df[
                ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
            ].apply(lambda row: row.nlargest(2).min(), axis=1),
        )

        # mutant_type_num
        df_control = df_control.assign(
            mutant_type_num=lambda df: df.groupby("barcode_id").transform("size")
        )

        df_control.to_feather(root_dir / "control" / "cluster" / f"{chip}.feather")


def stat_func_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/stat_func_control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(
            root_dir / "control" / "cluster" / f"{chip}.feather"
        )

        # number of mutant types per barcode_id
        df_control.groupby("barcode_id")["mutant_type_num"].first().clip(
            upper=100
        ).plot.hist(bins=np.linspace(0, 101, 102)).get_figure().savefig(
            save_dir / "mutant_type_num.pdf"
        )
        plt.close("all")

        # up_del_size
        df_stat = (
            df_control.assign(up_del_size=lambda df: utils.up_del_size(df))
            .groupby("up_del_size")["count"]
            .sum()
            .reset_index()
        )
        df_stat["up_del_size"].plot.hist(
            bins=np.linspace(0, 101, 102), weights=df_stat["count"]
        ).get_figure().savefig(save_dir / "up_del_size.pdf")
        plt.close("all")

        # down_del_size
        df_stat = (
            df_control.assign(down_del_size=lambda df: utils.down_del_size(df))
            .groupby("down_del_size")["count"]
            .sum()
            .reset_index()
        )
        df_stat["down_del_size"].plot.hist(
            bins=np.linspace(0, 101, 102), weights=df_stat["count"]
        ).get_figure().savefig(save_dir / "down_del_size.pdf")
        plt.close("all")

        # rand_ins_size
        df_stat = (
            df_control.assign(rand_ins_size=lambda df: utils.rand_ins_size(df))
            .groupby("rand_ins_size")["count"]
            .sum()
            .reset_index()
        )
        df_stat["rand_ins_size"].plot.hist(
            bins=np.linspace(0, 101, 102), weights=df_stat["count"]
        ).get_figure().savefig(save_dir / "rand_ins_size.pdf")
        plt.close("all")

        # count
        df_control["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302),
        ).get_figure().savefig(save_dir / "count.pdf")
        plt.close("all")

        # count_tot
        df_control.groupby("barcode_id")["count_tot"].first().clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302),
        ).get_figure().savefig(save_dir / "count_tot.pdf")
        plt.close("all")

        # freq_wt
        df_control.assign(freq_wt=lambda df: df["count_wt"] / df["count_tot"]).groupby(
            "barcode_id"
        )["freq_wt"].first().plot.hist(bins=100, logy=True).get_figure().savefig(
            save_dir / "freq_wt.pdf"
        )
        plt.close("all")

        # second_rel_first
        df_control.assign(
            second_rel_first=lambda df: df["second"] / (df["first"] + 1e-6)
        ).groupby("barcode_id")["second_rel_first"].first().plot.hist(
            bins=100, logy=True
        ).get_figure().savefig(
            save_dir / "second_rel_first.pdf"
        )
        plt.close("all")

        # type_max
        df_stat = (
            df_control.assign(
                type_max=lambda df: df[
                    ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
                ].idxmax(axis=1)
            )
            .groupby("barcode_id")
            .agg(
                count_tot=pd.NamedAgg(column="count", aggfunc="sum"),
                type_max=pd.NamedAgg(column="type_max", aggfunc="first"),
            )
            .reset_index()
            .groupby("type_max")
            .agg(
                num_type_max=pd.NamedAgg(column="type_max", aggfunc="size"),
                count_type_max=pd.NamedAgg(column="count_tot", aggfunc="sum"),
            )
        )
        df_stat["num_type_max"].plot.bar(logy=True).get_figure().savefig(
            save_dir / "num_type_max.pdf"
        )
        plt.close("all")
        df_stat["count_type_max"].plot.bar(logy=True).get_figure().savefig(
            save_dir / "count_type_max.pdf"
        )
        plt.close("all")


def filter_low_quality_barcode(
    root_dir: os.PathLike,
    max_mutant_type_num: int,
    min_count_tot: int,
    min_freq_wt: float,
    max_second_rel_first: float,
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "hq_bar", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/filter_low_quality_barcode/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(
            root_dir / "control" / "cluster" / f"{chip}.feather"
        )

        df_stat = pd.DataFrame(columns=["bar_num"], index=["full", "filter"])
        df_stat.loc["full", "bar_num"] = df_control["barcode_id"].unique().shape[0]

        df_control = df_control.query(
            """
                mutant_type_num <= @max_mutant_type_num and \
                (count_tot >= @min_count_tot or count_wt / count_tot >= @min_freq_wt) and \
                second / (first + 1e-6) <= @max_second_rel_first
            """
        ).reset_index(drop=True)

        df_stat.loc["filter", "bar_num"] = df_control["barcode_id"].unique().shape[0]
        df_stat["bar_num"].plot.bar().get_figure().savefig(save_dir / "bar_num.pdf")

        df_control.to_feather(root_dir / "control" / "hq_bar" / f"{chip}.feather")


def filter_low_quality_mutant(
    root_dir: os.PathLike,
    min_percentage: float,
    max_rank: int,
    max_up_del_size: int,
    max_down_del_size: int,
    max_rand_ins_size: int,
    min_count: int,
):
    assert max_down_del_size <= 3, "pam is resected"
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "hq_mut", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/filter_low_quality_mutant/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(
            root_dir / "control" / "hq_bar" / f"{chip}.feather"
        )

        df_stat = pd.DataFrame(columns=["mut_num"], index=["full", "filter"])
        df_stat.loc["full", "mut_num"] = df_control.shape[0]

        df_control = df_control.query(
            """
                percentage >= @min_percentage and \
                rank <= @max_rank and \
                cut1 - ref_end1 <= @max_up_del_size and \
                ref_start2 - cut2 <= @max_down_del_size and \
                random_insertion.str.len() <= @max_rand_ins_size and \
                count >= @min_count
            """
        ).reset_index(drop=True)

        df_stat.loc["filter", "mut_num"] = df_control.shape[0]
        df_stat["mut_num"].plot.bar().get_figure().savefig(save_dir / "mut_num.pdf")
        plt.close("all")

        df_control.groupby("barcode_id").size().rename("mutant_type_num").clip(
            upper=100
        ).plot.hist(bins=np.linspace(0, 101, 102)).get_figure().savefig(
            save_dir / "mutant_type_num.pdf"
        )

        df_control.to_feather(root_dir / "control" / "hq_mut" / f"{chip}.feather")


def generate_reference(
    root_dir: os.PathLike,
    ext: int,
):
    def generate_reference_row(row: pd.Series, ext: int) -> pd.Series:
        ref = (
            row["ref1"][: row["ref_end1"]]
            + row["random_insertion"]
            + row["ref2"][row["ref_start2"] - len(row["ref1"]) :]
        )
        ref1 = row["ref1"][: row["ref_end1"]] + row["random_insertion"]
        if row["ref_start2"] > row["cut2"]:
            assert len(ref1) >= row["ref_start2"] - row["cut2"], "upstream too short"
            ref1 = ref1[: (row["cut2"] - row["ref_start2"])]
        elif row["ref_start2"] < row["cut2"]:
            ref1 = (
                ref1
                + row["ref2"][
                    row["ref_start2"]
                    - len(row["ref1"]) : row["cut2"]
                    - len(row["ref1"])
                ]
            )
        cut = len(ref1)
        ref1 = ref[: cut + ext]
        ref2 = ref[cut - ext :]
        return pd.Series(
            {
                "ref1": ref1,
                "ref2": ref2,
                "cut": cut,
            }
        )

    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "ref", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(
            root_dir / "control" / "hq_mut" / f"{chip}.feather"
        )

        assert (
            df_control["ref_start2"] <= df_control["cut2"] + 3
        ).all(), "pam is mutated"

        df_ref = (
            df_control[
                ["ref_end1", "ref_start2", "random_insertion", "cut2", "ref1", "ref2"]
            ]
            .apply(lambda row, ext=ext: generate_reference_row(row, ext), axis=1)
            .assign(
                zero=0,
                ext=ext,
                ref2len=lambda df: df["ref2"].str.len(),
            )
        )

        df_ref[["zero", "ref1", "cut", "ext", "ref2", "ref2len"]].to_csv(
            root_dir / "ref" / f"{chip}.ref",
            sep="\t",
            header=False,
            index=False,
        )


def collect_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "treat" / "full", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse" / "sgRNA"):
        df_treat = pd.read_csv(
            root_dir / "parse" / "sgRNA" / parse_file,
            sep="\t",
            header=0,
            keep_default_na=False,
        ).drop(columns=["R1_barcode", "R2_barcode"])

        df_treat = agg(df_treat)
        df_treat.to_feather(
            root_dir / "treat" / "full" / f"{pathlib.Path(parse_file).stem}.feather"
        )


def stat_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = pathlib.Path("figures/align/stat_treat")
    save_dirs = []
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        df_treat = pd.read_feather(root_dir / "treat" / "full" / treat_file)
        stat(
            df=df_treat,
            save_dir=save_dir / pathlib.Path(treat_file).stem,
        )
        save_dirs.append(save_dir / pathlib.Path(treat_file).stem)

    draw(save_dirs, summary_dir=save_dir)


def filter_treat(
    root_dir: os.PathLike,
    range_R1_primer_length: list[int],
    range_R1_sgRNA_length: list[int],
    # min_R1_scaffold_prefix_length: int,
    # max_R1_tail_length: int,
    range_R2_primer_length: list[int],
    # range_barcode_head_length: int,
    # range_barcode_length: int,
    # range_CTG_target_prefix_length: int,
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
    # min_barcode_score: int,
    min_R1_sgRNA_bowtie2_score: int,
    # min_R2_sgRNA_bowtie2_score: int,
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
            f"figures/align/filter_treat/{pathlib.Path(treat_file).stem}"
        )
        os.makedirs(save_dir, exist_ok=True)

        df_treat = pd.read_feather(root_dir / "treat" / "full" / treat_file)

        df_stat = pd.DataFrame(columns=["row_num", "count"], index=["full", "filter"])
        df_stat.loc["full", "row_num"] = df_treat.shape[0]
        df_stat.loc["full", "count"] = df_treat["count"].sum()

        df_treat = df_treat.query(
            """
                @range_R1_primer_length[0] <= R1_primer.str.len() <= @range_R1_primer_length[1] and \
                @range_R1_sgRNA_length[0] <= R1_sgRNA.str.len() <= @range_R1_sgRNA_length[1] and \
                @range_R2_primer_length[0] <= R2_primer.str.len() <= @range_R2_primer_length[1] and \
                R1_primer_score >= @min_R1_primer_score and \
                R2_primer_score >= @min_R2_primer_score and \
                R1_sgRNA_bowtie2_score >= @min_R1_sgRNA_bowtie2_score and \
                G.isin(@G) and \
                count >= @min_count and \
                barcode_id >= 0 and \
                barcode_id == R1_sgRNA_id
            """
        ).reset_index(drop=True)
        df_stat.loc["filter", "row_num"] = df_treat.shape[0]
        df_stat.loc["filter", "count"] = df_treat["count"].sum()

        df_stat.to_csv(save_dir / "stat.csv", index=False)
        df_stat["row_num"].plot.bar().get_figure().savefig(save_dir / "row_num.pdf")
        df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")

        df_treat.to_feather(root_dir / "treat" / "filter" / treat_file)


def demultiplex(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "query" / "found", exist_ok=True)
    os.makedirs(root_dir / "query" / "not_found", exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "filter"):
        save_dir = pathlib.Path(
            f"figures/align/demultiplex/{pathlib.Path(treat_file).stem}"
        )
        os.makedirs(save_dir, exist_ok=True)

        chip = utils.infer_chip(treat_file)
        df_ref = pd.read_feather(root_dir / "control" / "hq_mut" / f"{chip}.feather")[
            ["barcode_id"]
        ].reset_index(names="ref_id")

        df_query = (
            pd.read_feather(root_dir / "treat" / "filter" / treat_file)
            .assign(
                query=lambda df: df["R2_primer"]
                + df["barcode_head"]
                + df["barcode"]
                + df["CTG_target_prefix"]
                + df["R2_sgRNA"]
                + df["pam"]
                + df["target_suffix"]
                + df["C"]
                + df["R2_scaffold_prefix"]
                + df["R2_tail"]
            )[["barcode_id", "query", "count"]]
            .merge(
                right=df_ref[["barcode_id", "ref_id"]],
                how="left",
                on=["barcode_id"],
            )
            .assign(
                ref_id=lambda df: df["ref_id"].fillna(-1).astype(int),
                count_distri=lambda df: df["count"]
                / df.groupby("query").transform("size"),
            )
        )

        # Get query number/count distribution for each ref. Count of query with multiple refs are distributed to each ref evenly.
        df_query.groupby("ref_id").size().rename("number").to_csv(
            save_dir / "query_number_per_ref.csv"
        )
        df_query.groupby("ref_id")["count_distri"].sum().rename("count").to_csv(
            save_dir / "query_count_per_ref.csv"
        )

        df_query.query("ref_id == -1")[["query", "count"]].to_csv(
            root_dir / "query" / "not_found" / f"{pathlib.Path(treat_file).stem}.query",
            sep="\t",
            header=False,
            index=False,
        )
        df_query.query("ref_id != -1")[["query", "count", "ref_id"]].to_csv(
            root_dir / "query" / "found" / f"{pathlib.Path(treat_file).stem}.query",
            sep="\t",
            header=False,
            index=False,
        )

    summary_demultiplex(save_dir="figures/align/demultiplex")


def summary_demultiplex(save_dir: os.PathLike):
    save_dir = pathlib.Path(os.fspath(save_dir))
    df_query_number_per_refs = []
    df_query_count_per_refs = []
    for stem in os.listdir(save_dir):
        if not os.path.isdir(save_dir / stem):
            continue
        df_query_number_per_refs.append(
            pd.read_csv(save_dir / stem / "query_number_per_ref.csv", header=0).assign(
                stem=stem
            )
        )
        df_query_count_per_refs.append(
            pd.read_csv(save_dir / stem / "query_count_per_ref.csv", header=0).assign(
                stem=stem
            )
        )

    pd.concat(df_query_number_per_refs).groupby("ref_id")[
        "number"
    ].sum().reset_index().plot.scatter(
        x="ref_id", y="number", logy=True
    ).get_figure().savefig(
        save_dir / "query_number_per_ref.pdf"
    )

    pd.concat(df_query_count_per_refs).groupby("ref_id")[
        "count"
    ].sum().reset_index().plot.scatter(
        x="ref_id", y="count", logy=True
    ).get_figure().savefig(
        save_dir / "query_count_per_ref.pdf"
    )
