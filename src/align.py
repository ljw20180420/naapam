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
            barcode_id=pd.NamedAgg(column="barcode_id", aggfunc="first"),
            count=pd.NamedAgg(column="count", aggfunc="sum"),
        )
        .sort_values("count", ascending=False)
        .reset_index()
    )


def collect_control(
    root_dir: os.PathLike,
):
    """
    R1_barcode,R1_primer,G,R1_sgRNA,R1_scaffold_prefix,R1_tail,R1_primer_score,R1_scaffold_prefix_score,R2_barcode,R2_primer,barcode_head,barcode,CTG_target_prefix,R2_sgRNA,pam,target_suffix,C,R2_scaffold_prefix,R2_tail,R2_primer_score,R2_scaffold_prefix_score,R2_sgRNA_score,barcode_score,barcode_id,count
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "control" / "full", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_controls = []
        for parse_file in os.listdir(root_dir / "parse_bar"):
            if (
                utils.infer_cas(parse_file) != "control"
                or utils.infer_chip(parse_file) != chip
            ):
                continue
            df_controls.append(
                pd.read_feather(
                    root_dir / "parse_bar" / parse_file,
                ).drop(columns=["R1_barcode", "R2_barcode"])
            )

        df_controls = pd.concat(df_controls)
        df_controls = agg(df_controls)
        df_controls.to_feather(root_dir / "control" / "full" / f"{chip}.feather")


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
        "R1_primer_score",
        "R1_scaffold_prefix_score",
        "R2_primer_score",
        "R2_sgRNA_score",
        "R2_scaffold_prefix_score",
        "barcode_score",
        "G",
        "C",
        "pam_tail",
        "barcode_id",
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
            "barcode_id",
            "count_full",
            "count_small",
            "R1_primer_length",
            "R2_primer_length",
        ]:
            logy = True
        else:
            logy = False

        if name in ["G", "C", "pam_tail"]:
            df_stat.groupby(name)["count"].sum().plot.bar(
                logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        elif name == "barcode_id":
            df_stat.groupby(name)["count"].sum().reset_index().plot.scatter(
                x=name, y="count", logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")
        else:
            bins = 300 if name.endswith("score") else 150
            df_stat[name].plot.hist(
                bins=bins, weights=df_stat["count"], logy=logy
            ).get_figure().savefig(summary_dir / f"{pathlib.Path(csv_file).stem}.pdf")

        plt.close("all")


def stat_control(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(root_dir / "control" / "full" / f"{chip}.feather")
        stat(df=df_control, save_dir=f"figures/align/stat/control/{chip}")
        draw(
            save_dirs=[f"figures/align/stat/control/{chip}"],
            summary_dir=f"figures/align/stat/control/{chip}",
        )


def filter_control(
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
                G.isin(@G) and \
                C.isin(@C) and \
                pam.str.slice(start=-2).isin(@pam_tail) and \
                count >= @min_count
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

        df_control.to_feather(root_dir / "control" / "filter" / f"{chip}.feather")


def group_stat_control(root_dir: os.PathLike, ext: int, percentage: float):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/stat/group/control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(
            root_dir / "control" / "filter" / f"{chip}.feather"
        )
        df_control = df_control.assign(
            ref1=lambda df: df["R2_primer"]
            + df["barcode_head"]
            + df["barcode"]
            + df["CTG_target_prefix"]
            + df["R2_sgRNA"].str.slice(stop=-3),
            ref2=lambda df: df["R2_sgRNA"].str.slice(start=-3)
            + df["pam"]
            + df["target_suffix"]
            + df["C"]
            + df["R2_scaffold_prefix"]
            + df["R2_tail"],
        )
        assert (
            df_control["ref1"].str.len() >= ext
        ).all(), f"ref1 too short (less than {ext})"
        assert (
            df_control["ref2"].str.len() >= ext
        ).all(), f"ref2 too short (less than {ext})"
        df_control = (
            df_control.assign(
                cut_around=lambda df: df["ref1"].str.slice(start=-ext)
                + df["ref2"].str.slice(stop=ext)
            )
            .groupby(["barcode_id", "cut_around"])["count"]
            .sum()
            .reset_index()
        )
        df_control = (
            df_control.sort_values("count")
            .reset_index(drop=True)
            .assign(
                enrich=lambda df: df.groupby("barcode_id")["count"].transform("cumsum")
                / df.groupby("barcode_id")["count"].transform("sum")
                >= 1 - percentage,
            )
        )

        df_stat = (
            df_control.query("enrich")
            .groupby("barcode_id")
            .size()
            .reset_index()
            .rename(columns={0: "cut_around_type_num"})
        )
        df_stat.to_csv(save_dir / "cut_around_type_num.csv", index=False)

        df_stat.plot.scatter(
            x="barcode_id", y="cut_around_type_num"
        ).get_figure().savefig(save_dir / "cut_around_type_num.pdf")
        plt.close("all")

        df_stat["cut_around_type_num"].clip(upper=100).plot.hist(
            bins=np.linspace(0, 101, 102)
        ).get_figure().savefig(save_dir / "cut_around_type_num_hist.pdf")
        plt.close("all")

        if chip in ["a1", "a2", "a3"]:
            df_plasmid = pd.read_csv(
                "plasmids/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv"
            )
        else:
            df_plasmid = pd.read_csv("plasmids/final_hgsgrna_libb_all_0811-NGG.csv")

        df_plasmid = (
            df_plasmid[["Target sequence"]]
            .rename(columns={"Target sequence": "target"})
            .reset_index(names="barcode_id")
            .assign(
                wt_cut_around=lambda df: df["target"]
                .str.slice(start=22 - ext, stop=22 + ext)
                .map(utils.rev_comp),
            )
        )

        for tem in range(1, 5):
            df_plasmid = df_plasmid.assign(
                **{
                    f"tem{tem}_cut_around": lambda df: (
                        df["target"].str.slice(start=22 - ext, stop=22 + tem)
                        + df["target"].str.slice(start=22 + tem + 1, stop=22 + ext + 1)
                    ).map(utils.rev_comp),
                }
            )

        df_control = df_control.merge(
            right=df_plasmid[
                [
                    "barcode_id",
                    "wt_cut_around",
                    "tem1_cut_around",
                    "tem2_cut_around",
                    "tem3_cut_around",
                    "tem4_cut_around",
                ]
            ],
            how="left",
            on=["barcode_id"],
            validate="many_to_one",
        )

        df_count = df_control.groupby("barcode_id")[["count"]].sum()
        cut_arounds = [
            "wt_cut_around",
            "tem1_cut_around",
            "tem2_cut_around",
            "tem3_cut_around",
            "tem4_cut_around",
        ]
        for i, cut_around in enumerate(cut_arounds):
            column_count = f"count_{cut_around.split('_', 1)[0]}"
            df_count[column_count] = df_control.query(
                f"cut_around == {cut_around}"
            ).set_index("barcode_id")["count"]
            df_count[column_count] = df_count[column_count].fillna(0)
            for j in range(i):
                df_count.loc[
                    df_plasmid.loc[df_count.index, cut_around]
                    == df_plasmid.loc[df_count.index, cut_arounds[j]],
                    column_count,
                ] = 0

        df_count = df_count.assign(
            count_re=lambda df: df["count"]
            - df["count_wt"]
            - df["count_tem1"]
            - df["count_tem2"]
            - df["count_tem3"]
            - df["count_tem4"],
            first=lambda df: df[
                ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
            ].max(axis=1),
            second=lambda df: df[
                ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
            ].apply(lambda row: row.nlargest(2).min(), axis=1),
        ).reset_index()

        df_count.to_csv(save_dir / "count.csv", index=False)

        df_count["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302),
        ).get_figure().savefig(save_dir / "count.pdf")
        plt.close("all")

        (df_count["count_wt"] / df_count["count"]).plot.hist(
            bins=100, logy=True
        ).get_figure().savefig(save_dir / "freq_wt.pdf")
        plt.close("all")

        (df_count["second"] / (df_count["first"] + 1e-6)).plot.hist(
            bins=100, logy=True
        ).get_figure().savefig(save_dir / "freq_second.pdf")
        plt.close("all")


def group_filter_control(
    root_dir: os.PathLike, max_freq_second: float, min_freq_wt: float, min_count: int
):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "ref" / "barcode", exist_ok=True)
    count_dir = pathlib.Path("figures/align/stat/group/control")
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/filter/group/control/{chip}")
        os.makedirs(save_dir, exist_ok=True)
        df_count = pd.read_csv(count_dir / chip / "count.csv", header=0)

        df_stat = pd.DataFrame(columns=["barcode_num"], index=["full", "filter"])
        df_stat.loc["full", "barcode_num"] = df_count.shape[0]

        df_count = df_count.query(
            """
                second / (first + 1e-6) <= @max_freq_second and \
                (count_wt / count >= @min_freq_wt or count >= @min_count)
            """
        ).reset_index(drop=True)

        df_stat.loc["filter", "barcode_num"] = df_count.shape[0]
        df_stat.to_csv(save_dir / "barcode_num.csv")
        df_stat["barcode_num"].plot.bar().get_figure().savefig(
            save_dir / "barcode_num.pdf"
        )

        df_count["dominant"] = df_count[
            ["count_wt", "count_tem1", "count_tem2", "count_tem3", "count_tem4"]
        ].idxmax(axis=1)

        df_dominant = df_count["dominant"].value_counts().reset_index()
        df_dominant.to_csv(save_dir / "dominant.csv", index=False)
        df_dominant.set_index("dominant").plot.bar().get_figure().savefig(
            save_dir / "dominant.pdf"
        )

        if chip in ["a1", "a2", "a3"]:
            df_plasmid = pd.read_csv(
                "plasmids/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv"
            )
        else:
            df_plasmid = pd.read_csv("plasmids/final_hgsgrna_libb_all_0811-NGG.csv")

        df_plasmid = df_plasmid.assign(
            ref=lambda df: (
                "AAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTG"
                + df["Target sequence"]
                + "CAG"
                + df["Barcode2"]
                + df["Primer binding sites 20nt"]
            ).map(utils.rev_comp),
        ).reset_index(names="barcode_id")

        df_ref = df_count[["barcode_id", "dominant", "first"]].merge(
            right=df_plasmid[["barcode_id", "ref"]],
            how="left",
            on=["barcode_id"],
            validate="one_to_one",
        )
        cut = 64
        for tem in range(1, 5):
            df_ref["ref"] = df_ref["ref"].where(
                df_ref["dominant"] != f"count_tem{tem}",
                df_ref["ref"].str.slice(stop=cut - tem - 1)
                + df_ref["ref"].str.slice(start=cut - tem),
            )
        df_ref.to_csv(root_dir / "ref" / "barcode" / f"{chip}.csv", index=False)


def generate_reference(
    root_dir: os.PathLike,
    ext: int = 10,
):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "ref" / "format", exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_ref = pd.read_csv(root_dir / "ref" / "barcode" / f"{chip}.csv", header=0)
        df_ref = df_ref.assign(
            cut=lambda df: df["ref"].str.len() - 79,
            ref1=lambda df: df["ref"]
            .str.slice(stop=64 + ext)
            .where(df["ref"].str.len() == 143, df["ref"].str.slice(stop=63 + ext)),
            ref2=lambda df: df["ref"].str.slice(start=-79 - ext),
            zero=0,
            ext=ext,
            ref2len=lambda df: df["ref2"].str.len(),
        )

        df_ref[["zero", "ref1", "cut", "ext", "ref2", "ref2len"]].to_csv(
            root_dir / "ref" / "format" / f"{chip}.ref",
            sep="\t",
            header=False,
            index=False,
        )


def collect_treat(root_dir: os.PathLike):
    """
    R1_barcode,R1_primer,G,R1_sgRNA,R1_scaffold_prefix,R1_tail,R1_primer_score,R1_scaffold_prefix_score,R2_barcode,R2_primer,barcode_head,barcode,CTG_target_prefix,R2_sgRNA,pam,target_suffix,C,R2_scaffold_prefix,R2_tail,R2_primer_score,R2_scaffold_prefix_score,R2_sgRNA_score,barcode_score,barcode_id,count
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "treat" / "full", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse_bar"):
        if utils.infer_cas(parse_file) == "control":
            continue
        df_treat = pd.read_feather(
            root_dir / "parse_bar" / parse_file,
        ).drop(columns=["R1_barcode", "R2_barcode"])

        df_treat = agg(df_treat)
        df_treat.to_feather(
            root_dir / "treat" / "full" / f"{pathlib.Path(parse_file).stem}.feather"
        )


def stat_treat(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dirs = []
    for treat_file in os.listdir(root_dir / "treat" / "full"):
        df_treat = pd.read_feather(root_dir / "treat" / "full" / treat_file)
        save_dir = f"figures/align/stat/treat/{pathlib.Path(treat_file).stem}"
        save_dirs.append(save_dir)
        stat(
            df=df_treat,
            save_dir=save_dir,
        )

    draw(save_dirs, summary_dir=f"figures/align/stat/treat")


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
                G.isin(@G) and \
                count >= @min_count
            """
        ).reset_index(drop=True)
        df_stat.loc["filter", "row_num"] = df_treat.shape[0]
        df_stat.loc["filter", "count"] = df_treat["count"].sum()

        df_stat.to_csv(save_dir / "stat.csv", index=False)
        df_stat["row_num"].plot.bar().get_figure().savefig(save_dir / "row_num.pdf")
        df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")

        df_treat.to_feather(root_dir / "treat" / "filter" / treat_file)


def demultiplex(
    root_dir: os.PathLike,
):
    """
    |<p5barcode|9-18><p5primer|19>G<sgRNA|20><scaffold|82-92>G<16|CCN|sgRNARC|5>CAG<barcodeRC|18><p7primerRC|21><p7barcodeRC|9-18>|
    """
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "query", exist_ok=True)
    os.makedirs(root_dir / "not_found", exist_ok=True)
    for treat_file in os.listdir(root_dir / "treat" / "filter"):
        chip = utils.infer_chip(treat_file)
        save_dir = pathlib.Path(
            f"figures/align/demultiplex/{chip}/{pathlib.Path(treat_file).stem}"
        )
        os.makedirs(save_dir, exist_ok=True)

        df_ref = pd.read_csv(
            root_dir / "ref" / "barcode" / f"{chip}.csv",
            header=0,
            usecols=["barcode_id"],
        ).reset_index(names="ref_id")

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
                validate="many_to_one",
            )
            .assign(ref_id=lambda df: df["ref_id"].fillna(-1).astype(int))
        )

        df_query.assign(found=lambda df: df["ref_id"] != -1).groupby("found")[
            "count"
        ].sum().reset_index().to_csv(save_dir / "found.csv", index=False)

        df_query["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302), logy=True
        ).get_figure().savefig(save_dir / "count.pdf")
        plt.close("all")

        df_query.query("ref_id == -1")["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302), logy=True
        ).get_figure().savefig(save_dir / "not_found_count.pdf")
        plt.close("all")

        df_query.query("ref_id != -1")["count"].clip(upper=300).plot.hist(
            bins=np.linspace(0, 301, 302), logy=True
        ).get_figure().savefig(save_dir / "found_count.pdf")
        plt.close("all")

        # Get query count distribution for each ref. Queries with multiple refs are distributed to each ref evenly.
        df_query.groupby("ref_id")["count"].sum().reset_index().to_csv(
            save_dir / "query_count_per_ref.csv", index=False
        )

        df_query.query("ref_id == -1")[["query", "count"]].reset_index(
            drop=True
        ).to_csv(
            root_dir / "not_found" / f"{pathlib.Path(treat_file).stem}.query",
            sep="\t",
            header=False,
            index=False,
        )
        df_query.query("ref_id != -1").reset_index(drop=True)[
            ["query", "count", "ref_id"]
        ].to_csv(
            root_dir / "query" / f"{pathlib.Path(treat_file).stem}.query",
            sep="\t",
            header=False,
            index=False,
        )


def draw_demultiplex():
    save_dir = pathlib.Path("figures/align/demultiplex")

    df_founds = []
    df_query_count_per_refs = []
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        for stem in os.listdir(save_dir / chip):
            df_founds.append(
                pd.read_csv(save_dir / chip / stem / "found.csv", header=0).assign(
                    stem=stem
                )
            )
            df_query_count_per_refs.append(
                pd.read_csv(
                    save_dir / chip / stem / "query_count_per_ref.csv", header=0
                ).assign(stem=stem)
            )

    df_found = (
        pd.concat(df_founds)
        .reset_index(drop=True)
        .pivot_table(values="count", index="stem", columns="found", aggfunc="sum")
    )
    df_found.set_index("stem").plot.bar(stacked=True).get_figure().savefig(
        save_dir / "found.pdf"
    )
    pd.concat(df_query_count_per_refs).groupby("ref_id")[
        "count"
    ].sum().reset_index().plot.scatter(x="ref_id", y="count").get_figure().savefig(
        save_dir / "query_count_per_ref.pdf"
    )
