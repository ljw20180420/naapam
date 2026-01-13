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


def collect_plasmid(plasmid_dir: os.PathLike) -> pd.DataFrame:
    plasmid_dir = pathlib.Path(os.fspath(plasmid_dir))
    df_plasmids = []
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_plasmids.append(
            pd.read_csv(f"plasmids/{chip}.csv", names=["xuhao", "plasmid"]).assign(
                chip=chip
            )
        )

    df_plasmid = pd.concat(df_plasmids).reset_index(drop=True)
    df_plasmid.to_feather(plasmid_dir / "plasmid.feather")


def collect_control(
    unique_dir: os.PathLike,
    control_dir: os.PathLike,
):
    unique_dir = pathlib.Path(os.fspath(unique_dir))
    control_dir = pathlib.Path(os.fspath(control_dir))
    os.makedirs(control_dir, exist_ok=True)
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_controls = []
        for unique_file in os.listdir(unique_dir):
            if (
                utils.infer_cas(unique_file) != "control"
                or utils.infer_chip(unique_file) != chip
            ):
                continue
            df_controls.append(
                pd.read_csv(
                    unique_dir / unique_file, sep="\t", names=["R1", "R2", "count"]
                )
            )

        df_control = (
            pd.concat(df_controls)
            .groupby(["R1", "R2"])["count"]
            .sum()
            .sort_values(ascending=False)
            .reset_index()
        )

        df_control.to_feather(control_dir / f"{chip}.feather")


def stat_control_step1(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/control/{chip}/step1")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(control_dir / f"{chip}.feather")
        df_control.query("count <= 100 and count >=2")["count"].plot.hist(
            bins=99
        ).get_figure().savefig(save_dir / "count.pdf")
        plt.close("all")
        df_control.query(
            "not R1.str.contains(@utils.p5primer()) and count <= 100 and count >=2"
        )["count"].plot.hist(bins=99).get_figure().savefig(save_dir / "no_u6_count.pdf")
        plt.close("all")
        df_control.query(
            "not R2.str.contains(@utils.p7primer()) and count <= 100 and count >=2"
        )["count"].plot.hist(bins=99).get_figure().savefig(save_dir / "no_p7_count.pdf")
        plt.close("all")


def filter_control_step1(
    control_dir: os.PathLike, min_count: int, min_no_u6_count: int, min_no_p7_count: int
):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = (
            pd.read_feather(control_dir / f"{chip}.feather")
            .query(
                """
                    count >= @min_count and \
                    (R1.str.contains(@utils.p5primer()) or count >= @min_no_u6_count) and \
                    (R2.str.contains(@utils.p7primer()) or count >= @min_no_p7_count)
                """
            )
            .reset_index(drop=True)
        )
        df_control.to_feather(control_dir / f"{chip}_step1.feather")


def strip_control_p5p7barcodes(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))

    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(control_dir / f"{chip}_step1.feather")
        # strip p5barcode
        R1sp = df_control["R1"].str.split(pat=utils.p5primer(), n=1, regex=False)
        R1sp = pd.DataFrame(
            R1sp.map(
                lambda segs, aligner=utils.get_embedding_aligner(), delimiters=[
                    utils.p5primer()
                ]: utils.fuzzy_split(segs, aligner, delimiters)
            ).tolist(),
            columns=["p5barcode", "p5primer", "R1", "p5primer_score"],
        )
        df_control = df_control.assign(
            p5primer=R1sp["p5primer"],
            R1=R1sp["R1"],
            p5primer_score=R1sp["p5primer_score"],
        )

        # strip p7barcode
        R2sp = df_control["R2"].str.split(pat=utils.p7primer(), n=1, regex=False)
        R2sp = pd.DataFrame(
            R2sp.map(
                lambda segs, aligner=utils.get_embedding_aligner(), delimiters=[
                    utils.p7primer()
                ]: utils.fuzzy_split(segs, aligner, delimiters)
            ).tolist(),
            columns=["p7barcode", "p7primer", "R2", "p7primer_score"],
        )
        df_control = df_control.assign(
            p7primer=R2sp["p7primer"],
            R2=R2sp["R2"],
            p7primer_score=R2sp["p7primer_score"],
        )

        # group p5primer, R1, p7primer, R2
        df_control = (
            df_control.groupby(["p5primer", "R1", "p7primer", "R2"])
            .agg(
                p5primer_score=pd.NamedAgg(column="p5primer_score", aggfunc="first"),
                p7primer_score=pd.NamedAgg(column="p7primer_score", aggfunc="first"),
                count=pd.NamedAgg(column="count", aggfunc="sum"),
            )
            .reset_index()
        )
        df_control.to_feather(control_dir / f"{chip}_strip_p5p7barcode.feather")


def stat_control_step2(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/control/{chip}/step2")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(control_dir / f"{chip}_strip_p5p7barcode.feather")
        df_control["p5primer_score"].plot.hist(
            bins=100, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "p5primer_score.pdf")
        plt.close("all")
        df_control["p7primer_score"].plot.hist(
            bins=100, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "p7primer_score.pdf")
        plt.close("all")
        df_control["R1"].str.len().plot.hist(
            bins=150, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R1_no_primer_length.pdf")
        plt.close("all")
        df_control["R2"].str.len().plot.hist(
            bins=150, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R2_no_primer_length.pdf")
        plt.close("all")


def filter_control_step2(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(control_dir / f"{chip}_strip_p5p7barcode.feather")
        df_control = df_control.query(
            """
                R1.str.len() >= 100 and \
                R2.str.len() >= 100 and \
                p5primer_score >= 33 and \
                p7primer_score >= 37
            """
        ).reset_index()
        df_control.to_feather(control_dir / f"{chip}_step2.feather")


def extract_sgRNA_barcode_target(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_control = pd.read_feather(control_dir / f"{chip}_step2.feather")
        R1sp = df_control["R1"].str.split(utils.scaffold(chip), n=1, regex=False)
        R1sp = pd.DataFrame(
            R1sp.map(
                lambda segs, aligner=utils.get_suffix_aligner(), delimiters=[
                    utils.scaffold(chip),
                    utils.scaffold_alt(chip),
                ]: utils.fuzzy_split(segs, aligner, delimiters)
            ).tolist(),
            columns=["sgRNA", "R1_scaffold", "R1", "R1_scaffold_score"],
        )
        df_control = df_control.assign(
            G_start=R1sp["sgRNA"].str.slice(stop=1),
            sgRNA=R1sp["sgRNA"].str.slice(start=1),
            R1_scaffold=R1sp["R1_scaffold"],
            R1=R1sp["R1"],
            R1_scaffold_score=R1sp["R1_scaffold_score"],
        )

        R2sp = df_control["R2"].str.split(utils.RCscaffold(chip), n=1, regex=False)
        R2sp = pd.DataFrame(
            R2sp.map(
                lambda segs, aligner=utils.get_suffix_aligner(), delimiters=[
                    utils.RCscaffold(chip),
                    utils.RCscaffold_alt(chip),
                ]: utils.fuzzy_split(segs, aligner, delimiters)
            ).tolist(),
            columns=[
                "barcode_CTG_target",
                "R2_RCscaffold",
                "R2",
                "R2_RCscaffold_score",
            ],
        )
        df_control = df_control.assign(
            barcode_CTG_target=R2sp["barcode_CTG_target"].str.slice(stop=-1),
            C_RCstart=R2sp["barcode_CTG_target"].str.slice(start=-1),
            R2_RCscaffold=R2sp["R2_RCscaffold"],
            R2=R2sp["R2"],
            R2_RCscaffold_score=R2sp["R2_RCscaffold_score"],
        )

        df_control.to_feather(
            control_dir / f"{chip}_extract_sgRNA_barcode_target.feather"
        )


def stat_control_step3(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        save_dir = pathlib.Path(f"figures/align/control/{chip}/step3")
        os.makedirs(save_dir, exist_ok=True)
        df_control = pd.read_feather(
            control_dir / f"{chip}_extract_sgRNA_barcode_target.feather"
        )
        df_control.groupby("G_start")["count"].sum().plot.bar().get_figure().savefig(
            save_dir / "G_start_count.pdf"
        )
        plt.close("all")
        df_control["sgRNA"].str.len().plot.hist(
            bins=30, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "sgRNA_length.pdf")
        plt.close("all")
        df_control["R1_scaffold"].str.len().plot.hist(
            bins=100, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R1_scaffold_length.pdf")
        plt.close("all")
        df_control["R1"].str.len().plot.hist(
            bins=30, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R1_length.pdf")
        plt.close("all")
        df_control["R1_scaffold_score"].plot.hist(
            bins=200, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R1_scaffold_score.pdf")
        plt.close("all")
        df_control["barcode_CTG_target"].str.len().plot.hist(
            bins=70, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "barcode_CTG_target_length.pdf")
        plt.close("all")
        df_control.groupby("C_RCstart")["count"].sum().plot.bar().get_figure().savefig(
            save_dir / "C_RCstart_count.pdf"
        )
        plt.close("all")
        df_control["R2_RCscaffold"].str.len().plot.hist(
            bins=100, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R2_RCscaffold_length.pdf")
        plt.close("all")
        df_control["R2"].str.len().plot.hist(
            bins=30, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R2_length.pdf")
        plt.close("all")
        df_control["R2_RCscaffold_score"].plot.hist(
            bins=200, weights=df_control["count"]
        ).get_figure().savefig(save_dir / "R2_RCscaffold_score.pdf")
        plt.close("all")


def filter_control_step3(control_dir: os.PathLike):
    control_dir = pathlib.Path(os.fspath(control_dir))
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        # TODO filter not G following U6
        # TODO filter sgRNA length
        # TODO filter R1_scaffold_score
        # TODO filter not G following scaffold
        # TODO filter R2_RCscaffold_score
        pass


def find_sgRNA_in_target():
    pass


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
