import os
import pathlib
import subprocess
import sys
from importlib import resources

import pandas as pd
import pysam
from Bio import Align, Seq

from . import utils


def get_embedding_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"
    aligner.open_left_deletion_score = 0
    aligner.extend_left_deletion_score = 0
    aligner.open_right_deletion_score = 0
    aligner.extend_right_deletion_score = 0

    return aligner


def get_suffix_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"
    aligner.open_left_deletion_score = 0
    aligner.extend_left_deletion_score = 0
    aligner.open_right_deletion_score = 0
    aligner.extend_right_deletion_score = 0
    aligner.open_right_insertion_score = 0
    aligner.extend_right_insertion_score = 0

    return aligner


def fuzzy_split(
    target: str, delimiter: str, aligner: Align.PairwiseAligner
) -> tuple[str | int]:
    if target == "":
        return "", "", "", 0
    if delimiter == "":
        return "", "", target, 0
    segs = target.split(delimiter, maxsplit=1)
    if len(segs) == 2:
        score = aligner.substitution_matrix[0, 0] * len(delimiter)
        return (
            segs[0],
            delimiter,
            segs[1],
            score,
        )
    alignment = aligner.align(target, delimiter)[0]
    if alignment.aligned.shape[1] == 0:
        start = 0
        end = 0
    else:
        start = alignment.aligned[0][0][0]
        end = alignment.aligned[0][-1][-1]
    return (
        alignment.target[:start],
        alignment.target[start:end],
        alignment.target[end:],
        alignment.score,
    )


def parse_read(
    read: str,
    primer_temp: str,
    scaffold_temps: list[str],
    embedding_aligner: Align.PairwiseAligner,
    suffix_aligner: Align.PairwiseAligner,
) -> tuple[str | list[str]]:
    barcode, primer, remain, primer_score = fuzzy_split(
        target=read, delimiter=primer_temp, aligner=embedding_aligner
    )
    scaffold_prefix_score = -float("inf")
    for idx, scaffold_temp in enumerate(scaffold_temps):
        new = fuzzy_split(
            target=remain, delimiter=scaffold_temp, aligner=suffix_aligner
        )
        if new[-1] > scaffold_prefix_score:
            variant, scaffold_prefix, tail, scaffold_prefix_score = new
            used_scaffold_idxs = [idx]
        elif new[-1] == scaffold_prefix_score:
            used_scaffold_idxs.append(idx)

    return (
        barcode,
        primer,
        variant,
        scaffold_prefix,
        tail,
        primer_score,
        scaffold_prefix_score,
        used_scaffold_idxs,
    )


def parse_R1(R1_variant: str) -> tuple[str]:
    G = R1_variant[:1]
    R1_sgRNA = R1_variant[1:]
    return G, R1_sgRNA


def parse_R2(
    R2_variant: str, R1_sgRNA: str, embedding_aliger: Align.PairwiseAligner
) -> tuple[str]:
    barcode_CTG_target = R2_variant[:-1]
    C = R2_variant[-1:]
    barcode_CTG_target_prefix, R2_sgRNA, pam_target_suffix, R2_sgRNA_score = (
        fuzzy_split(
            target=barcode_CTG_target, delimiter=R1_sgRNA, aligner=embedding_aliger
        )
    )
    pam = pam_target_suffix[:3]
    target_suffix = pam_target_suffix[3:]

    return barcode_CTG_target_prefix, R2_sgRNA, pam, target_suffix, C, R2_sgRNA_score


def main():
    unique_file = pathlib.Path(sys.argv[1])
    parse_file = (
        unique_file.parent.parent / "parse" / "nobar" / f"{unique_file.stem}.parse"
    )
    os.makedirs(parse_file.parent, exist_ok=True)
    chip = utils.infer_chip(unique_file)
    embedding_aligner = get_embedding_aligner()
    suffix_aligner = get_suffix_aligner()
    R1_primer_temp = utils.p5primer()
    R2_primer_temp = utils.p7primer()
    R1_scaffold_temps = [utils.scaffold(chip), utils.scaffold_alt(chip)]
    R2_scaffold_temps = [utils.RCscaffold(chip), utils.RCscaffold_alt(chip)]
    with open(unique_file, "r") as rd, open(parse_file, "w") as wd:
        wd.write(
            "R1_barcode\tR1_primer\tG\tR1_sgRNA\tR1_scaffold_prefix\tR1_tail\tR1_primer_score\tR1_scaffold_prefix_score\tR2_barcode\tR2_primer\tbarcode_CTG_target_prefix\tR2_sgRNA\tpam\ttarget_suffix\tC\tR2_scaffold_prefix\tR2_tail\tR2_primer_score\tR2_scaffold_prefix_score\tR2_sgRNA_score\tcount\n"
        )
        for line in rd:
            R1, R2, count = line.split()
            (
                R1_barcode,
                R1_primer,
                R1_variant,
                R1_scaffold_prefix,
                R1_tail,
                R1_primer_score,
                R1_scaffold_prefix_score,
                used_scaffold_idxs,
            ) = parse_read(
                read=R1,
                primer_temp=R1_primer_temp,
                scaffold_temps=R1_scaffold_temps,
                embedding_aligner=embedding_aligner,
                suffix_aligner=suffix_aligner,
            )
            G, R1_sgRNA = parse_R1(R1_variant)
            (
                R2_barcode,
                R2_primer,
                R2_variant,
                R2_scaffold_prefix,
                R2_tail,
                R2_primer_score,
                R2_scaffold_prefix_score,
                _,
            ) = parse_read(
                read=R2,
                primer_temp=R2_primer_temp,
                scaffold_temps=[R2_scaffold_temps[idx] for idx in used_scaffold_idxs],
                embedding_aligner=embedding_aligner,
                suffix_aligner=suffix_aligner,
            )
            (
                barcode_CTG_target_prefix,
                R2_sgRNA,
                pam,
                target_suffix,
                C,
                R2_sgRNA_score,
            ) = parse_R2(R2_variant, R1_sgRNA, embedding_aligner)

            wd.write(
                f"{R1_barcode}\t{R1_primer}\t{G}\t{R1_sgRNA}\t{R1_scaffold_prefix}\t{R1_tail}\t{R1_primer_score}\t{R1_scaffold_prefix_score}\t{R2_barcode}\t{R2_primer}\t{barcode_CTG_target_prefix}\t{R2_sgRNA}\t{pam}\t{target_suffix}\t{C}\t{R2_scaffold_prefix}\t{R2_tail}\t{R2_primer_score}\t{R2_scaffold_prefix_score}\t{R2_sgRNA_score}\t{count}\n"
            )


def build_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "barcode" / "index", exist_ok=True)
    with resources.as_file(
        resources.files(".plasmids")
        / "plasmids/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv"
    ) as pf:
        df_plasmid = pd.read_csv(pf, header=0)
    with open(root_dir / "barcode" / "index" / "barcode.fa", "w") as fd:
        for i, barcode in enumerate(df_plasmid["Barcode2"]):
            barcode = str(Seq.Seq(barcode).reverse_complement())
            fd.write(f">b{i}\n{barcode}\n")

        subprocess.run(
            args=[
                "bowtie2-build",
                (root_dir / "barcode" / "index" / "barcode.fa").as_posix(),
                (root_dir / "barcode" / "index" / "barcode").as_posix(),
            ]
        )


def prepare_barcode_CTG_target_prefix(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "barcode" / "read", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse" / "nobar"):
        df = pd.read_csv(
            root_dir / "parse" / "nobar" / parse_file,
            sep="\t",
            header=0,
            usecols=["barcode_CTG_target_prefix"],
            keep_default_na=False,
        )
        df["barcode_CTG_target_prefix"] = df["barcode_CTG_target_prefix"].where(
            df["barcode_CTG_target_prefix"] != "", "N"
        )
        df.to_csv(
            root_dir / "barcode" / "read" / f"{pathlib.Path(parse_file).stem}.csv",
            header=False,
            index=False,
        )


def map_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "barcode" / "align", exist_ok=True)
    for csv_file in os.listdir(root_dir / "barcode" / "read"):
        subprocess.run(
            args=[
                "bowtie2",
                "--quiet",
                "-p",
                "24",
                "--mm",
                "--norc",
                "--local",
                "-L",
                "7",
                "--ma",
                "1",
                "--mp",
                "2,2",
                "--rdg",
                "3,1",
                "--rfg",
                "3,1",
                "--score-min",
                "C,1",
                "-r",
                "-x",
                (root_dir / "barcode" / "index" / "barcode").as_posix(),
                "-U",
                (root_dir / "barcode" / "read" / csv_file).as_posix(),
                "-S",
                (
                    root_dir
                    / "barcode"
                    / "align"
                    / f"{pathlib.Path(csv_file).stem}.sam"
                ).as_posix(),
            ],
        )


def parse_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "parse" / "bar", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse" / "nobar"):
        sam_file = f"{pathlib.Path(parse_file).stem}.sam"
        sam = pysam.AlignmentFile(root_dir / "barcode" / "align" / sam_file)
        barcode_heads = []
        barcodes = []
        CTG_target_prefixs = []
        barcode_ids = []
        barcode_scores = []
        for align in sam.fetch():
            if (align.flag // 4) % 2 == 1:
                barcode_heads.append("")
                barcodes.append("")
                CTG_target_prefixs.append(align.seq if align.seq != "N" else "")
                barcode_ids.append(-1)
                barcode_scores.append(0)
            else:
                barcode_heads.append(align.seq[: align.qstart])
                barcodes.append(align.seq[align.qstart : align.qend])
                CTG_target_prefixs.append(align.seq[align.qend :])
                barcode_ids.append(align.rname)
                barcode_scores.append(align.get_tag("AS"))

        df_parse = (
            pd.read_csv(
                root_dir / "parse" / "nobar" / parse_file,
                sep="\t",
                header=0,
                keep_default_na=False,
            )
            .drop(columns="barcode_CTG_target_prefix")
            .assign(
                barcode_head=barcode_heads,
                barcode=barcodes,
                CTG_target_prefix=CTG_target_prefixs,
                barcode_id=barcode_ids,
                barcode_score=barcode_scores,
            )[
                [
                    "R1_barcode",
                    "R1_primer",
                    "G",
                    "R1_sgRNA",
                    "R1_scaffold_prefix",
                    "R1_tail",
                    "R1_primer_score",
                    "R1_scaffold_prefix_score",
                    "R2_barcode",
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
                    "R2_primer_score",
                    "R2_scaffold_prefix_score",
                    "R2_sgRNA_score",
                    "barcode_score",
                    "barcode_id",
                    "count",
                ]
            ]
        )

        df_parse.to_csv(
            root_dir / "parse" / "bar" / parse_file,
            sep="\t",
            index=False,
        )


def build_sgRNA(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "sgRNA" / "index", exist_ok=True)
    with resources.as_file(
        resources.files(".plasmids")
        / "plasmids/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv"
    ) as pf:
        df_plasmid = pd.read_csv(pf, header=0)
    with open(root_dir / "sgRNA" / "index" / "sgRNA.fa", "w") as fd:
        for i, sgRNA in enumerate(df_plasmid["sgRNA"]):
            fd.write(f">b{i}\n{sgRNA}\n")

        subprocess.run(
            args=[
                "bowtie2-build",
                (root_dir / "sgRNA" / "index" / "sgRNA.fa").as_posix(),
                (root_dir / "sgRNA" / "index" / "sgRNA").as_posix(),
            ]
        )


def prepare_R1_sgRNA_and_R2_sgRNA(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for column in ["R1_sgRNA", "R2_sgRNA"]:
        os.makedirs(root_dir / "sgRNA" / "read" / column, exist_ok=True)
        for parse_file in os.listdir(root_dir / "parse" / "bar"):
            df = pd.read_csv(
                root_dir / "parse" / "bar" / parse_file,
                sep="\t",
                header=0,
                usecols=[column],
                keep_default_na=False,
            )
            df[column] = df[column].where(df[column] != "", "N")
            df.to_csv(
                root_dir
                / "sgRNA"
                / "read"
                / column
                / f"{pathlib.Path(parse_file).stem}.csv",
                header=False,
                index=False,
            )


def map_sgRNA(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    for column in ["R1_sgRNA", "R2_sgRNA"]:
        os.makedirs(root_dir / "sgRNA" / "align" / column, exist_ok=True)
        for csv_file in os.listdir(root_dir / "sgRNA" / "read" / column):
            subprocess.run(
                args=[
                    "bowtie2",
                    "--quiet",
                    "-p",
                    "24",
                    "--mm",
                    "--norc",
                    "--local",
                    "-L",
                    "7",
                    "--ma",
                    "1",
                    "--mp",
                    "2,2",
                    "--rdg",
                    "3,1",
                    "--rfg",
                    "3,1",
                    "--score-min",
                    "C,1",
                    "-r",
                    "-x",
                    (root_dir / "sgRNA" / "index" / "sgRNA").as_posix(),
                    "-U",
                    (root_dir / "sgRNA" / "read" / column / csv_file).as_posix(),
                    "-S",
                    (
                        root_dir
                        / "sgRNA"
                        / "align"
                        / column
                        / f"{pathlib.Path(csv_file).stem}.sam"
                    ).as_posix(),
                ],
            )


def parse_sgRNA(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "parse" / "sgRNA", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse" / "bar"):
        sam_file = f"{pathlib.Path(parse_file).stem}.sam"

        R1_sam = pysam.AlignmentFile(
            root_dir / "sgRNA" / "align" / "R1_sgRNA" / sam_file
        )
        R1_sgRNA_ids = []
        R1_sgRNA_bowtie2_scores = []
        for align in R1_sam.fetch():
            if (align.flag // 4) % 2 == 1:
                R1_sgRNA_ids.append(-1)
                R1_sgRNA_bowtie2_scores.append(0)
            else:
                R1_sgRNA_ids.append(align.rname)
                R1_sgRNA_bowtie2_scores.append(align.get_tag("AS"))

        R2_sam = pysam.AlignmentFile(
            root_dir / "sgRNA" / "align" / "R2_sgRNA" / sam_file
        )
        R2_sgRNA_ids = []
        R2_sgRNA_bowtie2_scores = []
        for align in R2_sam.fetch():
            if (align.flag // 4) % 2 == 1:
                R2_sgRNA_ids.append(-1)
                R2_sgRNA_bowtie2_scores.append(0)
            else:
                R2_sgRNA_ids.append(align.rname)
                R2_sgRNA_bowtie2_scores.append(align.get_tag("AS"))

        df_parse = pd.read_csv(
            root_dir / "parse" / "bar" / parse_file,
            sep="\t",
            header=0,
            keep_default_na=False,
        ).assign(
            R1_sgRNA_id=R1_sgRNA_ids,
            R1_sgRNA_bowtie2_score=R1_sgRNA_bowtie2_scores,
            R2_sgRNA_id=R2_sgRNA_ids,
            R2_sgRNA_bowtie2_score=R2_sgRNA_bowtie2_scores,
        )[
            [
                "R1_barcode",
                "R1_primer",
                "G",
                "R1_sgRNA",
                "R1_scaffold_prefix",
                "R1_tail",
                "R1_primer_score",
                "R1_scaffold_prefix_score",
                "R2_barcode",
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
                "R2_primer_score",
                "R2_scaffold_prefix_score",
                "R2_sgRNA_score",
                "barcode_score",
                "R1_sgRNA_bowtie2_score",
                "R2_sgRNA_bowtie2_score",
                "barcode_id",
                "R1_sgRNA_id",
                "R2_sgRNA_id",
                "count",
            ]
        ]

        df_parse.to_csv(
            root_dir / "parse" / "sgRNA" / parse_file,
            sep="\t",
            index=False,
        )
