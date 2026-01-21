import os
import pathlib
import subprocess

import pandas as pd
import pysam
from Bio import Align, Seq


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


def read(
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


def R1(R1_variant: str) -> tuple[str]:
    G = R1_variant[:1]
    R1_sgRNA = R1_variant[1:]
    return G, R1_sgRNA


def R2(
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


def build_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "barcode")
    df_plasmid = pd.read_csv(
        "plasmids/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt.csv",
        header=0,
    )
    with open(root_dir / "barcode" / "barcode.fa", "w") as fd:
        for i, barcode in enumerate(df_plasmid["Barcode2"]):
            barcode = str(Seq.Seq(barcode).reverse_complement())
            fd.write(f">b{i}\n{barcode}\n")

        subprocess.run(
            args=[
                "bowtie2-build",
                (root_dir / "barcode" / "barcode.fa").as_posix(),
                (root_dir / "barcode" / "barcode").as_posix(),
            ]
        )


def prepare_barcode_CTG_target_prefix(root_dir: os.PathLike):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "bowtie2", exist_ok=True)
    for parse_file in os.listdir(root_dir / "parse"):
        df = pd.read_csv(
            root_dir / "parse" / parse_file,
            sep="\t",
            header=0,
            usecols=["barcode_CTG_target_prefix"],
            keep_default_na=False,
        )
        df["barcode_CTG_target_prefix"] = df.where(
            df["barcode_CTG_target_prefix"] != "", "N"
        )
        df.to_csv(
            root_dir / "bowtie2" / f"{pathlib.Path(parse_file).stem}.csv",
            header=False,
            index=False,
        )


def map_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(root_dir)
    for csv_file in os.listdir(root_dir / "bowtie2"):
        if not csv_file.endswith(".csv"):
            continue
        subprocess.run(
            args=[
                "bowtie2",
                "--quiet",
                "--mm",
                "--norc",
                "--local",
                "-L",
                "15",
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
                (root_dir / "barcode" / "barcode").as_posix(),
                "-U",
                (root_dir / "bowtie2" / csv_file).as_posix(),
                "-S",
                (
                    root_dir / "bowtie2" / f"{pathlib.Path(csv_file).stem}.sam"
                ).as_posix(),
            ],
        )


def parse_barcode(root_dir: os.PathLike):
    root_dir = pathlib.Path(root_dir)
    os.makedirs(root_dir / "parse_bar", exist_ok=True)
    for sam_file in os.listdir(root_dir / "bowtie2"):
        if not sam_file.endswith(".sam"):
            continue

        sam = pysam.AlignmentFile(root_dir / "bowtie2" / sam_file)
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
                root_dir / "parse" / f"{pathlib.Path(sam_file).stem}.parse",
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

        df_parse.to_feather(
            root_dir / "parse_bar" / f"{pathlib.Path(sam_file).stem}.parse"
        )
