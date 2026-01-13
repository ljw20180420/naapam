import os
import pathlib
import re

from Bio import Align, Seq


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


def infer_time(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"-(?:a1|a2|a3|g1n|g2n|g3n)-(\d)", name).group(1)


def infer_sample(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"^(.+?)-", name).group(1)


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
    segs: list[str], aligner: Align.PairwiseAligner, delimiters: list[str]
) -> list[str]:
    if len(segs) == 2:
        match_score = aligner.substitution_matrix[0, 0]
        prefix, suffix = segs
        primer_score = match_score * len(delimiters[0])
        return [prefix, delimiters[0], suffix, primer_score]

    max_score = -float("inf")
    for delimiter in delimiters:
        alignment_new = aligner.align(segs[0], delimiter)[0]
        if alignment_new.score > max_score:
            alignment = alignment_new
            max_score = alignment_new.score

    if alignment.aligned.shape[1] == 0:
        start = 0
        end = 0
    else:
        start = alignment.aligned[0][0][0]
        end = alignment.aligned[0][-1][-1]
    return [
        alignment.target[:start],
        alignment.target[start:end],
        alignment.target[end:],
        alignment.score,
    ]
