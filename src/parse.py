from Bio import Align


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
