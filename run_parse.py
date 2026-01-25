#!/usr/bin/env python

import os
import pathlib
import sys

from src import parse, utils

unique_file = pathlib.Path(sys.argv[1])
parse_file = unique_file.parent.parent / "parse" / "nobar" / f"{unique_file.stem}.parse"
os.makedirs(parse_file.parent, exist_ok=True)
chip = utils.infer_chip(unique_file)
embedding_aligner = parse.get_embedding_aligner()
suffix_aligner = parse.get_suffix_aligner()
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
        ) = parse.read(
            read=R1,
            primer_temp=R1_primer_temp,
            scaffold_temps=R1_scaffold_temps,
            embedding_aligner=embedding_aligner,
            suffix_aligner=suffix_aligner,
        )
        G, R1_sgRNA = parse.R1(R1_variant)
        (
            R2_barcode,
            R2_primer,
            R2_variant,
            R2_scaffold_prefix,
            R2_tail,
            R2_primer_score,
            R2_scaffold_prefix_score,
            _,
        ) = parse.read(
            read=R2,
            primer_temp=R2_primer_temp,
            scaffold_temps=[R2_scaffold_temps[idx] for idx in used_scaffold_idxs],
            embedding_aligner=embedding_aligner,
            suffix_aligner=suffix_aligner,
        )
        barcode_CTG_target_prefix, R2_sgRNA, pam, target_suffix, C, R2_sgRNA_score = (
            parse.R2(R2_variant, R1_sgRNA, embedding_aligner)
        )

        wd.write(
            f"{R1_barcode}\t{R1_primer}\t{G}\t{R1_sgRNA}\t{R1_scaffold_prefix}\t{R1_tail}\t{R1_primer_score}\t{R1_scaffold_prefix_score}\t{R2_barcode}\t{R2_primer}\t{barcode_CTG_target_prefix}\t{R2_sgRNA}\t{pam}\t{target_suffix}\t{C}\t{R2_scaffold_prefix}\t{R2_tail}\t{R2_primer_score}\t{R2_scaffold_prefix_score}\t{R2_sgRNA_score}\t{count}\n"
        )
