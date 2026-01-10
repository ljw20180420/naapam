#!/usr/bin/env python

import re

import pandas as pd

df = pd.read_csv(
    "~/sdc1/sx_data/SX/algs/wt1-g2n-1.R2.fq.alg.gz",
    sep="\t",
    names=[
        "index",
        "count",
        "score",
        "ref_id",
        "updangle",
        "ref_start1",
        "query_start1",
        "ref_end1",
        "query_end1",
        "random_insertion",
        "ref_start2",
        "query_start2",
        "ref_end2",
        "query_end2",
        "downdangle",
        "cut1",
        "cut2",
        "ref",
        "query",
    ],
    keep_default_na=False,
)

# total = df.groupby("ref_id")["count"].sum()
# wt = (
#     df.query("ref_end1 == 100 and ref_start2 == 154 and random_insertion == ''")
#     .groupby("ref_id")["count"]
#     .sum()
# )
# rs = (wt / total).fillna(0)
# (rs == 1).sum()


def write_reads(df: pd.DataFrame):
    ref_id = (
        df.query("ref_end1 == 95 and ref_start2 == 150 and random_insertion == ''")
        .groupby("ref_id")["count"]
        .sum()
        .sort_values(ascending=False)
        .reset_index()
        .loc[0, "ref_id"]
    )
    df = df.query("ref_id == @ref_id").sort_values("count", ascending=False)

    with open("temp.txt", "w") as fd:
        output_ref = True
        for ref, query, count in zip(df["ref"], df["query"], df["count"]):
            cutoff = re.search("[acgtn]", ref).span()[1]
            ref = ref[cutoff:]
            query = query[cutoff:]
            mid = re.search("[acgtn]", ref).span()[1]
            ref = ref[mid - 47 : mid - 17] + "|" + ref[mid + 17 : mid + 47]
            query = query[mid - 47 : mid - 17] + "|" + query[mid + 17 : mid + 47]
            if output_ref:
                fd.write(f"{ref}\n")
                output_ref = False
            fd.write(f"{query}\t{count}\n")


write_reads(df)
