import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import special

from . import utils


def correct_alg(root_dir: os.PathLike, temperature: float):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "align" / "correct", exist_ok=True)
    for alg_file in os.listdir(root_dir / "align" / "raw"):
        df_alg = utils.read_alg(
            root_dir / "align" / "raw" / alg_file,
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
            ],
        )
        df_query = pd.read_csv(
            root_dir / "query" / "found" / f"{pathlib.Path(alg_file).stem}.query",
            sep="\t",
            names=["query", "count", "ref_id"],
            keep_default_na=False,
        )
        df_alg["index"] = df_query.groupby("query").transform("ngroup")

        chip = utils.infer_chip(alg_file)
        df_ref = pd.read_feather(
            root_dir / "control" / "hq_mut" / f"{chip}.feather"
        ).reset_index(names="ref_id")
        df_alg = df_alg.merge(
            right=df_ref[["ref_id", "count"]].rename(columns={"count": "count_ref"}),
            how="left",
            on=["ref_id"],
            validate="many_to_one",
        ).assign(
            weight=lambda df: (
                df["count_ref"]
                * df.groupby("index")["score"].transform(
                    lambda score, temperature=temperature: special.softmax(
                        score / temperature
                    )
                )
            ),
            count_distri=lambda df: df["count"]
            * df.groupby("index")["weight"].transform(lambda se: se / se.sum()),
        )

        df_alg = df_alg.astype(str)
        df_output = pd.DataFrame(
            {
                "meta": df_alg["index"].str.cat(
                    df_alg[
                        [
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
                            "count_ref",
                            "count_distri",
                        ]
                    ],
                    sep="\t",
                ),
                "ref": df_alg["ref"],
                "query": df_alg["query"],
            },
        )
        df_output.stack().to_csv(
            root_dir / "align" / "correct" / alg_file,
            header=False,
            index=False,
        )


def stat_read(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = root_dir / "figures" / "analyze" / "stat_read"
    os.makedirs(save_dir, exist_ok=True)
    df_algs = []
    for alg_file in os.listdir(root_dir / "align" / "correct"):
        df_algs.append(
            utils.read_alg(
                root_dir / "align" / "correct" / alg_file,
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
                    "count_ref",
                    "count_distri",
                ],
            )
            .groupby("score")["count_distri"]
            .sum()
            .rename("count")
            .reset_index()
        )

    df_alg = pd.concat(df_algs)
    df_alg["score"].plot.hist(bins=300, weights=df_alg["count"]).get_figure().savefig(
        save_dir / "score.pdf"
    )
    plt.close("all")


def collect_data(
    root_dir: os.PathLike,
    min_score: int,
):
    """
    Only collect data. Do not apply any annotation. Only apply read-wise filter such as score.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "full", exist_ok=True)
    os.makedirs(root_dir / "analyze" / "control" / "full", exist_ok=True)
    os.makedirs(root_dir / "figures" / "analyze" / "collect_data", exist_ok=True)

    df_stat = pd.DataFrame(
        {
            "row_num": [0, 0],
            "count": [0, 0],
        },
        index=["full", "filter"],
    )
    df_algs = []
    for alg_file in os.listdir(root_dir / "align" / "correct"):
        df_alg = utils.read_alg(
            root_dir / "align" / "correct" / alg_file,
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
                "count_ref",
                "count_distri",
            ],
        )

        df_stat.loc["full", "row_num"] += df_alg.shape[0]
        df_stat.loc["full", "count"] += df_alg["count_distri"].sum()

        df_alg = df_alg.query("score >= @min_score")

        df_stat.loc["filter", "row_num"] += df_alg.shape[0]
        df_stat.loc["filter", "count"] += df_alg["count_distri"].sum()

        df_alg = (
            df_alg.groupby(
                [
                    "ref_id",
                    "cut1",
                    "cut2",
                    "ref_end1",
                    "random_insertion",
                    "ref_start2",
                ]
            )["count_distri"]
            .sum()
            .rename("count")
            .reset_index()
            .astype(
                {
                    "ref_end1": "int8",
                    "ref_start2": "int8",
                    "cut1": "int8",
                    "cut2": "int8",
                }
            )
            .assign(stem=pathlib.Path(alg_file).stem)
        )
        df_algs.append(df_alg)

    df_stat.to_csv(root_dir / "figures" / "analyze" / "collect_data" / "stat.csv")
    df_stat["row_num"].plot.bar().get_figure().savefig(
        root_dir / "figures" / "analyze" / "collect_data" / "row_num.pdf"
    )
    plt.close("all")
    df_stat["count"].plot.bar().get_figure().savefig(
        root_dir / "figures" / "analyze" / "collect_data" / "count.pdf"
    )
    plt.close("all")

    df_alg = pd.concat(df_algs).assign(cas=lambda df: df["stem"].map(utils.infer_cas))
    df_alg.query("cas != 'control'").drop(columns="cas").reset_index(
        drop=True
    ).to_feather(root_dir / "analyze" / "treat" / "full" / "treat.feather")
    df_alg.query("cas == 'control'").drop(columns="cas").reset_index(
        drop=True
    ).to_feather(root_dir / "analyze" / "control" / "full" / "control.feather")


def stat_ref(root_dir: os.PathLike, min_count_tot: int, max_up_del_size: int):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = root_dir / "figures" / "analyze" / "stat_ref"
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "full" / "treat.feather"
    )
    df_treat.groupby(["stem", "ref_id"])["count"].sum().clip(upper=300).plot.hist(
        bins=np.linspace(0, 301, 302)
    ).get_figure().savefig(save_dir / "count_tot.pdf")
    plt.close("all")

    df_treat = df_treat.assign(
        count_tot=lambda df: df.groupby(["stem", "ref_id"])["count"].transform("sum")
    )

    df_treat.assign(freq_nowt=utils.freq_nowt(df_treat)).query(
        "count_tot >= @min_count_tot"
    ).groupby(["stem", "ref_id"])["freq_nowt"].first().plot.hist(
        bins=100
    ).get_figure().savefig(
        save_dir / "freq_nowt.pdf"
    )
    plt.close("all")

    for tem in range(1, 5):
        df_treat_tem = df_treat.assign(
            count_blunt=lambda df: utils.count_tem_dele(df, tem, 0),
            freq_blunt=lambda df: df["count_blunt"] / (df["count_tot"] + 1e-6),
        )
        for dele in range(1, max_up_del_size + 1):
            df_treat_dele = (
                df_treat_tem.assign(
                    count_dele=lambda df: utils.count_tem_dele(df, tem, dele),
                    freq_dele=lambda df: df["count_dele"] / (df["count_tot"] + 1e-6),
                    freq_dele_rel_blunt=lambda df: df["count_dele"]
                    / (df["count_blunt"] + 1e-6),
                )
                .query("count_tot >= @min_count_tot")
                .groupby(["stem", "ref_id"])
                .agg(
                    freq_blunt=pd.NamedAgg(column="freq_blunt", aggfunc="first"),
                    freq_dele=pd.NamedAgg(column="freq_dele", aggfunc="first"),
                    freq_dele_rel_blunt=pd.NamedAgg(
                        column="freq_dele_rel_blunt", aggfunc="first"
                    ),
                )
            )

            for column in ["freq_blunt", "freq_dele", "freq_dele_rel_blunt"]:
                if column == "freq_blunt":
                    name = f"freq_blunt_tem{tem}"
                elif column == "freq_dele":
                    name = f"freq_dele{dele}_tem{tem}"
                elif column == "freq_dele_rel_blunt":
                    name = f"freq_dele{dele}_tem{tem}_rel_blunt"
                df_treat_dele[column].clip(upper=10).plot.hist(
                    bins=300, logy=True
                ).get_figure().savefig(save_dir / f"{name}.pdf")
                plt.close("all")


def fit_ref(root_dir: os.PathLike, max_up_del_size: int):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "fit", exist_ok=True)
    save_dir = root_dir / "figures" / "analyze" / "fit_ref"
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "full" / "treat.feather"
    )
    df_treat = (
        df_treat.assign(
            up_del_size=lambda df: utils.up_del_size(df),
        )
        .query("ref_end1 <= cut1 and up_del_size <= @max_up_del_size")
        .groupby(["up_del_size"])["count"]
        .sum()
        .reset_index()
        .assign(
            freq_rel_blunt=lambda df: df["count"]
            / df.loc[df["up_del_size"] == 0, "count"].item()
        )
    )

    log_freq_rel_blunt, const_up_del_size = patsy.dmatrices(
        "np.log(freq_rel_blunt) ~ up_del_size", data=df_treat, return_type="dataframe"
    )
    result = sm.OLS(
        endog=log_freq_rel_blunt,
        exog=const_up_del_size,
    ).fit()
    prediction = result.get_prediction(const_up_del_size)
    summary_frame = prediction.summary_frame(alpha=0.05)
    summary_frame = summary_frame.assign(
        up_del_size=const_up_del_size["up_del_size"],
    )
    summary_frame.to_csv(root_dir / "fit" / "summary_frame.csv", index=False)

    plt.scatter(x=df_treat["up_del_size"], y=np.log(df_treat["freq_rel_blunt"]))
    plt.plot(
        summary_frame["up_del_size"],
        summary_frame["mean"],
        color="red",
        label="fitted freq_rel_blunt",
    )
    plt.fill_between(
        summary_frame["up_del_size"],
        summary_frame["mean_ci_lower"],
        summary_frame["mean_ci_upper"],
        color="red",
        alpha=0.3,
        label="95% CI",
    )
    plt.fill_between(
        summary_frame["up_del_size"],
        summary_frame["obs_ci_lower"],
        summary_frame["obs_ci_upper"],
        color="blue",
        alpha=0.2,
        label="95% PI",
    )
    plt.savefig(save_dir / "ci.pdf")
    plt.close("all")


def filter_ref(
    root_dir: os.PathLike,
    min_count_tot: int,
    max_freq_nowt: float,
):
    """
    Use positive mask because nan compare always return False.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "filter" / "ref", exist_ok=True)
    save_dir = root_dir / "figures" / "analyze" / "filter_ref"
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "full" / "treat.feather"
    )

    df_stat = pd.DataFrame(index=["full", "filter"], columns=["ref_num", "count"])
    df_stat.loc["full", "ref_num"] = (
        df_treat[["stem", "ref_id"]].drop_duplicates().shape[0]
    )
    df_stat.loc["full", "count"] = df_treat["count"].sum()

    count_tot = df_treat.groupby(["stem", "ref_id"])["count"].transform("sum")
    mask = count_tot >= min_count_tot
    mask = mask & (utils.freq_nowt(df_treat) <= max_freq_nowt)

    df_treat = df_treat.loc[mask].reset_index(drop=True)

    df_stat.loc["filter", "ref_num"] = (
        df_treat[["stem", "ref_id"]].drop_duplicates().shape[0]
    )
    df_stat.loc["filter", "count"] = df_treat["count"].sum()
    df_stat.to_csv(save_dir / "stat.csv")
    df_stat["ref_num"].plot.bar().get_figure().savefig(save_dir / "ref_num.pdf")
    plt.close("all")
    df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")
    plt.close("all")

    df_treat.to_feather(
        root_dir / "analyze" / "treat" / "filter" / "ref" / "treat.feather"
    )


def stat_mutant(root_dir: os.PathLike, min_count_tot: int):
    root_dir = pathlib.Path(os.fspath(root_dir))
    save_dir = root_dir / "figures" / "analyze" / "stat_mutant"
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "filter" / "ref" / "treat.feather"
    )

    utils.up_del_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "up_del_size.pdf")
    plt.close("all")

    utils.down_del_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "down_del_size.pdf")
    plt.close("all")

    utils.rand_ins_size(df_treat).clip(upper=30).plot.hist(
        bins=np.linspace(0, 31, 32), weights=df_treat["count"]
    ).get_figure().savefig(save_dir / "rand_ins_size.pdf")
    plt.close("all")

    df_treat.assign(
        count_tot=lambda df: df.groupby(["stem", "ref_id"])["count"].transform("sum"),
        freq_mutant=lambda df: utils.freq_mutant(df),
    ).query("count_tot >= @min_count_tot")["freq_mutant"].plot.hist(
        bins=100
    ).get_figure().savefig(
        save_dir / "freq_mutant.pdf"
    )
    plt.close("all")


def filter_mutant(
    root_dir: os.PathLike,
    max_up_del_size: int,
    max_down_del_size: int,
    max_rand_ins_size: int,
    max_freq_mutant: float,
    fit_tem_max: int,
):
    """
    Do not filter mutant because missing mutant are treated as count 0. Use a column legal to mark mutant passing the filter.
    """
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "filter" / "mutant", exist_ok=True)
    save_dir = root_dir / "figures" / "analyze" / "filter_mutant"
    os.makedirs(save_dir, exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "filter" / "ref" / "treat.feather"
    )

    df_stat = pd.DataFrame(index=["full", "filter"], columns=["mutant_num", "count"])
    df_stat.loc["full", "mutant_num"] = df_treat.shape[0]
    df_stat.loc["full", "count"] = df_treat["count"].sum()

    outlier_ratio = np.exp(
        pd.read_csv(root_dir / "fit" / "summary_frame.csv", header=0)["obs_ci_upper"]
    )
    df_treat["tem_indicator"] = df_treat["cut2"] - df_treat["ref_start2"] > fit_tem_max
    for tem in range(1, fit_tem_max + 1):
        df_treat["tem_indicator"] = df_treat["tem_indicator"] | (
            (df_treat["cut2"] - df_treat["ref_start2"] == tem)
            & (df_treat["cut1"] - df_treat["ref_end1"] == 0)
        )
        count_blunt = utils.count_tem_dele(df_treat, tem, 0)
        for up_del_size in range(1, max_up_del_size + 1):
            df_treat["tem_indicator"] = df_treat["tem_indicator"] | (
                (df_treat["cut2"] - df_treat["ref_start2"] == tem)
                & (df_treat["cut1"] - df_treat["ref_end1"] == up_del_size)
                & (
                    utils.count_tem_dele(df_treat, tem, up_del_size)
                    <= outlier_ratio[up_del_size] * count_blunt
                )
            )

    df_treat = df_treat.assign(
        legal=(
            (utils.up_del_size(df_treat) <= max_up_del_size)
            & (utils.down_del_size(df_treat) <= max_down_del_size)
            & (utils.rand_ins_size(df_treat) <= max_rand_ins_size)
            & ((utils.freq_mutant(df_treat) <= max_freq_mutant) | utils.is_wt(df_treat))
        )
    )

    df_stat.loc["filter", "mutant_num"] = df_treat["legal"].sum()
    df_stat.loc["filter", "count"] = df_treat.query("legal")["count"].sum()
    df_stat.to_csv(save_dir / "stat.csv")
    df_stat["mutant_num"].plot.bar().get_figure().savefig(save_dir / "mutant_num.pdf")
    plt.close("all")
    df_stat["count"].plot.bar().get_figure().savefig(save_dir / "count.pdf")
    plt.close("all")

    df_treat.to_feather(
        root_dir / "analyze" / "treat" / "filter" / "mutant" / "treat.feather"
    )


def kim_correct(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "correct", exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "merge" / "treat.feather"
    )
    count_kim, count_tot_kim, freq_mut_kim = utils.kim(df_treat)
    df_treat = df_treat.assign(
        count_kim=count_kim,
        count_tot_kim=count_tot_kim,
        freq_mut_kim=freq_mut_kim,
    )

    df_treat.to_feather(
        root_dir / "analyze" / "treat" / "correct" / "treat.feather",
    )


def annote_columns(root_dir: os.PathLike):
    root_dir = pathlib.Path(os.fspath(root_dir))
    os.makedirs(root_dir / "analyze" / "treat" / "annote", exist_ok=True)
    os.makedirs(root_dir / "analyze" / "control" / "annote", exist_ok=True)

    df_treat = pd.read_feather(
        root_dir / "analyze" / "treat" / "correct" / "treat.feather"
    )
    df_control = pd.read_feather(
        root_dir / "analyze" / "control" / "dup" / "control.feather"
    )

    # Read barcode and sgRNA from the generated fasta file to prevent the inconsistency when use custom plasmid file.
    df_bar = (
        pd.read_csv(root_dir / "barcode" / "index" / "barcode.fa", names=["barcode"])
        .loc[1::2, :]
        .reset_index(drop=True)
        .reset_index(names="barcode_id")
    )
    df_sgRNA = (
        pd.read_csv(root_dir / "sgRNA" / "index" / "sgRNA.fa", names=["sgRNA"])
        .loc[1::2, :]
        .reset_index(drop=True)
        .reset_index(names="barcode_id")
    )
    df_ref = []
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_ref.append(
            pd.read_csv(
                root_dir / "ref" / f"{chip}.ref",
                sep="\t",
                names=["zero", "ref1", "cut1", "ext", "ref2", "ref2len"],
            )[["ref1", "ref2"]]
            .reset_index(names="ref_id")
            .assign(chip=chip)
        )
    df_ref = pd.concat(df_ref).reset_index(drop=True)

    df_column = []
    for chip in ["a1", "a2", "a3", "g1n", "g2n", "g3n"]:
        df_column.append(
            pd.read_feather(root_dir / "control" / "hq_mut" / f"{chip}.feather")[
                ["barcode_id"]
            ]
            .reset_index(names="ref_id")
            .assign(chip=chip)
        )
    df_column = pd.concat(df_column).reset_index(drop=True)

    df_column = (
        df_column.merge(
            right=df_bar,
            how="left",
            on=["barcode_id"],
            validate="many_to_one",
        )
        .merge(
            right=df_sgRNA,
            how="left",
            on=["barcode_id"],
            validate="many_to_one",
        )
        .merge(
            right=df_ref,
            how="left",
            on=["chip", "ref_id"],
            validate="many_to_one",
        )
    )

    df_treat = df_treat.merge(
        right=df_column,
        how="left",
        on=["chip", "ref_id"],
        validate="many_to_one",
    ).assign(mh=lambda df: df.apply(utils.get_mh, axis=1))

    df_treat.to_csv(
        root_dir / "analyze" / "treat" / "annote" / "treat.csv",
        index=False,
        na_rep="NA",
    )

    df_control = df_control.merge(
        right=df_column,
        how="left",
        on=["chip", "ref_id"],
        validate="many_to_one",
    ).assign(mh=lambda df: df.apply(utils.get_mh, axis=1))

    df_control.to_csv(
        root_dir / "analyze" / "control" / "annote" / "control.csv",
        index=False,
        na_rep="NA",
    )
