import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.genmod.families.links as links


def fit_model(df: pd.DataFrame, cut1: int, cut2: int) -> pd.DataFrame:
    def exp_fit(df: pd.DataFrame, cut1: int, cut2: int):
        zero_del_df = df.query("ref_end1 == @cut1")[
            ["rep", "chip", "sample", "barcode", "ref_start2", "true_frequency"]
        ].reset_index(drop=True)
        df = df.merge(
            right=zero_del_df,
            how="inner",
            on=["rep", "chip", "sample", "barcode", "ref_start2"],
            suffixes=["", "_zero_del"],
            validate="many_to_one",
        ).assign(
            delete_size=lambda df: cut1 - df["ref_end1"],
            normalized_true_frequency=lambda df: df["true_frequency"]
            / df["true_frequency_zero_del"],
        )
        df.plot.scatter(
            x="delete_size",
            y="normalized_true_frequency",
            # s="weight",
        ).get_figure().savefig("shit.png")
        plt.close("all")
        results = smf.ols(
            "np.log(normalized_true_frequency) ~ delete_size - 1",
            data=df,
        ).fit()
        plt.scatter(df["delete_size"], np.log(df["normalized_true_frequency"]))
        plt.plot(df["delete_size"], results.predict(df["delete_size"]))
        plt.savefig("shit2.png")

    fitted_exp = (
        df.query("ref_start2 < @cut2 and ref_end1 <= @cut1")
        .query(
            "ref_end1 - @cut1 != ref_start2 - @cut2 and ref_end1 - @cut1 != ref_start2 - @cut2 - 1"
        )
        .groupby(
            [
                "cas",
                "polq",
                "rep",
                "chip",
                "sample",
                "barcode",
                "ref_end1",
                "ref_start2",
            ]
        )
        .agg(
            weight=pd.NamedAgg(column="count", aggfunc="sum"),
            true_frequency=pd.NamedAgg(column="true_frequency", aggfunc="sum"),
        )
        .reset_index()
        .groupby(["cas", "polq"])
        .apply(
            lambda df, cut1=cut1, cut2=cut2: exp_fit(df, cut1, cut2),
            include_groups=True,
        )
    )
