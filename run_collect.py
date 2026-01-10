#!/usr/bin/env python

import os
import pathlib

from src.collect import collect_data

table_dir = pathlib.Path("/home/ljw/sdc1/alltablefiles_table")
treat_tables = [
    table_dir / table
    for table in os.listdir(table_dir)
    if table.endswith(".table") and not table.startswith("wt")
]
df_treat = collect_data(treat_tables, cut1=50, cut2=60, min_score=-1)
df_treat.to_feather("/home/ljw/sdc1/naapam.treat.feather")


table_dir = pathlib.Path("/home/ljw/sdc1/alltablefiles_table")
control_tables = [
    table_dir / table
    for table in os.listdir(table_dir)
    if table.endswith(".table") and table.startswith("wt")
]
df_control = (
    collect_data(control_tables, cut1=50, cut2=60, min_score=-1)
    .drop(columns=["cas", "polq"])
    .rename(columns={"sample": "rep"})
)
df_control.to_feather("/home/ljw/sdc1/naapam.control.feather")
