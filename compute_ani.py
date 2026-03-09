#!/usr/bin/env python3
import sqlite3
import pandas as pd

def switch_cols(row):
    if row["cluster_1"] > row["cluster_2"]:
        tmp = row["cluster_1"]
        row["cluster_1"] = row["cluster_2"]
        row["cluster_2"] = tmp
    return row

def compute_ani(res_df, mld_db_path):
    """Adds the average ANI for each comparison."""
    con = sqlite3.connect(mld_db_path)
    ani_df = pd.read_sql_query("""select average_divergence, t1.cluster as cluster_1, t2.cluster as cluster_2
    from lastz
    left join taxon t1 on t1.genome = lastz.genome1
    left join taxon t2 on t2.genome = lastz.genome2;""", con)
    ani_df = ani_df.apply(lambda x: switch_cols(x), axis=1)
    avg_ani_df = ani_df.groupby(["cluster_1", "cluster_2"], dropna=False).mean().reset_index()
    res_df = pd.merge(res_df, avg_ani_df, "left", left_on=["species_1", "species_2"], right_on=["cluster_1", "cluster_2"]).drop(["cluster_1", "cluster_2"], axis=1)
    return res_df
