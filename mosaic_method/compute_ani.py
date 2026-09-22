#!/usr/bin/env python3
import sqlite3
import pandas as pd


def compute_ani(res_df, mld_db_path, taxon_df, cluster_level):
    """Adds the average ANI for each comparison."""
    con = sqlite3.connect(mld_db_path)
    cur = con.cursor()
    cur.execute("CREATE TEMP TABLE temp_level (cluster TEXT, genome TEXT);")
    taxon_level = taxon_df[[cluster_level, "genome"]].values
    cur.executemany("INSERT INTO temp_level(cluster, genome) VALUES (?, ?)", taxon_level)
    ani_df = pd.read_sql_query("""select average_divergence, t1.cluster as cluster_1, t2.cluster as cluster_2
    from lastz
    left join temp_level t1 on t1.genome = lastz.genome1
    left join temp_level t2 on t2.genome = lastz.genome2;""", con)
    cur.execute("DROP TABLE temp_level;")
    ani_df["cluster_1"] = ani_df["cluster_1"].astype(str)
    ani_df["cluster_2"] = ani_df["cluster_2"].astype(str)
    # drop inner comparisons
    ani_df = ani_df[ani_df["cluster_1"] != ani_df["cluster_2"]]
    mask = ani_df["cluster_1"] > ani_df["cluster_2"]
    ani_df.loc[mask, ["cluster_1", "cluster_2"]] = ani_df.loc[mask, ["cluster_2", "cluster_1"]].values
    avg_ani_df = ani_df.groupby(["cluster_1", "cluster_2"], dropna=False).mean().reset_index()
    # group by cluster level if different from database
    res_df = pd.merge(res_df, avg_ani_df, "left", left_on=["species_1", "species_2"], right_on=["cluster_1", "cluster_2"]).drop(["cluster_1", "cluster_2"], axis=1)
    return res_df
