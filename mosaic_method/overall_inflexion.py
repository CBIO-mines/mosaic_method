#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd

from mosaic_method.fitting import theoretical_mld


def inflexions(res_df, muc, mus, delta, L0_df, smal_dif, min_r_infl=1):
    """
    Computes if mc > mh for the fitted parameters.
    """
    res_df["r_infl"] = 0
    res_df["infl_exist"] = "no"
    match_lengths = np.logspace(0, 4, 1000)
    L0_df["bac1"] = L0_df["bac1"].astype(str)
    L0_df["bac2"] = L0_df["bac2"].astype(str)
    if not "L0" in res_df:
        res_L0_df = pd.merge(
            res_df, L0_df, left_on=["species_1", "species_2"], right_on=["bac1", "bac2"]
        )
        res_L0_df = res_L0_df.drop(["bac1", "bac2"], axis=1)
    else:
        res_L0_df = res_df

    for i, row in res_L0_df.iterrows():
        mh, mc = theoretical_mld(
            [row["log10tau"], row["log10rho"]],
            smal_dif,
            match_lengths,
            mus,
            muc,
            delta,
            row["L0"],
            False,
        )
        if np.any(mc - mh > 0):
            r_inflexion = match_lengths[np.max(np.where(mc - mh > 0))]
            r_inflexion = int(r_inflexion)
            if r_inflexion > min_r_infl:
                res_L0_df.at[i, "infl_exist"] = "yes"
        else:
            r_inflexion = 0
        res_L0_df.at[i, "r_infl"] = r_inflexion
    return res_L0_df
