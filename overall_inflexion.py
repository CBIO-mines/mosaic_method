#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd

from fit.fit import theoretical_mld

def inflexions(res_df, muc, mus, delta, L0_df, smal_dif, min_r_infl = 1):
    """
    Computes if mc > mh for the fitted parameters.
    """
    res_df["r_infl"] = 0
    res_df["infl_exist"] = "no"
    match_lengths = np.logspace(0, 4, 1000)
    for i, row in res_df.iterrows():
        current_L0 = L0_df.loc[
            (L0_df["bac1"] == row["species_1"]) & (L0_df["bac2"] == row["species_2"])
            , "L0"
        ].values
        mh, mc = theoretical_mld(
            [row["log10tau"], row["log10rho"]],
            smal_dif,
            match_lengths,
            mus,
            muc,
            delta,
            current_L0,
            False
        )
        if np.any(mc - mh > 0):
            r_inflexion = match_lengths[np.max(np.where(mc - mh > 0))]
            r_inflexion = int(r_inflexion)
            if r_inflexion > min_r_infl:
                res_df.at[i, "infl_exist"] = "yes"
        else:
            r_inflexion = 0
        res_df.at[i, "r_infl"] = r_inflexion
    return res_df
