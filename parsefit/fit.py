#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import dual_annealing


def sum_mlds(mld_comp_df):
    """
    Takes a df of mlds by comparison and sums it.
    """
    summed_colmuns = mld_comp_df.drop(labels="comp", axis=1).sum(axis=0, skipna=True)
    summed_df = summed_colmuns.reset_index()
    summed_df.columns = ["match_length", "freq"]
    summed_df = summed_df.astype({"match_length" : "int64"})
    return summed_df.sort_values(by=["match_length"]).reset_index(drop=True)


def bin_mld(summed_df, linear_bin_width, limit_size, power_increment, ncomp):
    """
    Bins and normalizes a summed mld according to a specific pattern.

    Parameters
    ----------
    summed_df: pd.DataFrame
    an unbinned dataframe with a column "match_length" and a column "freq"
    linear_bin_width: float
    the width of the bin in the linear part of bin vector
    limit_size: float
    the limit at which the bin vector switches from linear to log
    power_increment: float
    the "bin width" of the log part
    ncomp: float
    the number of comparison summed here

    Returns
    -------
    a binned pd.Dataframe with columns :
    - "match_length" containing the geometric mean of the bin
    - "freq" containing the counts corresponding to the bin
    """
    match_bin = list(np.arange(0.5, limit_size, linear_bin_width))
    initial_len = len(match_bin)
    cur_power = 0.1
    while match_bin[-1] < max(summed_df["match_length"]):
        match_bin += [match_bin[initial_len - 1]*10**(cur_power)]
        cur_power += power_increment

    res = pd.DataFrame.from_dict(
        {"match_length" : match_bin,
         "freq" : [0] * len(match_bin)}
    )
    summed_row = 0
    binned_row = 0
    while(summed_row < summed_df.shape[0]):
        if(summed_df.loc[summed_row, "match_length"] <= res.loc[binned_row, "match_length"]):
            res.loc[binned_row, "freq"] += summed_df.loc[summed_row, "freq"]
            summed_row += 1
        else:
            binned_row += 1
            if binned_row > res.shape[0] - 1:
                break

    for binned_row in range(res.shape[0] - 1):
        len_bin = res.loc[binned_row + 1, "match_length"] - res.loc[binned_row, "match_length"]
        res.loc[binned_row + 1, "freq"] /= len_bin * ncomp

    gmean_match_length = np.sqrt(np.array(match_bin)[1:]*np.array(match_bin)[:-1])
    res.drop([0], inplace=True)
    res["match_length"] = gmean_match_length
    res.reset_index(drop=True)

    return res




def theoretical_mld(opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit = False):
    """
    Computes the theoretical match length distribution according to the paper.
    """
    ml_low = match_lengths - smal_dif
    ml_hi = match_lengths + smal_dif
    tau = np.double(10.0**opt_pars[0])
    rho = np.double(10.0**opt_pars[1])
    if L0_fit:
        L0 = np.double(10.0**opt_pars[2])
    mua = min(delta/tau, mus)

    mc = 2*((1 + match_lengths*mua*tau)/np.exp(match_lengths*mua*tau) - (1 + match_lengths*muc*tau)/np.exp(match_lengths*muc*tau))/(match_lengths**2*(muc**2 - mus**2)*tau**2)
    mc_low = 2*((1 + ml_low*mua*tau)/np.exp(ml_low*mua*tau) - (1 + ml_low*muc*tau)/np.exp(ml_low*muc*tau))/(ml_low**2*(muc**2 - mus**2)*tau**2)
    mc_hi = 2*((1 + ml_hi*mua*tau)/np.exp(ml_hi*mua*tau) - (1 + ml_hi*muc*tau)/np.exp(ml_hi*muc*tau))/(ml_hi**2*(muc**2 - mus**2)*tau**2)

    mc = L0*(mc_low + mc_hi - 2*mc)/smal_dif**2
    np.nan_to_num(mc, copy=False)

    if tau < delta/mus:
        mh = (2*(-np.exp(-(match_lengths*muc*tau)) + np.exp(-(match_lengths*mus*tau)) + match_lengths*(-muc + mus)*tau))/(match_lengths**2*(-muc**2 + mus**2)*tau)
        mh_low = (2*(-np.exp(-(ml_low*muc*tau)) + np.exp(-(ml_low*mus*tau)) + ml_low*(-muc + mus)*tau))/(ml_low**2*(-muc**2 + mus**2)*tau)
        mh_hi = (2*(-np.exp(-(ml_hi*muc*tau)) + np.exp(-(ml_hi*mus*tau)) + ml_hi*(-muc + mus)*tau))/(ml_hi**2*(-muc**2 + mus**2)*tau)
    else:
        mh = (-2*(-np.exp(-(match_lengths*muc*tau)) + match_lengths*(-muc + mus)*tau + (1 + match_lengths*(delta - mus*tau))/np.exp(match_lengths*delta)))/(match_lengths**2*(muc**2 - mus**2)*tau)
        mh_low = (-2*(-np.exp(-(ml_low*muc*tau)) + ml_low*(-muc + mus)*tau + (1 + ml_low*(delta - mus*tau))/np.exp(ml_low*delta)))/(ml_low**2*(muc**2 - mus**2)*tau)
        mh_hi = (-2*(-np.exp(-(ml_hi*muc*tau)) + ml_hi*(-muc + mus)*tau + (1 + ml_hi*(delta - mus*tau))/np.exp(ml_hi*delta)))/(ml_hi**2*(muc**2 - mus**2)*tau)

    mh = L0*rho*(mh_low + mh_hi - 2*mh)/smal_dif**2

    return mh, mc


def Lllocal(opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit = False):
    """
    The squared relative difference to minimize.
    """
    mh_calc, mc_calc = theoretical_mld(opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit)
    mt_calc = mh_calc + mc_calc
    return np.mean(((mt_calc - empirical_mld)/(mt_calc + empirical_mld))**2.)


def fit_params(opt_method, init_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0):
    """
    Interface to minimize from scipy
    """
    if L0:
        L0_fit = False
    else:
        L0_fit = True

    if opt_method in ["nelder-mead", "Nelder-Mead", "BFGS", "L-BFGS-B"]:
        res_opt = minimize(
            Lllocal,
            init_pars,
            method=opt_method,
            args=(
                empirical_mld,
                smal_dif,
                match_lengths,
                mus,
                muc,
                delta,
                L0,
                L0_fit
            ),
            options={'xatol': 1e-8, 'disp': True}
        )
    elif opt_method == "dual-annealing":
        res_opt = dual_annealing(
            Lllocal,
            bounds = [(4, 10), (-15, -4)],
            args=(
                empirical_mld,
                smal_dif,
                match_lengths,
                mus,
                muc,
                delta,
                L0,
                L0_fit
            )
        )
    else:
        sys.exit("Unexistent/unimplemented optimization method requested")

    return res_opt


def write_results(binned_mld, opted_pars, out_mld, out_pars, L0):
    """
    Writes the results of the fit and parsing to specified files.
    """
    binned_mld.to_csv(out_mld, index=False)
    with open(out_pars, "w") as outfile:
        outfile.write("log10tau,log10rho,L0\n")
        outfile.write(f"{opted_pars[0]}")
        for par in opted_pars[1:]:
            outfile.write(f",{par}")
        if L0:
            outfile.write(f",{L0}")
        outfile.write("\n")



def plot_surface(min_logtau, max_logtau, min_logrho, max_logrho, num_points, output_file, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, fitted_params=None):
    """Plots the Lllocal surface in a given region of the parameters to optimize."""

    x_range = np.linspace(min_logtau, max_logtau, num_points)
    y_range = np.linspace(min_logrho, max_logrho, num_points)
    x_vals, y_vals = np.meshgrid(x_range, y_range)

    # Calculate the corresponding Z values using Lllocal
    z_vals = np.zeros_like(x_vals)
    for i in range(len(x_range)):
        for j in range(len(y_range)):
            opt_pars = (x_vals[i, j], y_vals[i, j])
            z_vals[i, j] = Lllocal(opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0)

    # Create a 3D surface plot
    fig = plt.figure(figsize=(10, 8))

    for i, azim in enumerate(range(0, 280, 90)):
        ax = fig.add_subplot(221 + i, projection="3d", computed_zorder=False)
        if fitted_params is not None:
            ax.scatter(fitted_params[0], fitted_params[1],
                       Lllocal(fitted_params, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0),
                       color='red', s=100, label='Fitted Parameters', zorder=10)
            ax.legend()
        ax.plot_surface(x_vals, y_vals, z_vals, cmap="viridis", zorder=1)
        ax.view_init(azim=azim, elev=60)
        ax.set_title(f"View with Azimuth = {azim}")
        # Set labels for the axes
        ax.set_xlabel("logtau")
        ax.set_ylabel("logrho")
        ax.set_zlabel("Lllocal")

    fig.tight_layout()

    # Save the plot to the specified output file
    plt.savefig(output_file, dpi=300)
