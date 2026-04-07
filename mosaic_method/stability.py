#!/usr/bin/env python3


from fit import Lllocal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def around_Llocal(distance, opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit = False):
    """
    Explores the vicinity of optimal parameters fitted with fit.fit_params.
    """
    theta = np.linspace(0, 2*np.pi, 100)
    tau_values = distance*np.cos(theta) + opt_pars[0]
    rho_values = distance*np.sin(theta) + opt_pars[1]
    opt_Llocal = Lllocal(opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit = False)
    Llvalues = []
    for i in range(len(tau_values)):
        Llvalues.append(Lllocal((tau_values[i], rho_values[i]), empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit = False))
    Llvalues = np.array(Llvalues)
    fig, ax = plt.subplots()
    ax.plot(theta, Llvalues)
    ax.axhline(y=opt_Llocal, c="r")
    ax.axvline(x = 0, c="olivedrab")
    ax.axvline(x = np.pi, c="yellowgreen")
    ax.axvline(x = 2*np.pi, c="olivedrab")
    ax.axvline(x = np.pi/2, c="navy")
    ax.axvline(x = 3*np.pi/2, c="lightblue")
    ax.text(0.2, np.max(Llvalues), f"radius (but in log scale, so not a circle) : {distance}")
    return fig
    # fig.savefig("../Bac_Esch_aroundfit.png")
