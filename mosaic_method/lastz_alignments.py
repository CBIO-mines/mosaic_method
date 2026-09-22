#!/usr/bin/env python3

import io
import re
import subprocess as sp

import pandas as pd
import numpy as np


def parse_cigarx_line(line):
    tokens = re.findall(r"\d*=", line)
    if not tokens:
        return np.zeros(1, dtype=np.int64)
    matches = [int(t[:-1]) if len(t) > 1 else 1 for t in tokens]
    count_array = np.zeros(max(matches), dtype=np.int64)
    for match_length in matches:
        count_array[match_length - 1] += 1
    return count_array


def run_lastz(target, query):
    """
    Runs lastz, adds cigarx matches length counts together and returns them as a counter.
    """
    target = target + "[multiple]"
    query = query  # +"[multiple]"
    lastz_cmd = [
        "lastz",
        target,
        query,
        "--format=general:length1,idfrac,cigarx",
        "--allocate:traceback=2000M",
    ]
    try:
        res_lastz = sp.run(lastz_cmd, capture_output=True, check=True, encoding="utf-8")
        lastz_status = "ok"
    except sp.CalledProcessError as e:
        lastz_output = e.stdout
        lastz_status = "error"
        lastz_stderr = e.stderr
        raise RuntimeError("LASTZ failed")

    summed_count_array = np.zeros(1, np.int64)
    summed_matches = 0
    summed_aligned = 0
    average_divergence = 0
    res = {}
    if lastz_status != "error":
        res["status"] = lastz_status
        res_lastz_df = pd.read_csv(io.StringIO(res_lastz.stdout), sep="\t")
        for _, row in res_lastz_df.iterrows():
            # cigarx counting
            count_array = parse_cigarx_line(row["cigarx"])
            # Pad the arrays with zeroes if they have different sizes
            if len(count_array) < len(summed_count_array):
                count_array = np.pad(
                    count_array,
                    (0, len(summed_count_array) - len(count_array)),
                    mode="constant",
                )
            elif len(count_array) > len(summed_count_array):
                summed_count_array = np.pad(
                    summed_count_array,
                    (0, len(count_array) - len(summed_count_array)),
                    mode="constant",
                )
            summed_count_array += count_array

            # idfrac calculations
            summed_matches += int(row["idfrac"].split("/")[0])
            summed_aligned += int(row["idfrac"].split("/")[1])
        if summed_aligned == 0:
            average_divergence = 1
        else:
            average_divergence = 1 - summed_matches / summed_aligned
    else:
        res["status"] = "error"
        res["stderr"] = lastz_stderr
        res["stdout"] = lastz_output
    res_tuple = (summed_count_array, average_divergence, summed_matches, summed_aligned)
    res["result"] = res_tuple
    return res
