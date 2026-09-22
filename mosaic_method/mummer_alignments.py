#!/usr/bin/env python3
#
import re
import subprocess as sp
import tempfile
import os

import pandas as pd
import numpy as np


def conc_show_aligns_matches(mummer_out):
    """
    Concatenates all matches and mismatches lines from a mummer show-aligns alignment.
    """

    # concatenate all lines
    misandmatches = ""

    # is the previous line an alignment line (with atgc)
    prev_align = False
    cur_align = False

    # calc the padding before the alignment
    counted_start_letter = False
    padding = 0


    lines = mummer_out.splitlines()
    for line in lines:
        if line == "":
            continue
        prev_align = cur_align
        if line[0].isdigit():
            cur_align = True
        else:
            cur_align = False
        if not counted_start_letter and cur_align:
            counted_start_letter = True
            detect_bases = re.compile(r"[atgcATGC]")
            base_match = detect_bases.search(line)
            if base_match:
                padding = base_match.start()

        if prev_align and not cur_align:
            pattern = re.compile(r"[^ \^]")
            if pattern.search(line):
                print("Error in line : ")
                print(line)
                break
            misandmatches += line[padding:]
    return misandmatches


def mummer_count_matches(misandmatches):
    """
    Counts the length of matches in a string of matches and mismatches.
    """
    matches = []
    cur_match = 0
    for char in misandmatches:
        if char == "^":
            if cur_match:
                matches.append(cur_match)
                cur_match = 0
        elif char == " ":
            cur_match += 1
        else:
            print(f"Error with char {char}")
    if cur_match:
        matches.append(cur_match)
    count_array = np.zeros(max(matches), dtype=np.int64)
    for match_length in matches:
        # first element is matches of length 1
        count_array[match_length - 1] += 1
    return count_array


def calc_count_array_len(summed_count_array):
    total_match = 0
    for i, count in enumerate(summed_count_array):
        total_match += (i+1)*count
    return total_match


def parse_show_align_output(mummer_out):
    """
    Splits a mummer show-aligns output in multiple alignments.
    """
    count_arrays = []
    length1s = []
    sum_matches = []
    alignment = ""
    split_align_pattern = r'(?s)(BEGIN.*?\])(.*?)(END.*?)(?=\n|$)'
    for m in re.finditer(split_align_pattern, mummer_out):
        alignment = m.group(2).strip("\n")
        misandmatches = conc_show_aligns_matches(alignment)
        count_array = mummer_count_matches(misandmatches)
        sum_matches.append(calc_count_array_len(count_array))
        count_arrays.append(count_array)
        stats = m.group(1).split(" ")
        length1 = abs(int(stats[6]) - int(stats[4]))
        length1s.append(length1)
    return count_arrays, length1s, sum_matches


def run_mummer(target, query, prefix):
    """
    Runs mummer, adds matches length counts together and returns them as a counter.
    """
    # flag to know if a temporary file was created
    temp_create = False
    sequence_names = re.compile(r">[^\n]*")
    with open(target, "r") as f:
        target_names = sequence_names.findall(f.read())
    with open(query, "r") as f:
        query_names = sequence_names.findall(f.read())
    if len(target_names) > 1:
        temp_create = True
        # concatenate all target sequences in one file
        with open(target, "r") as f:
            target_seq = f.read()
        # remove all sequence names except the first one
        target_seq_lines = target_seq.splitlines()
        target_seq_lines = [target_seq_lines[0]] + [line for line in target_seq_lines[1:] if not line.startswith(">")]
        target_seq = "\n".join(target_seq_lines)
        # tempfile with target sequences concatenated
        target_f = tempfile.NamedTemporaryFile(mode="w", delete=False)
        target_f.write(target_seq)
        target_f.close()
        target = target_f.name

    mummer_status = "ok"
    try:
        cmd_out_mummer = sp.run(
        ['nucmer', '--mum', "--prefix", prefix, target, query],
        check=True,
        capture_output=True
        )
    except sp.CalledProcessError as e:
        mummer_output = e.stdout
        mummer_status = "error"
        mummer_stderr = e.stderr

    summed_count_array = None
    summed_matches = 0
    summed_aligned = 0
    average_divergence = 1
    res = {}
    if mummer_status != "error":
        # check for empty delta file (no alignment found)
        with open(f"{prefix}.delta", "r") as f:
            if len(f.readlines()) <= 2:
                if temp_create:
                    os.remove(target)
                os.remove(f"{prefix}.delta")
                res_tuple = (np.zeros(1, dtype=np.int64), 1, 0, 0)
                res["result"] = res_tuple
                res["status"] = mummer_status
                return res
        target_names = [tn[1:] for tn in target_names]
        query_names = [qn[1:] for qn in query_names]

        for query_name in query_names:
            shal_command = ['show-aligns', '-r', f"{prefix}.delta", target_names[0].split(" ")[0], query_name.split(" ")[0]]
            try:
                res_show_aligns = sp.run(
                    shal_command,
                    capture_output=True,
                    check=True,
                    encoding="utf-8"
                )
            except sp.CalledProcessError:
                # case where a given contig does not have a single match in the target
                # another possibility would be to scan for them beforehand
                # but better ask for forgiveness
                continue
            count_arrays, length1s, sum_matches = parse_show_align_output(res_show_aligns.stdout)
            for count_array, length1, sum_match in zip(count_arrays, length1s, sum_matches):
                if summed_count_array is None:
                    summed_count_array = count_array
                    summed_matches = sum_match
                    summed_aligned = length1
                    continue
                if len(count_array) < len(summed_count_array):
                    count_array = np.pad(count_array, (0, len(summed_count_array) - len(count_array)), mode='constant')
                elif len(count_array) > len(summed_count_array):
                    summed_count_array = np.pad(summed_count_array, (0, len(count_array) - len(summed_count_array)), mode='constant')
                summed_count_array += count_array
                summed_matches += sum_match
                summed_aligned += length1
        if summed_aligned == 0:
            average_divergence = 1
            summed_count_array = np.zeros(1, dtype=np.int64)
        else:
            average_divergence = 1 - summed_matches / summed_aligned
        if temp_create:
            os.remove(target)
        os.remove(f"{prefix}.delta")
    else:
        res["stderr"] = mummer_stderr
        res["stdout"] = mummer_output
    if summed_count_array is None:
        summed_count_array = np.zeros(1, dtype=np.int64)
    res_tuple = (summed_count_array, average_divergence, summed_matches, summed_aligned)
    res["result"] = res_tuple
    res["status"] = mummer_status
    return res
