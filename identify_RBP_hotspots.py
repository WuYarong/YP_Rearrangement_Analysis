#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
import argparse
import sys



def poisson_test(observed_count, expected_count):
    p_value = 1 - poisson.cdf(observed_count, mu=expected_count)
    return p_value

def replace_small_p_value(p_value, threshold=1e-16):
    return max(threshold, p_value)


def main(input_file, significance_level, apply_bonferroni, prefix):
    significant_results = []
    corrected_p_values = []

    output_file_all = f"{prefix}_all_results.xls"
    output_file_significant = f"{prefix}_significant_results.xls"


    output_fa = open(output_file_all, 'w')
    output_fs = open(output_file_significant, 'w')

    tmp_bnd_info_dict = dict()

    with open(input_file, 'r') as f:
        total_sv_count = 0
        for line in f.readlines():
            line = line.replace('\n', '')
            if line.startswith('#'):
                header = line.split('\t')
                output_fa.write('{}\tuncorrected_p_value\tcorrected_p_value\n'.format(line))
                output_fs.write('{}\tuncorrected_p_value\tcorrected_p_value\n'.format(line))
            else:
                tmp = line.split('\t')
                bnd_id = tmp[1]
                total_sv_count += int(tmp[-1])
                tmp_bnd_info_dict[bnd_id] = tmp

    total_bnd_count = len(tmp_bnd_info_dict.keys())
    expected_count = total_sv_count / float(total_bnd_count)
    print('total_sv_bnd_count:{}\ttotal_bnd_count:{}\texpected_count:{}'.format(total_sv_count, total_bnd_count, expected_count))

    for key in tmp_bnd_info_dict.keys():
        observed_count = int(tmp_bnd_info_dict[key][-1])            
        p_value = poisson_test(observed_count, expected_count)
        p_value = replace_small_p_value(p_value)  # Replace small p-values
        significant_results.append([key, observed_count, p_value])


    if apply_bonferroni:
        # Collect all p-values for Bonferroni correction
        all_p_values = [p_value for _, _, p_value in significant_results]
        # _, corrected_p_values, _, _ = multipletests(all_p_values, alpha=significance_level, method='bonferroni')
        _, corrected_p_values, _, _ = multipletests(all_p_values, alpha=significance_level, method='fdr_bh')
    else:
        corrected_p_values = [p_value for _, _, p_value in significant_results]


    # if apply_bonferroni:
    #     # Collect all p-values for Bonferroni correction
    #     all_p_values = [p_value for _, _, p_value in significant_results]
    #     num_tests = len(all_p_values)
    #     corrected_p_values = [min(1.0, p * num_tests) for p in all_p_values]
    # else:
    #     corrected_p_values = [p_value for _, _, p_value in significant_results]


    for i in range(0, len(significant_results)):
        bnd_id = significant_results[i][0]
        p_value = float(significant_results[i][2])
        corrected_p_value = corrected_p_values[i]
        output_fa.write('{}\t{}\t{}\n'.format('\t'.join(tmp_bnd_info_dict[bnd_id]), p_value, corrected_p_value))
        
        if corrected_p_value < significance_level:
             output_fs.write('{}\t{}\t{}\n'.format('\t'.join(tmp_bnd_info_dict[bnd_id]), p_value, corrected_p_value))


def print_usage():
    print("Usage: python script.py input_file [--significance_level SIGNIFICANCE_LEVEL] [--apply_fdr] [--prefix PREFIX]")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Identification of RBP (rearrangement-related breakpoint) hotspots")
    parser.add_argument("input_file", help="tab-delimited file")
    parser.add_argument("--significance_level", type=float, default=0.05, help="significance level")
    parser.add_argument("--apply_fdr", action="store_true", help="whether to perform DFR correction")
    parser.add_argument("--prefix", default="output", help="prefix of output fiels")

    args = parser.parse_args()
    main(args.input_file, args.significance_level, args.apply_bonferroni, args.prefix)
