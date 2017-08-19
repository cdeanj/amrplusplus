#!/usr/bin/env python3

import argparse

amr_level_names = {0: 'Class', 1: 'Mechanism', 2: 'Group'}


def amr_load_data(file_name_list, threshold):
    samples = {}
    labels = set()
    for file in file_name_list:
        with open(file, 'r') as f:
            data = f.read().split('\n')[1:]
            for entry in data:
                if not entry:
                    continue
                entry = entry.split('\t')
                coverage = float(entry[3])
                if coverage < threshold:
                    continue
                sample = entry[0].split('.')[0]
                count = float(entry[2])
                gene_name = entry[1]
                try:
                    samples[sample][gene_name] = count
                except KeyError:
                    try:
                        samples[sample].setdefault(gene_name, count)
                    except KeyError:
                        samples.setdefault(sample, {gene_name: count})
                labels.add(gene_name)
    return samples, labels


def output_amr_analytic_data(outdir, S, L):
    with open(outdir + '/AMR_analytic_matrix.csv', 'w') as amr:
        local_sample_names = []
        for sample, dat in S.items():
            local_sample_names.append(sample)
        amr.write(','.join(local_sample_names) + '\n')
        for label in L:
            local_counts = []
            amr.write(label + ',')
            for local_sample in local_sample_names:
                if label in S[local_sample]:
                    local_counts.append(str(S[local_sample][label]))
                else:
                    local_counts.append(str(0))
            amr.write(','.join(local_counts) + '\n')


parser = argparse.ArgumentParser()
parser.add_argument('--input_file_list', '-i', nargs='+',
                    help='Use globstar to pass a list of files, (Ex: *.tsv)')
parser.add_argument('--output_directory', '-o',
                    help='Output directory for writing the AMR_analytic_matrix.csv file')
parser.add_argument('--threshold', '-t', type=float, default=80,
                    help='Gene coverage threshold under which to discard counts [0, 100]')


if __name__ == '__main__':
    args = parser.parse_args()
    S, L = amr_load_data(args.input_file_list, args.threshold)
    output_amr_analytic_data(args.output_directory, S, L)
