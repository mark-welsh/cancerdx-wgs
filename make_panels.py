import argparse
import pandas as pd
import numpy as np
import sys
from multiprocessing import Pool


def get_snp_info(df, snp_pos):
    new_df = df[(df["start"] < snp_pos) & (df["stop"] > snp_pos)]
    try:
        idx = new_df.index.values[0]
    except IndexError:
        return np.nan, np.nan, np.nan

    gene = new_df["gene"].tolist()[0]
    category = new_df["category"].tolist()[0]

    return idx, gene, category


def snp_map(dfs):
    chrm_top_data = dfs[0] 
    chrm_mid_df = dfs[1]
    chrm_mid_df["cnv_number"], chrm_mid_df["gene"], chrm_mid_df["category"] = zip(*chrm_mid_df.apply(lambda mid_row: get_snp_info(chrm_top_data, mid_row["pos"]), axis=1))
    return chrm_mid_df


def calc_color(ratio):
    if ratio > 0.4:
        color = "gain"
    elif ratio < -0.7:
        color = "loss"
    else:
        color = "normal"
    return color


def calc_limits(df):
    df["category"] = df.apply(lambda row: calc_color(row["log2ratio"]), axis=1) 
    return df


def main():
    col_dtypes = {
        'chrom': str, 'start': np.int32, \
        'stop': np.int32, 'pos': np.int32, \
        'vaf': np.float64, 'log2ratio': np.float64, \
        'category': str, 'gene': str, \
        'arm': str, 'transcript': str
    }

    parser = argparse.ArgumentParser()
    parser.add_argument("--num_cores", dest="num_cores", type=int, default=1)
    parser.add_argument("--top_panel", dest="top_panel", required=True)
    parser.add_argument("--mid_panel", dest="mid_panel", required=True)
    parser.add_argument("--sample_id", dest="sample_id", required=True)
    args = parser.parse_args()

    top_df = pd.read_csv(args.top_panel, sep='\t', header=0, names=["chrom", "start", "stop", "log2ratio", "type", "gene", "arm", "transcript"], dtype=col_dtypes)
    mid_df = pd.read_csv(args.mid_panel, sep='\t', header=0, names=["chrom", "pos", "vaf"], dtype=col_dtypes)

    vertical_lines = [0]
    chrom_changes_mask = top_df["chrom"].ne(top_df["chrom"].shift().bfill()).astype(int)
    chrom_changes = chrom_changes_mask.where(chrom_changes_mask == 1).dropna().index.values.tolist()
    vertical_lines.extend(chrom_changes)
    vertical_lines.append(len(top_df))

    print('Parallelizing limit calculation...')
    top_df_split = np.array_split(top_df, 10)
    pool = Pool(args.num_cores)
    new_top_df = pd.concat(pool.map(calc_limits, top_df_split))
    pool.close()
    pool.join()
    final_top_df = new_top_df[["chrom", "start", "stop", "log2ratio", "category", "gene", "arm", "transcript"]]

    chr_pairs = list()
    for chrom in mid_df["chrom"].unique():
        chr_pairs.append([final_top_df[final_top_df["chrom"] == chrom], mid_df[mid_df["chrom"] == chrom]])

    print('Parallelizing SNP lookup...')
    pool = Pool(args.num_cores)
    new_mid_df = pd.concat(pool.map(snp_map, chr_pairs))
    pool.close()
    pool.join()

    final_mid_df = new_mid_df[["cnv_number", "vaf", "pos", "gene", "category", "chrom"]]
    final_mid_df.dropna(inplace=True)

    sample_n = int(round(len(final_mid_df)/3000))
    print('Sampling every {} SNP...'.format(sample_n))
    final_mid_df = final_mid_df.iloc[::sample_n, :]

    final_top_df.to_csv(args.sample_id + '.panels.txt', sep='\t', index=False)
    final_mid_df.to_csv(args.sample_id + '.middle.panels.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
